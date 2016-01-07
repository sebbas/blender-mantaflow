#include "MANTA.h"
#include "WTURBULENCE.h"
#include "scenarios/smoke.h"

extern "C" bool manta_check_grid_size(struct FLUID_3D *fluid, int dimX, int dimY, int dimZ)
{
	/*Y and Z axes are swapped in manta and blender*/
	if (!(dimX == fluid->xRes() && dimY == fluid->yRes() && dimZ == fluid->zRes())) {
		for (int cnt(0); cnt < fluid->_totalCells; cnt++)
			fluid->_density[cnt] = 0.0f;
		return false;
	}
	return true;
}

extern "C" bool manta_check_wavelets_size(struct WTURBULENCE *wt, int dimX, int dimY, int dimZ)
{
	/*Y and Z axes are swapped in manta and blender*/
	if (!(dimX == wt->_xResBig && dimY == wt->_yResBig && dimZ == wt->_zResBig)) {
		for (int cnt(0); cnt < wt->_totalCellsBig; cnt++)
			wt->_densityBig[cnt] = 0.0f;
		return false;
	}
	return true;
}

void read_rotated_grid(gzFile gzf, float *data, int size_x, int size_y, int size_z)
{
	assert(size_x > 1 && size_y > 1 && size_z > 1);
	float* temp_data = (float*)malloc(sizeof(float) * size_x * size_y * size_z);
//	data = (float*)malloc(sizeof(float) * size_x * size_y * size_z);
	gzread(gzf, temp_data, sizeof(float)* size_x * size_y * size_z);
	for (int cnt_x(0); cnt_x < size_x; ++cnt_x)
{
		for (int cnt_y(0); cnt_y < size_y; ++cnt_y)
		{
			for (int cnt_z(0); cnt_z < size_z; ++cnt_z)
			{
				data[cnt_x + size_x * cnt_y + size_x*size_y * cnt_z] = temp_data[cnt_x + size_x * cnt_y + size_x*size_y * cnt_z];
			}			
		}
	}
}

static void wavelets_add_lowres_density(SmokeDomainSettings *sds)
{
	assert(sds != NULL);
	for (int cnt_x(0); cnt_x < sds->wt->_xResBig; ++cnt_x)
	{
		for (int cnt_y(0); cnt_y < sds->wt->_yResBig; ++cnt_y)
		{
			for (int cnt_z(0); cnt_z < sds->wt->_zResBig; ++cnt_z)
			{
				//scale down to domain res
				float x_sc = 1. * sds->base_res[0] * cnt_x / sds->wt->_xResBig;
				float y_sc = 1. * sds->base_res[1] * cnt_y / sds->wt->_yResBig;
				float z_sc = 1. * sds->base_res[2] * cnt_z / sds->wt->_zResBig;
				//finding cells to interpolate from
				int start_x = int(x_sc / 1);
				int start_y = int(y_sc / 1);
				int start_z = int(z_sc / 1);
				int end_x = ((x_sc - start_x > 0.001) && (start_x + 1 < sds->base_res[0]))? start_x + 1: start_x;
				int end_y = ((y_sc - start_y > 0.001) && (start_y + 1 < sds->base_res[1]))? start_y + 1: start_y;
				int end_z = ((z_sc - start_z > 0.001) && (start_z + 1 < sds->base_res[2]))? start_z + 1: start_z;
				//interpolation
				float add_value = 0;
				int cnt=0;
				for(int x(start_x); x <= end_x; ++x)
				{
					for(int y(start_y); y <= end_y; ++y)
					{
						for(int z(start_z); z <= end_z; ++z)
						{
							cnt++;
							add_value += sds->fluid->_density[x + y*sds->base_res[0] + z * sds->base_res[0]*sds->base_res[1]];	
						}			
					}	
				}
				add_value /= float(cnt);
				sds->wt->_densityBig[cnt_x + cnt_y *sds->wt->_xResBig + cnt_z*sds->wt->_xResBig*sds->wt->_yResBig] += add_value;
			}
		}
	}
}

//PR need SMD data here for wavelets 
extern "C" int read_mantaflow_sim(struct SmokeDomainSettings *sds, char *name, bool reading_wavelets)
{
    /*! l /*! legacy headers for reading old files */
	typedef struct {
		int dimX, dimY, dimZ;
		int frames, elements, elementType, bytesPerElement, bytesPerFrame;
	} UniLegacyHeader;
	
	typedef struct {
		int dimX, dimY, dimZ;
		int gridType, elementType, bytesPerElement;
	} UniLegacyHeader2;
	
	/* uni file header - currently used */ 
	typedef struct {
		int dimX, dimY, dimZ;
		int gridType, elementType, bytesPerElement;
		char info[256]; /* mantaflow build information */
		unsigned long long timestamp; /* creation time */
	} UniHeader;
	
#	if NO_ZLIB!=1
    gzFile gzf = gzopen(name, "rb");
//    if (!gzf) {
//		if(reading_wavelets){
//			for (int cnt(0); cnt < sds->wt->_totalCellsBig; cnt++)
//				sds->wt->_densityBig[cnt] = 0.0f;
//		}
//		else{
//			for (int cnt(0); cnt < sds->fluid->_totalCells; cnt++)
//				sds->fluid->_density[cnt] = 0.0f;
//		}
//		return 0;
//	}
	
    char ID[5] = {0,0,0,0,0};
	gzread(gzf, ID, 4);
	/* legacy file format */
    if (!strcmp(ID, "DDF2")) {
        UniLegacyHeader head;
		gzread(gzf, &head, sizeof(UniLegacyHeader));
        int numEl = head.dimX*head.dimY*head.dimZ;
        gzseek(gzf, numEl, SEEK_CUR);
        /* actual grid read */
        if ( ! reading_wavelets){
//			if (!manta_check_grid_size(sds->fluid, head.dimX, head.dimY, head.dimZ))	return 0;
			gzread(gzf, sds->fluid->_density, sizeof(float)*numEl);
		}
		else {
			if (!manta_check_wavelets_size(sds->wt, head.dimX, head.dimY, head.dimZ))	return 0;
			gzread(gzf, sds->wt->_densityBig, sizeof(float)*numEl);
    	} 
	}
	/* legacy file format 2 */
    else if (!strcmp(ID, "MNT1")) {
        UniLegacyHeader2 head;
        gzread(gzf, &head, sizeof(UniLegacyHeader2));
		/* actual grid read*/
        if ( ! reading_wavelets){
//			if (!manta_check_grid_size(sds->fluid, head.dimX, head.dimY, head.dimZ))	return 0;
        	gzread(gzf, sds->fluid->_density, sizeof(float)*head.dimX*head.dimY*head.dimZ);
    	}
		else{
			if (!manta_check_wavelets_size(sds->wt, head.dimX, head.dimY, head.dimZ))	return 0;
        	gzread(gzf, sds->wt->_densityBig, sizeof(float)*head.dimX*head.dimY*head.dimZ);
    	}
	}
	/* current file format*/
    else if (!strcmp(ID, "MNT2")) {
        UniHeader head;
        gzread(gzf, &head, sizeof(UniHeader));
		/* actual grid read */
        if ( ! reading_wavelets){
//			if (!manta_check_grid_size(sds->fluid, head.dimX, head.dimY, head.dimZ))	return 0;
			/*Y and Z axes are swapped in manta and blender*/
			gzread(gzf,sds->fluid->_density, sizeof(float)*head.dimX*head.dimY*head.dimZ);
    		
		}
		else{
			if (!manta_check_wavelets_size(sds->wt, head.dimX, head.dimY, head.dimZ))	return 0;
			/*Y and Z axes are swapped in manta and blender*/
			gzread(gzf,sds->wt->_densityBig, sizeof(float)*head.dimX*head.dimY*head.dimZ);
			gzread(gzf,sds->wt->_densityBigOld, sizeof(float)*head.dimX*head.dimY*head.dimZ);
//			wavelets_add_lowres_density(sds);
		}
	}
    gzclose(gzf);
	return 1;
#	endif	/*zlib*/
	return 0;
}


void Manta_API::indent_ss(stringstream& ss, int indent)
{
	/*two-spaces indent*/
	if (indent < 0) return;
	std::string indentation = ""; 
	for (size_t cnt(0); cnt < indent; ++cnt) {
		indentation += "  ";
	}
	ss << indentation;
}

void Manta_API::manta_gen_noise(stringstream& ss, char* solver, int indent, char *noise, int seed, bool load, bool clamp, float clampNeg, float clampPos, float valScale, float valOffset, float timeAnim)
{
	if (ss == NULL)/*should never be here*/
	{
		return;
	}
	indent_ss(ss, indent);
	ss << noise << " = "<<solver<<".create(NoiseField, fixedSeed=" << seed << ", loadFromFile="<< (load?"True":"False") <<") \n";
	ss << noise << ".posScale = vec3(20) \n";
	ss << noise << ".clamp = " << ((clamp)?"True":"False") << " \n";
	ss << noise << ".clampNeg = " << clampNeg << " \n";
	ss << noise << ".clampPos = " << clampPos << " \n";
	ss << noise << ".valScale = " << valScale << " \n";
	ss << noise << ".valOffset = " << valOffset << " \n";
	ss << noise << ".timeAnim = " << timeAnim << " \n";
}

void Manta_API::manta_solve_pressure(stringstream& ss, char *flags, char *vel, char *pressure, bool useResNorms, int openBound, int solver_res,float cgMaxIterFac, float cgAccuracy)
{
	/*open:0 ; vertical : 1; closed:2*/
	ss << "  solvePressure(flags=" << flags << ", vel=" << vel << ", pressure=" << pressure << ", useResNorm=" << (useResNorms?"True":"False") << ", openBound='";	
	
	if(openBound == 1) /*vertical*/
	{
		ss << "yY'";
	}
	else if (openBound == 0) /*open*/
	{
		if(solver_res == 2)
			ss << "xXyY";
		else
			ss << "xXyYzZ";
	}
	ss << "'";	/*empty for closed bounds*/ 
	
	ss << ", cgMaxIterFac=" << cgMaxIterFac << ", cgAccuracy=" << cgAccuracy << ") \n";
}

void Manta_API::manta_advect_SemiLagr(stringstream& ss, int indent, char *flags, char *vel, char *grid, int order)
{
	if((order <=1) || (flags == NULL) || (vel == NULL) || (grid == NULL)){return;}
	indent_ss(ss, indent);
	ss << "advectSemiLagrange(flags=" << flags << ", vel=" << vel \
	<< ", grid=" << grid << ", order=" << order << ") \n"; 
}

/*create solver, handle 2D case*/
void Manta_API::manta_create_solver(stringstream& ss, char *name, char *nick, char *grid_size_name, int x_res, int y_res, int z_res, int dim)
{
	if ((dim != 2) && (dim != 3))
	{ return; }
	if (dim == 2)
	{ y_res = 1; }
	ss << grid_size_name << " = vec3(" << x_res << ", " << y_res << ", " << z_res << ")" << " \n";
	ss << name << " = Solver(name = '" << nick << "', gridSize = " << grid_size_name << ", dim = " << dim << ") \n";
}

inline bool file_exists (const std::string& name) {
    return ( access( name.c_str(), F_OK ) != -1 );
}

/*blender transforms obj coords to [-1,1]. This method transforms them back*/
void Manta_API::add_mesh_transform_method(stringstream& ss)
{
	ss << "def transform_back(obj, gs):\n" <<
	"  obj.scale(gs/2)\n" <<
	"  obj.offset(gs/2)\n\n";
}

void Manta_API::stop_manta_sim()
{
	pthread_cancel(manta_thread);
}

string Manta_API::gridNameFromType(const string &type)
{
	if (type == "float")
	{
		return "RealGrid";
	}
	else if (type == "Vec3")
	{
		return "MACGrid";	
	}
	else
	{
		cout<<"ERROR: can not create grid from type: "<< type << endl;
		return "";
	}
}

void Manta_API::addGrid(void * data, string name, string type, int x, int y, int z, bool is2D = false)
{
	if (data == NULL || name == "" || gridNameFromType(type) == "") return;
	cout << "Adding Grid: " << name << endl;
	std::ostringstream stringStream;
	
	/* Temporary gridname */
	stringStream << "temp_" << name;
	std::string grid_name = stringStream.str();
	
	/* Generate command that sets up grid variable */
	stringStream.str("");
	stringStream << grid_name << " = s.create(" << gridNameFromType(type) << ")";
	const std::string command_1 = stringStream.str();

	/* Generate command that makes the grid read in the data */
	stringStream.str("");
	if (is2D)
	{
		/* For 2D case, Y and Z axes are switched, Y axis is '1' for Mantaflow */
		stringStream << grid_name << ".readGridFromMemory(\'"<< data << "\', " << x << "," << z << "," << 1 << ")";
	}
	else
	{
		stringStream << grid_name << ".readGridFromMemory(\'"<< data << "\', " << x << "," << y << "," << z << ")";
	}
	const std::string command_2 = stringStream.str();
	
	/* Generate command that maps temp grid to our 'real' grid */
	const std::string command_3 = name + ".add(" + grid_name + ")";
	
	/* Execute all commands */
	PyGILState_STATE gilstate = PyGILState_Ensure();
	PyRun_SimpleString(command_1.c_str());
	PyRun_SimpleString(command_2.c_str());
	PyRun_SimpleString(command_3.c_str());
	PyGILState_Release(gilstate);		
}

void Manta_API::addAdaptiveGrid(void * data, string gridName, string solverName, string type, int minX, int minY, int minZ, int maxX, int maxY, int maxZ, bool is2D = false)
{
	if (data == NULL || gridName == "" || gridNameFromType(type) == "") return;
	{
		cout << "NULL values passed to grid addAdaptiveGrid for grid " << gridName <<endl;
		return;
	}
	std::ostringstream stringStream;
	stringStream << "temp_" <<gridName;
	std::string temp_grid_name = stringStream.str();
	stringStream.str("");
	stringStream << temp_grid_name << " = "<< solverName << ".create(" << gridNameFromType(type) << ")";
	const std::string command_1 = stringStream.str();
	stringStream.str("");
	
	if (is2D){
		stringStream << temp_grid_name << ".readAdaptiveGridFromMemory(\'"<< data << "\', vec3(" << minX << "," << minZ << "," << 1 << 
		"), vec3(" << maxX << "," << maxZ << "," << 1 << ") )";
	}
	else{
		stringStream << temp_grid_name << ".readAdaptiveGridFromMemory(\'"<< data << "\', vec3(" << minX << "," << minY << "," << minZ << 
		"), vec3(" << maxX << "," << maxY << "," << maxZ << ") )";
	}
	const std::string command_2 = stringStream.str();
	const std::string command_3 = gridName + ".add(" + temp_grid_name + ")";
	PyGILState_STATE gilstate = PyGILState_Ensure();
	PyRun_SimpleString("print('Reading Adaptive grid from memory')");
	PyRun_SimpleString("print (s)");
	PyRun_SimpleString(command_1.c_str());
	PyRun_SimpleString(command_2.c_str());
	PyRun_SimpleString(command_3.c_str());
	PyGILState_Release(gilstate);		
}

void Manta_API::export_obstacles(float *data, int x, int y, int z, bool is2D = false)
{
	if (data == NULL){
		cout << "NULL passed to grid export_obstacles " <<endl;  return;
	}
	std::ostringstream stringStream;
	std::string grid_name = "obs_sdf";
	stringStream.str("");
	stringStream << grid_name << " = s.create(RealGrid)";
	const std::string command_1 = stringStream.str();
	stringStream.str("");
	cout<<"Exporting obstacles"<<endl;
	if (is2D){
		stringStream << grid_name << ".readGridFromMemory(\'"<< data << "\', " << x << "," << z << "," << 1 << ")";
	}
	else{
		stringStream << grid_name << ".readGridFromMemory(\'"<< data << "\', " << x << "," << y << "," << z << ")";
	}
	const std::string command_2 = stringStream.str();
	const std::string command_3 = grid_name + ".applyToGrid(grid = flags, value = FlagObstacle)";
	PyGILState_STATE gilstate = PyGILState_Ensure();
	PyRun_SimpleString(command_1.c_str());
	PyRun_SimpleString(command_2.c_str());
	PyRun_SimpleString(command_3.c_str());
	PyGILState_Release(gilstate);		
}

std::string Manta_API::get_manta_smoke_script(SmokeModifierData *smd)
{
    std::string smoke_script = "";
	
	// Check if high res is enabled
	// Need to check if wt exists and NOT just check if high res flag is set (smd->domain->flags & MOD_SMOKE_HIGHRES)
	// because wt might not exist, i.e. when FLUID_3D constructor is called before WTURBULENCE constructor
	if (smd->domain->wt) {
		smoke_script = smoke_setup_high + smoke_import_high + smoke_step_high;
	} else {
		// TODO: Need to figure out how to handle liquids when high resolution grids are enabled, not just for low res grids
		if (smd->domain->flags & MOD_SMOKE_MANTA_USE_LIQUID)
			smoke_script = smoke_setup_low  + liquid_step_low;
		else
			smoke_script = smoke_setup_low + smoke_import_low + smoke_inflow_low + smoke_step_low ;
	}
	return smoke_script;
}

void Manta_API::run_manta_sim_file_lowRes(SmokeModifierData *smd)
{
	std::string smoke_script = get_manta_smoke_script(smd);
	std::string final_script = parseScript(smoke_script, smd);
	
	PyGILState_STATE gilstate = PyGILState_Ensure();
	PyRun_SimpleString(final_script.c_str());
	PyGILState_Release(gilstate);
}

void Manta_API::run_manta_sim_file_highRes(SmokeModifierData *smd)
{
	std::string smoke_script = get_manta_smoke_script(smd);
	std::string final_script = parseScript(smoke_script, smd);
	
	PyGILState_STATE gilstate = PyGILState_Ensure();
	PyRun_SimpleString(final_script.c_str());
	PyGILState_Release(gilstate);
}

std::string Manta_API::getRealValue( const std::string& varName, SmokeModifierData *smd)
{
	ostringstream ss;
	cout << "name is " << varName << endl;
	bool is2D = smd->domain->fluid->manta_resoution == 2;
	if (varName == "UVS_CNT")
		ss << smd->domain->manta_uvs_num ;
	else if (varName == "UPRES")
		ss << smd->domain->amplify+1;
	else if (varName == "WLT_STR")
		ss << smd->domain->strength ;
	else if (varName == "RES")
		ss <<  smd->domain->maxres;
	else if (varName == "RESX")
		ss <<  smd->domain->fluid->_xRes;
	else if (varName == "RESY")
		if (is2D) {	ss <<  smd->domain->fluid->_zRes;}
		else { 		ss <<  smd->domain->fluid->_yRes;}
	else if (varName == "RESZ")
		if (is2D){	ss << 1;}
		else { 		ss << smd->domain->fluid->_zRes;}
	else if (varName == "SOLVER_DIM")
		ss <<  smd->domain->manta_solver_res;
	else if (varName == "NOISE_CN")
		ss <<  smd->domain->noise_clamp_neg;
	else if (varName == "NOISE_CP")
		ss <<  smd->domain->noise_clamp_pos;
	else if (varName == "NOISE_VALSCALE")
		ss <<  smd->domain->noise_val_scale;
	else if (varName == "NOISE_VALOFFSET")
		ss << smd->domain->noise_val_offset;
	else if (varName == "NOISE_TIMEANIM")
		ss << smd->domain->noise_time_anim;
	else if (varName == "HRESX")
		ss << smd->domain->wt->getResBig()[0];
	else if (varName == "HRESY")
		if (is2D) {	ss << smd->domain->wt->getResBig()[2];}
		else { 		ss << smd->domain->wt->getResBig()[1];}
	else if (varName == "HRESZ")
		if (is2D) {	ss << 1;}
		else { 		ss << smd->domain->wt->getResBig()[2];}
	else if (varName == "TIMESTEP")
		ss << smd->domain->time_scale * 0.1f;
	else if (varName == "XL_TIMESTEP")
		ss << smd->domain->time_scale * 0.1f;
	else if (varName == "USE_WAVELETS")
		ss << ((smd->domain->flags & MOD_SMOKE_HIGHRES) ? "True" : "False");
	else if (varName == "BUYO_X")
		ss << 0.;
	else if (varName == "BUYO_Y")
		ss << 0.;
	else if (varName == "BUYO_Z")
		ss << (-smd->domain->beta);
	else if (varName == "ALPHA")
		ss << (-smd->domain->alpha);
	else if (varName == "BETA")
		ss << (-smd->domain->beta);
	else if (varName == "ADVECT_ORDER")
		ss << 2;
	else if (varName == "MANTA_EXPORT_PATH") {
		char parent_dir[1024];
		BLI_split_dir_part(smd->domain->_manta_filepath, parent_dir, sizeof(parent_dir));
		ss << parent_dir;
	} else if (varName == "VORTICITY") {
		ss << smd->domain->vorticity / smd->domain->fluid->_constantScaling;
	} else if (varName == "BOUNDCONDITIONS") {
		// OLD SETUP. WHY LIKE THAT??
		/*if(smd->domain->border_collisions == SM_BORDER_OPEN) ss << "xXyY";
		else if (smd->domain->border_collisions == SM_BORDER_VERTICAL) ss << "xXyY";
		else if (smd->domain->border_collisions == SM_BORDER_CLOSED) ss << "xXyY";
		
		if (smd->domain->manta_solver_res == 3){
			if(smd->domain->border_collisions == SM_BORDER_OPEN) ss << "z";
			else if (smd->domain->border_collisions == SM_BORDER_VERTICAL) ss << "z";
			else if (smd->domain->border_collisions == SM_BORDER_CLOSED) ss << "zZ";
		}*/
		if(smd->domain->border_collisions == SM_BORDER_OPEN) ss << "xXyY";
		else if (smd->domain->border_collisions == SM_BORDER_VERTICAL) ss << "zZ";
		else if (smd->domain->border_collisions == SM_BORDER_CLOSED) ss << "";
		
		if (smd->domain->manta_solver_res == 3) {
			if(smd->domain->border_collisions == SM_BORDER_OPEN) ss << "zZ";
			else if (smd->domain->border_collisions == SM_BORDER_VERTICAL) ss << "";
			else if (smd->domain->border_collisions == SM_BORDER_CLOSED) ss << "";
		}
	}
	else if (varName == "DO_OPEN")
		ss << ((smd->domain->border_collisions == SM_BORDER_CLOSED) ? "False" : "True");
	else if (varName == "GRAVITY")
		ss << "vec3(0,0,-0.981)";
	else if (varName == "ABS_FLOW")
		ss << ((smd->flow->flags & MOD_SMOKE_FLOW_ABSOLUTE)? "True" : "False");
	else if (varName == "DENSITY_MEM")
		ss << smd->domain->fluid->_density;
	else if (varName == "DENSITY_SIZE")
		ss << sizeof(float) * smd->domain->total_cells;
	else if (varName == "XL_DENSITY_MEM")
		ss << smd->domain->wt->_densityBig;
	else if (varName == "XL_DENSITY_SIZE")
		ss << sizeof(float) * smd->domain->wt->_xResBig * smd->domain->wt->_yResBig * smd->domain->wt->_zResBig;
	else if (varName == "BURNING_RATE")
		ss << (smd->domain->burning_rate);
	else if (varName == "FLAME_SMOKE")
		ss << (smd->domain->flame_smoke);
	else if (varName == "IGNITION_TEMP")
		ss << (smd->domain->flame_ignition);
	else if (varName == "MAX_TEMP")
		ss << (smd->domain->flame_max_temp);
	else if (varName == "DT")
		ss << smd->domain->fluid->_dt;
	else if (varName == "FLAME_SMOKE_COLOR_X")
		ss << smd->domain->flame_smoke_color[0];
	else if (varName == "FLAME_SMOKE_COLOR_Y")
		ss << smd->domain->flame_smoke_color[1];
	else if (varName == "FLAME_SMOKE_COLOR_Z")
		ss << smd->domain->flame_smoke_color[2];
	else if (varName == "USING_COLORS")
		ss << ((smd->domain->fluid->using_colors) ? "True" : "False");
	else if (varName == "USING_HEAT")
		ss << ((smd->domain->fluid->using_heat) ? "True" : "False");
	else if (varName == "USING_FIRE")
		ss << ((smd->domain->fluid->using_fire) ? "True" : "False");
	else 
		cout << "ERROR: Unknown option:" << varName <<endl;
	return ss.str();
}

std::string Manta_API::parseLine(const string& line, SmokeModifierData *smd)
{
	if (line.size() == 0) return "";
	string res = "";
	int currPos = 0, start_del = 0, end_del = -1;
	bool readingVar = false;
	const char delimiter = '$';
	while (currPos < line.size()){
		if(line[currPos] == delimiter && ! readingVar) {
			readingVar 	= true;
			start_del	= currPos + 1;
			res 		+= line.substr(end_del + 1, currPos - end_del -1);
		}
		else if(line[currPos] == delimiter && readingVar){
			readingVar 	= false;
			end_del 	= currPos;
			res 		+= getRealValue(line.substr(start_del, currPos - start_del), smd);
		}
		currPos ++;
	}
	res += line.substr(end_del+1, line.size()- end_del);
	return res;
}

std::string Manta_API::parseScript(const string& setup_string, SmokeModifierData *smd)
{
	std::istringstream f(setup_string);
	ostringstream res;
	string line = "";
	while(getline(f, line)) {
		res << parseLine(line, smd) << "\n";
	}
	return res.str();
}

void Manta_API::manta_export_grids(SmokeModifierData *smd)
{
	// Export the scene file
	std::string smoke_script = get_manta_smoke_script(smd);
	
	std::string final_script = "";
	if (smd->domain->flags & MOD_SMOKE_HIGHRES) {
		final_script = Manta_API::parseScript(smoke_script, smd) + standalone_high;
	} else {
		final_script = Manta_API::parseScript(smoke_script, smd) + standalone_low;
	}
	
	ofstream myfile;
	myfile.open(smd->domain->_manta_filepath);
	myfile << final_script;
	myfile.close();
	
	// Run python environment to export grids, that is, create the grid files
	PyGILState_STATE gilstate = PyGILState_Ensure();
	if (smd->domain->flags & MOD_SMOKE_HIGHRES) {
		PyRun_SimpleString(Manta_API::parseScript(smoke_export_high, smd).c_str());
	} else {
		PyRun_SimpleString(Manta_API::parseScript(smoke_export_low, smd).c_str());
	}
	PyGILState_Release(gilstate);
}

string Manta_API::getGridPointer(std::string gridName, std::string solverName)
{
	if ((gridName == "") && (solverName == "")) {
		return "";
	}
	cout << "getting grid pointer " << gridName<< " , " << solverName <<endl;
	PyGILState_STATE gilstate = PyGILState_Ensure();
	PyObject *main = PyImport_AddModule("__main__");
	if (main == NULL){cout << "null" << 1 << endl;return "";}
	PyObject *globals = PyModule_GetDict(main);
	if (globals == NULL){cout << "null" << 12 << endl;return "";}
	PyObject *grid_object = PyDict_GetItemString(globals, gridName.c_str());
	if (grid_object == NULL){cout << "null" << 13 << endl;return "";}
	PyObject* func = PyObject_GetAttrString(grid_object,(char*)"getDataPointer");
	if (func == NULL){cout << "null" << 14 << endl;return "";}
	PyObject* retured_value = PyObject_CallObject(func, NULL);
	PyObject* encoded = PyUnicode_AsUTF8String(retured_value);
	if (retured_value == NULL){cout << "null" << 15 << endl;return "";}
	std::string res = PyBytes_AsString(encoded);
	cout << "Pointer on "<< gridName << " " << res << endl;
	PyGILState_Release(gilstate);		
	return res;
}

// init direct access functions from blender
void Manta_API::initBlenderRNA(float *alpha, float *beta, float *dt_factor, float *vorticity, int *borderCollision, float *burning_rate,
							  float *flame_smoke, float *flame_smoke_color, float *flame_vorticity, float *flame_ignition_temp, float *flame_max_temp)
{
	_alpha = alpha;
	_beta = beta;
	_dtFactor = dt_factor;
	_vorticityRNA = vorticity;
	_borderColli = borderCollision;
	_burning_rate = burning_rate;
	_flame_smoke = flame_smoke;
	_flame_smoke_color = flame_smoke_color;
	_flame_vorticity = flame_vorticity;
	_ignition_temp = flame_ignition_temp;
	_max_temp = flame_max_temp;
}

void * Manta_API::pointerFromString(const std::string& s){
	stringstream ss(s);
	void *gridPointer = NULL;
	ss >> gridPointer;
	return gridPointer;
}


void Manta_API::updatePointers(FLUID_3D *fluid)
{
	cout << "Updating pointers" << endl;
	if (fluid->manta_resoution == 2)
	{
		cout << "2D" << endl;
		float* manta_fluid_density = (float* )pointerFromString(getGridPointer("density", "s")); 
		int* manta_fluid_flags = (int* )pointerFromString(getGridPointer("flags", "s"));
		if (fluid->_density != NULL){
			for (int cnt(0); cnt < fluid->xRes() * fluid->yRes() * fluid->zRes(); ++cnt) {
				fluid->_density[cnt] = 0.;
				fluid->_manta_flags[cnt] = 2;
			}
		}
		int step = 0;
		for (int cnty(0);cnty<fluid->yRes(); ++cnty)
			for(int cntz(0);cntz<fluid->zRes(); ++cntz)
			{
				step = fluid->xRes() + cnty * fluid->xRes() + cntz * fluid->xRes()*fluid->yRes(); 
				if ((step < 0) || (step > fluid->_totalCells)){
					cout << "UpdatePointers: step is larger than cell dim" << step << endl;
				}
				fluid->_density[step] = manta_fluid_density[cnty + cntz*fluid->xRes()];
				fluid->_manta_flags[step] = manta_fluid_flags[cnty + cntz*fluid->xRes()];
			}		
	}
	else {
		cout << "3D" << endl;
		fluid->_density = (float* )pointerFromString(getGridPointer("density", "s"));	
		fluid->_manta_flags = (int* )pointerFromString(getGridPointer("flags", "s"));
	}
	
	fluid->_manta_inflow = (float* )pointerFromString(getGridPointer("inflow_grid", "s"));
	fluid->_fuel_inflow = (float* )pointerFromString(getGridPointer("fuel_inflow", "s"));
	if (fluid-> manta_resoution == 2) {return;}
	if (fluid->using_colors) {
		fluid->_color_r = (float* )pointerFromString(getGridPointer("color_r_low", "s"));
		fluid->_color_g = (float* )pointerFromString(getGridPointer("color_g_low", "s"));
		fluid->_color_b = (float* )pointerFromString(getGridPointer("color_b_low", "s"));
	}
	if (fluid->using_heat) {
		fluid->_heat = (float* )pointerFromString(getGridPointer("heat_low", "s"));
	}
	if (fluid->using_fire) {
		fluid->_flame = (float* )pointerFromString(getGridPointer("flame_low", "s"));
		fluid->_fuel = (float* )pointerFromString(getGridPointer("fuel_low", "s"));
		fluid->_react = (float* )pointerFromString(getGridPointer("react_low", "s"));
	}
	fluid->_xVelocity = (float* )pointerFromString(getGridPointer("x_vel", "s"));
	fluid->_yVelocity = (float* )pointerFromString(getGridPointer("y_vel", "s"));
	fluid->_zVelocity = (float* )pointerFromString(getGridPointer("z_vel", "s"));
}

void Manta_API::updateHighResPointers(WTURBULENCE *wt)
{
	cout << "Updating pointers high res" << endl;
	wt->_densityBig = (float* )pointerFromString(getGridPointer("xl_density", "xl"));;
	if (wt->using_colors) {
		wt->_color_rBig = (float* )pointerFromString(getGridPointer("color_r_high", "xl"));
		wt->_color_gBig = (float* )pointerFromString(getGridPointer("color_g_high", "xl"));
		wt->_color_bBig = (float* )pointerFromString(getGridPointer("color_b_high", "xl"));
	}
	if (wt->using_fire) {
		wt->_flameBig = (float* )pointerFromString(getGridPointer("flame_high", "xl"));
		wt->_fuelBig = (float* )pointerFromString(getGridPointer("fuel_high", "xl"));
		wt->_reactBig = (float* )pointerFromString(getGridPointer("react_high", "xl"));
	}
}

