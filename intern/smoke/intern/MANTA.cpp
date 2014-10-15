#include "MANTA.h"
#include "WTURBULENCE.h"
#include "scenarios/smoke.h"

Manta_API* Manta_API::_instance = 0;
Manta_API* Manta_API::instance(){
	if (_instance == 0){
		_instance = new Manta_API;
	}
	return _instance;
}

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
    if (!gzf) {
		if(reading_wavelets){
			for (int cnt(0); cnt < sds->wt->_totalCellsBig; cnt++)
				sds->wt->_densityBig[cnt] = 0.0f;
		}
		else{
			for (int cnt(0); cnt < sds->fluid->_totalCells; cnt++)
				sds->fluid->_density[cnt] = 0.0f;
		}
		return 0;
	}
	
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
			if (!manta_check_grid_size(sds->fluid, head.dimX, head.dimY, head.dimZ))	return 0;
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
			if (!manta_check_grid_size(sds->fluid, head.dimX, head.dimY, head.dimZ))	return 0;
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
			if (!manta_check_grid_size(sds->fluid, head.dimX, head.dimY, head.dimZ))	return 0;
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

void Manta_API::manta_cache_path(char *filepath)
{
	char *name="manta";
	BLI_make_file_string("/", filepath, BLI_temporary_dir(), name);
}

//void BLI_dir_create_recursive(const char *filepath);
void Manta_API::create_manta_folder()
{
	char* filepath=NULL;
	manta_cache_path(filepath);
	//BLI_dir_create_recursive(filepath);
	
}

void *Manta_API::run_manta_scene_thread(void *arguments)
{
//	struct manta_arg_struct *args = (struct manta_arg_struct *)arguments;
//	//create_manta_folder();
//	//PyInterpreterState *st = PyThreadState_GET()->interp;
//	//PyThreadState *ts = Py_NewInterpreter();
//	
//	vector<string> a;
//	a.push_back(args->filepath);
//	//a.push_back("manta_scene.py");
//	//args.push_back("test_1.py");
//	
//	runMantaScript(a);
//	
//	//system("./manta manta_scene.py");
//	pthread_exit(NULL);
	return NULL;
}

void Manta_API::run_manta_scene(Scene *s, SmokeModifierData *smd)
{
//	smd->domain->manta_sim_frame = 0;
//	PyGILState_STATE gilstate = PyGILState_Ensure();
////	for (int fr=0; fr< 1; ++fr)
//	int fr = s->r.cfra;
//	{
////		if(smd->domain->manta_sim_frame == -1)
////			break;
//		printf("Simulation Step");
//		manta_write_effectors(s, smd);
//		smd->domain->manta_sim_frame = fr;
//		std::string frame_str = static_cast<ostringstream*>( &(ostringstream() << fr) )->str();
//		std::string py_string_0 = string("sim_step(").append(frame_str);
//		std::string py_string_1 = py_string_0.append(")\0");
//		//		std::string py_string_1 = string("sim_step()\0");
//		PyRun_SimpleString(py_string_1.c_str());
//		//		frame_num ++;
//	}
//	PyGILState_Release(gilstate);
	//returning simulation state to "not simulating" aka -1
//	smd->domain->manta_sim_frame = -1;
//
//	
//	
	struct manta_arg_struct *args = (struct manta_arg_struct*)malloc(sizeof(struct manta_arg_struct));
	args->smd = *smd;
	args->s = *s;
//	args.frame_num = smd->domain->manta_end_frame - smd->domain->manta_start_frame;
//	int rc = pthread_create(&manta_thread, NULL, run_manta_sim_thread, (void *)args);
//	pthread_join(manta_thread,NULL);
//	pthread_detach(manta_thread);
	run_manta_sim_thread((void*) args);
}

void Manta_API::stop_manta_sim()
{
	pthread_cancel(manta_thread);
}

void Manta_API::addGrid(float * data, string name, int x, int y, int z)
{
	std::ostringstream stringStream;
	stringStream << "temp_" << name;
	std::string grid_name = stringStream.str();
	stringStream.str("");
	stringStream << grid_name << " = s.create(RealGrid)";
	const std::string command_1 = stringStream.str();
	stringStream.str("");
	stringStream << grid_name << ".readGridFromMemory("<< data << "," << x << "," << z << "," << y << ")";
	const std::string command_2 = stringStream.str();
	const std::string command_3 = name + ".add(" + grid_name + ")";
	PyGILState_STATE gilstate = PyGILState_Ensure();
	PyRun_SimpleString(command_1.c_str());
	PyRun_SimpleString(command_2.c_str());
	PyRun_SimpleString(command_3.c_str());
	PyGILState_Release(gilstate);		
}

void Manta_API::addAdaptiveGrid(float * data, string name, int minX, int minY, int minZ, int maxX, int maxY, int maxZ)
{
	if (data == NULL)
	{
		cout << "NULL pointer passed to grid addAdaptiveGrid for grid " << name <<endl;
		return;
	}
	std::ostringstream stringStream;
	stringStream << "temp_" << name;
	std::string grid_name = stringStream.str();
	stringStream.str("");
	stringStream << grid_name << " = s.create(RealGrid)";
	const std::string command_1 = stringStream.str();
	stringStream.str("");
	stringStream << grid_name << ".readAdaptiveGridFromMemory(\'"<< data << "\',\'" << name << "\', vec3(" << minX << "," << minY << "," << minZ << 
	"), vec3(" << maxX << "," << maxY << "," << maxZ << ") )";
	const std::string command_2 = stringStream.str();
	const std::string command_3 = name + ".add(" + grid_name + ")";
	PyGILState_STATE gilstate = PyGILState_Ensure();
	PyRun_SimpleString("print('Reading Adaptive grid from memory')");
	PyRun_SimpleString("print (s)");
	PyRun_SimpleString(command_1.c_str());
	PyRun_SimpleString(command_2.c_str());
	PyRun_SimpleString(command_3.c_str());
	PyGILState_Release(gilstate);		
}

void Manta_API::run_manta_sim_thread(void *arguments)
{
	struct manta_arg_struct *args = (struct manta_arg_struct *)arguments;
	SmokeModifierData *smd = &args->smd;
	Scene *s = &args->s;
	int num_sim_steps = smd->domain->manta_end_frame - smd->domain->manta_start_frame + 1;
	smd->domain->manta_sim_frame = 0;
	PyGILState_STATE gilstate = PyGILState_Ensure();
//	for (int fr=0; fr< num_sim_steps; ++fr) {
//		if(smd->domain->manta_sim_frame == -1)
//			break;
		printf("Simulation Step");
		manta_write_effectors(s, smd);
		smd->domain->manta_sim_frame = s->r.cfra;
		std::string frame_str = static_cast<ostringstream*>( &(ostringstream() << s->r.cfra) )->str();
		std::string py_string_0 = string("sim_step(").append(frame_str);
		std::string py_string_1 = py_string_0.append(")\0");
		PyRun_SimpleString(py_string_1.c_str());
		cout<< "done"<<manta_sim_running<<endl;
	//}
	//returning simulation state to "not simulating" aka -1
	smd->domain->manta_sim_frame = -1;
	PyGILState_Release(gilstate);
}

void Manta_API::generate_manta_sim_file(SmokeModifierData *smd)
{
//	 /*create python file with 2-spaces indentation*/
	bool wavelets = smd->domain->flags & MOD_SMOKE_HIGHRES;
//	bool noise_clamp = smd->domain->flags & MOD_SMOKE_NOISE_CLAMP; 
//	float noise_clamp_neg = smd->domain->noise_clamp_neg;
//	float noise_clamp_pos = smd->domain->noise_clamp_pos;
//	float noise_val_scale = smd->domain->noise_val_scale;
//	float noise_val_offset = smd->domain->noise_val_offset;
//	float noise_time_anim = smd->domain->noise_time_anim;
//	int num_sim_frames = smd->domain->manta_end_frame - smd->domain->manta_start_frame + 1;
//	if(num_sim_frames < 1)
//		return;

	/*constrcting final setup*/
	string smoke_script = smoke_setup_low + ((wavelets)?smoke_setup_high:"") + smoke_step_low + ((wavelets)?smoke_step_high:"");	
	ofstream manta_setup_file;
	manta_setup_file.open("manta_scene.py", std::fstream::trunc);
	manta_setup_file << smoke_script ;
	manta_setup_file.close();

	parseFile(smoke_script, smd);
	vector<string> a;
	a.push_back("manta_scene.py");
	runMantaScript("",a);
}

std::string Manta_API::getRealValue( const std::string& varName, SmokeModifierData *smd)
{
	ostringstream ss;
	if (varName == "UVS_CNT")
		ss << smd->domain->manta_uvs_num ;
	else if (varName == "UPRES")
		ss << smd->domain->amplify+1;
	else if (varName == "WLT_STR")
		ss << smd->domain->strength ;
	else if (varName == "RES")
		ss <<  smd->domain->maxres;
	else if (varName == "RESX")
		ss <<  smd->domain->fluid->xRes();
	else if (varName == "RESY")
		ss <<  smd->domain->fluid->yRes();
	else if (varName == "RESZ")
		ss <<  smd->domain->fluid->zRes();
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
		ss << smd->domain->wt->getResBig()[1];
	else if (varName == "HRESZ")
		ss << smd->domain->wt->getResBig()[2];
	else if (varName == "TIMESTEP")
		ss << smd->domain->time_scale * 0.1f;
	else if (varName == "XL_TIMESTEP")
		ss << smd->domain->time_scale * 0.1f;
	else if (varName == "USE_WAVELETS")
		ss << (smd->domain->flags & MOD_SMOKE_HIGHRES)?"True":"False";
	else if (varName == "BUYO_X")
		ss << 0.;
	else if (varName == "BUYO_Y")
		ss << 0.;
	else if (varName == "BUYO_Z")
		ss << (-smd->domain->beta);
	else if (varName == "ADVECT_ORDER")
		ss << 2;
	else if (varName == "ABS_FLOW")
		ss << (smd->flow->flags & MOD_SMOKE_FLOW_ABSOLUTE)?"True":"False";
	else if (varName == "DENSITY_MEM")
		ss << smd->domain->fluid->_density;
	else if (varName == "DENSITY_SIZE")
		ss << sizeof(float) * smd->domain->total_cells;
	else 
		cout<< "ERROR: Unknown option:"<< varName <<endl; 
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
		if(line[currPos] == delimiter && ! readingVar){
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

void Manta_API::parseFile(const string & setup_string, SmokeModifierData *smd)
{
//	ifstream f (file);
std::istringstream f(setup_string);
	ofstream of;
	of.open("manta_scene.py", std::fstream::trunc);
	string line="";
//	if (f.is_open()){
		while(getline(f,line)){
			of << parseLine(line,smd) << "\n";
		}
//		f.close();
//	}
//	else{
//		printf ("Error: No scenario file found");
//	}
	of.close();
}


