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

void Manta_API::manta_cache_path(char *filepath)
{
	char *name="manta";
	BLI_make_file_string("/", filepath, BLI_temp_dir_session(), name);
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

void Manta_API::run_manta_scene(Manta_API * fluid)
{
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

void Manta_API::addGrid(void * data, string name, string type, int x, int y, int z)
{
	std::ostringstream stringStream;
	stringStream << "temp_" << name;
	std::string grid_name = stringStream.str();
	stringStream.str("");
	stringStream << grid_name << " = s.create(" << gridNameFromType(type) << ")";
	const std::string command_1 = stringStream.str();
	stringStream.str("");
	stringStream << grid_name << ".readGridFromMemory(\'"<< data << "\', " << x << "," << y << "," << z << ")";
	const std::string command_2 = stringStream.str();
	const std::string command_3 = name + ".add(" + grid_name + ")";
	PyGILState_STATE gilstate = PyGILState_Ensure();
	PyRun_SimpleString(command_1.c_str());
	PyRun_SimpleString(command_2.c_str());
	PyRun_SimpleString(command_3.c_str());
	PyGILState_Release(gilstate);		
}

void Manta_API::addAdaptiveGrid(void * data, string gridName, string solverName, string type, int minX, int minY, int minZ, int maxX, int maxY, int maxZ)
{
	if (data == NULL)
	{
		cout << "NULL pointer passed to grid addAdaptiveGrid for grid " << gridName <<endl;
		return;
	}
	std::ostringstream stringStream;
	stringStream << "temp_" <<gridName;
	std::string temp_grid_name = stringStream.str();
	stringStream.str("");
	stringStream << temp_grid_name << " = "<< solverName << ".create(" << gridNameFromType(type) << ")";
	const std::string command_1 = stringStream.str();
	stringStream.str("");
	stringStream << temp_grid_name << ".readAdaptiveGridFromMemory(\'"<< data << "\', vec3(" << minX << "," << minY << "," << minZ << 
	"), vec3(" << maxX << "," << maxY << "," << maxZ << ") )";
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

void Manta_API::export_obstacles(float *data, int x, int y, int z)
{
	std::ostringstream stringStream;
	std::string grid_name = "obs_sdf";
	stringStream.str("");
	stringStream << grid_name << " = s.create(RealGrid)";
	const std::string command_1 = stringStream.str();
	stringStream.str("");
	stringStream << grid_name << ".readGridFromMemory(\'"<< data << "\', " << x << "," << y << "," << z << ")";
	const std::string command_2 = stringStream.str();
	const std::string command_3 = grid_name + ".applyToGrid(grid = flags, value = FlagObstacle)";
	PyGILState_STATE gilstate = PyGILState_Ensure();
	PyRun_SimpleString(command_1.c_str());
	PyRun_SimpleString(command_2.c_str());
	PyRun_SimpleString(command_3.c_str());
	PyGILState_Release(gilstate);		
}

void Manta_API::run_manta_sim_highRes(WTURBULENCE *wt)
{
	PyGILState_STATE gilstate = PyGILState_Ensure();
	int sim_frame = 1;
//	manta_write_effectors(fluid);
	std::string frame_str = static_cast<ostringstream*>( &(ostringstream() << sim_frame) )->str();
	std::string py_string_0 = string("sim_step_high(").append(frame_str);
	std::string py_string_1 = py_string_0.append(")\0");
	cout << "Debug C++: densityPointer:" << Manta_API::getGridPointer("density", "s")<<endl;
	PyRun_SimpleString("print ('pyhton density pointer:' + density.getDataPointer())");
	PyRun_SimpleString(py_string_1.c_str());
	cout<< "done"<<manta_sim_running<<endl;
	PyGILState_Release(gilstate);
	updateHighResPointers(wt);
}

void Manta_API::generate_manta_sim_file_highRes(SmokeModifierData *smd)
{
	string smoke_script = smoke_setup_high + smoke_step_high;		
	std::string final_script = parseScript(smoke_script, smd);
	PyGILState_STATE gilstate = PyGILState_Ensure();
	PyRun_SimpleString(final_script.c_str());
	PyGILState_Release(gilstate);
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
		ss <<  smd->domain->fluid->_xRes;
	else if (varName == "RESY")
		ss <<  smd->domain->fluid->_yRes;
	else if (varName == "RESZ")
		ss <<  smd->domain->fluid->_zRes;
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
	else if (varName == "XL_DENSITY_MEM")
		ss << smd->domain->wt->_densityBig;
	else if (varName == "XL_DENSITY_SIZE")
		ss << sizeof(float) * smd->domain->wt->_xResBig * smd->domain->wt->_yResBig * smd->domain->wt->_zResBig;
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

std::string Manta_API::parseScript(const string & setup_string, SmokeModifierData *smd)
{
//	ifstream f (file);
	std::istringstream f(setup_string);
//	ofstream of; /*PR: for Debug*/
	ostringstream res;
//	of.open("manta_scene.py", std::fstream::trunc);
	string line="";
//	if (f.is_open()){
		while(getline(f,line)){
//			of << parseLine(line,smd) << "\n";
			res << parseLine(line,smd) << "\n"; 
		}
//		f.close();
//	}
//	else{
//		printf ("Error: No scenario file found");
//	}
//	of.close();
	return res.str();
}

string Manta_API::getGridPointer(std::string gridName, std::string solverName)
{
	if ((gridName == "") && (solverName == "")){
		return "";
	}
	cout << "pointer to grid " << gridName << endl; 
#ifdef WITH_MANTA
	cout << "MANTA_DEFINED_________" << endl;
#else
	cout << "MANTA_NOT_DEFINED_________" << endl;
#endif
	
	PyGILState_STATE gilstate = PyGILState_Ensure();
	PyObject *main = PyImport_AddModule("__main__");
	if (main == NULL){cout << "null" << 1 << endl;}
    PyObject *globals = PyModule_GetDict(main);
    if (globals == NULL){cout << "null" << 12 << endl;}
    PyObject *grid_object = PyDict_GetItemString(globals, gridName.c_str());
    if (grid_object == NULL){cout << "null" << 13 << endl;}
    PyObject* func = PyObject_GetAttrString(grid_object,(char*)"getDataPointer");
    if (func == NULL){cout << "null" << 14 << endl;}
    PyObject* retured_value = PyObject_CallObject(func, NULL);
	PyObject* encoded = PyUnicode_AsUTF8String(retured_value);
	if (retured_value == NULL){cout << "null" << 15 << endl;}
	std::string res = strdup(PyBytes_AsString(encoded));
	cout << "RESRES" << res << "___" << endl;
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


void Manta_API::updatePointers(FLUID_3D *fluid)
{
	stringstream ss(getGridPointer("density", "s"));
	void *gridPointer = NULL;
	ss >> gridPointer;
	fluid->_density = (float* )gridPointer;
	ss.str("");
}

void Manta_API::updateHighResPointers(WTURBULENCE *wt)
{
	stringstream ss(getGridPointer("xl_density", "xl"));
	void *gridPointer = NULL;
	ss >> gridPointer;
	wt->_densityBig = (float* )gridPointer;
	ss.str("");
}

Manta_API::Manta_API(int *res, float dx, float dtdef, int init_heat, int init_fire, int init_colors,SmokeDomainSettings *sds): _xRes(res[0]), _yRes(res[1]), _zRes(res[2]), _res(0.0f)
{
	/*Here, we assume Python script has initalized the solver and all fields*/	
	
	//	// set simulation consts
	_dt = dtdef;	// just in case. set in step from a RNA factor
	_dx = dx;
	_totalCells   = _xRes * _yRes * _zRes;
	_slabSize = _xRes * _yRes;

//	
//	_iterations = 100;
//	_tempAmb = 0; 
//	_heatDiffusion = 1e-3;
//	_totalTime = 0.0f;
//	_totalSteps = 0;
//	_res = Vec3Int(_xRes,_yRes,_zRes);
//	_maxRes = MAX3(_xRes, _yRes, _zRes);
//	
//	// initialize wavelet turbulence
//	/*
//	 if(amplify)
//	 _wTurbulence = new WTURBULENCE(_res[0],_res[1],_res[2], amplify, noisetype);
//	 else
//	 _wTurbulence = NULL;
//	 */
//	
//	// scale the constants according to the refinement of the grid
//	if (!dx)
//		_dx = 1.0f / (float)_maxRes;
//	else
//		_dx = dx;
//	_constantScaling = 64.0f / _maxRes;
//	_constantScaling = (_constantScaling < 1.0f) ? 1.0f : _constantScaling;
//	_vorticityEps = 2.0f / _constantScaling; // Just in case set a default value
//	
//	// allocate arrays
	_xVelocity    = new float[_totalCells];
	_yVelocity    = new float[_totalCells];
	_zVelocity    = new float[_totalCells];
	_xVelocityOb  = new float[_totalCells];
	_yVelocityOb  = new float[_totalCells];
	_zVelocityOb  = new float[_totalCells];
//	_xVelocityOld = new float[_totalCells];
//	_yVelocityOld = new float[_totalCells];
//	_zVelocityOld = new float[_totalCells];
	_xForce       = new float[_totalCells];
	_yForce       = new float[_totalCells];
	_zForce       = new float[_totalCells];
	_density      = NULL ; //new float[_totalCells];
//	_densityOld   = new float[_totalCells];
	_obstacles    = new unsigned char[_totalCells]; // set 0 at end of step
//	
//	// For threaded version:
//	_xVelocityTemp = new float[_totalCells];
//	_yVelocityTemp = new float[_totalCells];
//	_zVelocityTemp = new float[_totalCells];
//	_densityTemp   = new float[_totalCells];
//	
//	// DG TODO: check if alloc went fine
//	
	for (int x = 0; x < _totalCells; x++)
	{
//		_densityOld[x]   = 0.0f;
		_xVelocity[x]    = 0.0f;
		_yVelocity[x]    = 0.0f;
		_zVelocity[x]    = 0.0f;
		_xVelocityOb[x]  = 0.0f;
		_yVelocityOb[x]  = 0.0f;
		_zVelocityOb[x]  = 0.0f;
//		_xVelocityOld[x] = 0.0f;
//		_yVelocityOld[x] = 0.0f;
//		_zVelocityOld[x] = 0.0f;
		_xForce[x]       = 0.0f;
		_yForce[x]       = 0.0f;
		_zForce[x]       = 0.0f;
		_obstacles[x]    = false;
	}
//	
//	/* heat */
//	_heat = _heatOld = _heatTemp = NULL;
//	if (init_heat) {
//		initHeat();
//	}
//	// Fire simulation
//	_flame = _fuel = _fuelTemp = _fuelOld = NULL;
//	_react = _reactTemp = _reactOld = NULL;
//	if (init_fire) {
//		initFire();
//	}
//	// Smoke color
//	_color_r = _color_rOld = _color_rTemp = NULL;
//	_color_g = _color_gOld = _color_gTemp = NULL;
//	_color_b = _color_bOld = _color_bTemp = NULL;
//	if (init_colors) {
//		initColors(0.0f, 0.0f, 0.0f);
//	}
//	
//	// boundary conditions of the fluid domain
//	// set default values -> vertically non-colliding
//	_domainBcFront = true;
//	_domainBcTop = false;
//	_domainBcLeft = true;
//	_domainBcBack = _domainBcFront;
//	_domainBcBottom = _domainBcTop;
//	_domainBcRight	= _domainBcLeft;
//	
//	_colloPrev = 1;	// default value
	
//	sds->fluid = this;
	generate_manta_sim_file_lowRes(sds->smd);
}

Manta_API::~Manta_API()
{
	if (_xVelocity) delete[] _xVelocity;
	if (_yVelocity) delete[] _yVelocity;
	if (_zVelocity) delete[] _zVelocity;
	if (_xVelocityOb) delete[] _xVelocityOb;
	if (_yVelocityOb) delete[] _yVelocityOb;
	if (_zVelocityOb) delete[] _zVelocityOb;
//	if (_xVelocityOld) delete[] _xVelocityOld;
//	if (_yVelocityOld) delete[] _yVelocityOld;
//	if (_zVelocityOld) delete[] _zVelocityOld;
	if (_xForce) delete[] _xForce;
	if (_yForce) delete[] _yForce;
	if (_zForce) delete[] _zForce;
	if (_density) delete[] _density;
//	if (_densityOld) delete[] _densityOld;
//	if (_heat) delete[] _heat;
//	if (_heatOld) delete[] _heatOld;
	if (_obstacles) delete[] _obstacles;
	
//	if (_xVelocityTemp) delete[] _xVelocityTemp;
//	if (_yVelocityTemp) delete[] _yVelocityTemp;
//	if (_zVelocityTemp) delete[] _zVelocityTemp;
//	if (_densityTemp) delete[] _densityTemp;
//	if (_heatTemp) delete[] _heatTemp;
	
//	if (_flame) delete[] _flame;
//	if (_fuel) delete[] _fuel;
//	if (_fuelTemp) delete[] _fuelTemp;
//	if (_fuelOld) delete[] _fuelOld;
//	if (_react) delete[] _react;
//	if (_reactTemp) delete[] _reactTemp;
//	if (_reactOld) delete[] _reactOld;
//	
//	if (_color_r) delete[] _color_r;
//	if (_color_rOld) delete[] _color_rOld;
//	if (_color_rTemp) delete[] _color_rTemp;
//	if (_color_g) delete[] _color_g;
//	if (_color_gOld) delete[] _color_gOld;
//	if (_color_gTemp) delete[] _color_gTemp;
//	if (_color_b) delete[] _color_b;
//	if (_color_bOld) delete[] _color_bOld;
//	if (_color_bTemp) delete[] _color_bTemp;
	
    // printf("deleted fluid\n");
}
