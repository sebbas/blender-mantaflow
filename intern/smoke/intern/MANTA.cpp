#include "MANTA.h"
#include "WTURBULENCE.h"
#include "scenarios/smoke.h"
void runMantaScript(vector<string>& args);//defined in manta_pp/pwrapper/pymain.cpp

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

void wavelets_add_lowres_density(SmokeDomainSettings *sds)
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
    /*! legacy headers for reading old files */
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

void indent_ss(stringstream& ss, int indent)
{
	/*two-spaces indent*/
	if (indent < 0) return;
	std::string indentation = ""; 
	for (size_t cnt(0); cnt < indent; ++cnt) {
		indentation += "  ";
	}
	ss << indentation;
}

void manta_gen_noise(stringstream& ss, char* solver, int indent, char *noise, int seed, bool load, bool clamp, float clampNeg, float clampPos, float valScale, float valOffset, float timeAnim)
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

void manta_solve_pressure(stringstream& ss, char *flags, char *vel, char *pressure, bool useResNorms, int openBound, int solver_res,float cgMaxIterFac, float cgAccuracy)
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

void manta_advect_SemiLagr(stringstream& ss, int indent, char *flags, char *vel, char *grid, int order)
{
	if((order <=1) || (flags == NULL) || (vel == NULL) || (grid == NULL)){return;}
	indent_ss(ss, indent);
	ss << "advectSemiLagrange(flags=" << flags << ", vel=" << vel \
	<< ", grid=" << grid << ", order=" << order << ") \n"; 
}

/*create solver, handle 2D case*/
void manta_create_solver(stringstream& ss, char *name, char *nick, char *grid_size_name, int x_res, int y_res, int z_res, int dim)
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
void add_mesh_transform_method(stringstream& ss)
{
	ss << "def transform_back(obj, gs):\n" <<
	"  obj.scale(gs/2)\n" <<
	"  obj.offset(gs/2)\n\n";
}

void manta_cache_path(char *filepath)
{
	char *name="manta";
	BLI_make_file_string("/", filepath, BLI_temporary_dir(), name);
}

//void BLI_dir_create_recursive(const char *filepath);
void create_manta_folder()
{
	char* filepath=NULL;
	manta_cache_path(filepath);
	//BLI_dir_create_recursive(filepath);
	
}

void *run_manta_scene_thread(void *arguments)
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

void run_manta_scene(Scene *s, SmokeModifierData *smd)
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
	int rc = pthread_create(&manta_thread, NULL, run_manta_sim_thread, (void *)args);
//	pthread_join(manta_thread,NULL);
//	pthread_detach(manta_thread);
}

void stop_manta_sim()
{
	pthread_cancel(manta_thread);
}


void *run_manta_sim_thread(void *arguments)
{
	struct manta_arg_struct *args = (struct manta_arg_struct *)arguments;
	SmokeModifierData *smd = &args->smd;
	Scene *s = &args->s;
	int num_sim_steps = smd->domain->manta_end_frame - smd->domain->manta_start_frame + 1;
	smd->domain->manta_sim_frame = 0;
	PyGILState_STATE gilstate = PyGILState_Ensure();
	for (int fr=0; fr< num_sim_steps; ++fr) {
		if(smd->domain->manta_sim_frame == -1)
			break;
		printf("Simulation Step");
		manta_write_effectors(s, smd);
		smd->domain->manta_sim_frame = fr;
		std::string frame_str = static_cast<ostringstream*>( &(ostringstream() << fr) )->str();
		std::string py_string_0 = string("sim_step(").append(frame_str);
		std::string py_string_1 = py_string_0.append(")\0");
		PyRun_SimpleString(py_string_1.c_str());
		cout<< "done"<<manta_sim_running<<endl;
	}
	//returning simulation state to "not simulating" aka -1
	smd->domain->manta_sim_frame = -1;
	PyGILState_Release(gilstate);
}

void generate_manta_sim_file(Scene *scene, SmokeModifierData *smd)
{
	 /*create python file with 2-spaces indentation*/
	bool wavelets = smd->domain->flags & MOD_SMOKE_HIGHRES;
	bool noise_clamp = smd->domain->flags & MOD_SMOKE_NOISE_CLAMP; 
	float noise_clamp_neg = smd->domain->noise_clamp_neg;
	float noise_clamp_pos = smd->domain->noise_clamp_pos;
	float noise_val_scale = smd->domain->noise_val_scale;
	float noise_val_offset = smd->domain->noise_val_offset;
	float noise_time_anim = smd->domain->noise_time_anim;
	int num_sim_frames = smd->domain->manta_end_frame - smd->domain->manta_start_frame + 1;
	if(num_sim_frames < 1)
		return;
	ofstream manta_setup_file;
//	manta_setup_file.open("manta_scene.py", std::fstream::trunc);
	stringstream ss; /*setup contents*/
	
	/*header*/
	ss << "from manta import * \n";
	ss << "import os, shutil, math, sys \n";
	if (!file_exists("manta_flow.obj")){
		return;
	}
	add_mesh_transform_method(ss);
	/*Data Declaration*/
	/*Wavelets variables*/
	int upres = smd->domain->amplify+1;
	ss << "uvs = " << smd->domain->manta_uvs_num << "\n";
	ss << "velInflow = vec3(0, 0, 1)"<< "\n";	/*TODO:add UI*/
	if (wavelets) {
		ss << "upres = " << upres << "\n";
		ss << "wltStrength = " << smd->domain->strength << "\n";
		if(upres > 0)/*TODO:add UI*/
		{	ss << "octaves = int( math.log(upres)/ math.log(2.0) + 0.5 ) \n";	}
		else
		{	ss << "octaves = 0"<< "\n";	}
	}
	else upres = 0;
	
	/*Solver Resolution*/
	ss << "res = " << smd->domain->maxres << " \n";
	manta_create_solver(ss, "s", "main", "gs", smd->domain->fluid->xRes(), smd->domain->fluid->yRes(), smd->domain->fluid->zRes(), smd->domain->manta_solver_res);
	ss << "s.timestep = " << smd->domain->time_scale << " \n";
	
	/*Noise Field*/
	manta_gen_noise(ss, "s", 0, "noise", 256, true, noise_clamp, noise_clamp_neg, noise_clamp_pos, noise_val_scale, noise_val_offset, noise_time_anim);
	
	/*Inflow sources*/
	ss << "source = s.create(Mesh)\n";
	ss << "source.load('manta_flow.obj')\n";
	ss << "transform_back(source, gs)\n";
	ss << "sourceVel = s.create(Mesh)\n";
	ss << "sourceVel.load('manta_flow.obj')\n";
	ss << "transform_back(sourceVel, gs)\n";
	
	/*Wavelets: larger solver*/
	if(wavelets && upres>0)
	{
		manta_create_solver(ss, "xl", "larger", "xl_gs", smd->domain->fluid->xRes() * upres, smd->domain->fluid->yRes()* upres, smd->domain->fluid->zRes() * upres, smd->domain->manta_solver_res);
		ss << "xl.timestep = " << smd->domain->time_scale * 0.5f<< " \n";
		
		ss << "xl_vel = xl.create(MACGrid) \n";
		ss << "xl_density = xl.create(RealGrid) \n";/*smoke simulation*/
		ss << "xl_flags = xl.create(FlagGrid) \n";
		ss << "xl_flags.initDomain() \n";
		ss << "xl_flags.fillGrid() \n";		
		ss << "xl_source = s.create(Mesh)\n";
		ss << "xl_source.load('manta_flow.obj')\n";
		ss << "transform_back(xl_source, gs)\n";
		
		/*High-res obstacle handling*/
		if (file_exists("manta_coll.obj"))
		{
			ss << "xl_obs = s.create(Mesh)\n";
			ss << "xl_obs.load('manta_coll.obj')\n";
			ss << "transform_back(xl_obs, gs)\n";
			ss << "xl_obs.applyToGrid(grid=xl_flags, value=FlagObstacle,cutoff=3)\n";
		}
		manta_gen_noise(ss, "xl", 0, "xl_noise", 256, true, noise_clamp, noise_clamp_neg, noise_clamp_pos, noise_val_scale, noise_val_offset, noise_time_anim * (float)upres);
	}
	/*Flow setup*/
	ss << "flags = s.create(FlagGrid) \n";/*must always be present*/
	ss << "flags.initDomain() \n";
	ss << "flags.fillGrid() \n";
	/*Obstacle handling*/
	if (file_exists("manta_coll.obj"))
	{
		ss << "obs = s.create(Mesh)\n";
		ss << "obs.load('manta_coll.obj')\n";
		ss << "transform_back(obs, gs)\n";
		ss << "obs.applyToGrid(grid=flags, value=FlagObstacle, cutoff=3)\n";
		ss << "sdf_obs  = s.create(LevelsetGrid)\n";
		ss << "obs.meshSDF(obs, sdf_obs, 1.1)\n";
	}
	/*Create the array of UV grids*/
	if(wavelets){
		ss << "uv = [] \n";
		ss << "for i in range(uvs): \n";
		ss << "  uvGrid = s.create(VecGrid) \n";
		ss << "  uv.append(uvGrid) \n";
		ss << "  resetUvGrid( uv[i] ) \n";
	}
	/*Grids setup*/
	/*For now, only one grid of each kind is needed*/
	ss << "vel = s.create(MACGrid) \n";
	ss << "density = s.create(RealGrid) \n";/*smoke simulation*/
	ss << "pressure = s.create(RealGrid) \n";/*must always be present*/
	ss << "energy = s.create(RealGrid) \n";
	ss << "tempFlag  = s.create(FlagGrid)\n";
	ss << "sdf_flow  = s.create(LevelsetGrid)\n";
	ss << "forces = s.create(MACGrid)\n";
	
	ss << "source.meshSDF(source, sdf_flow, 1.1)\n";
	ss << "source_shape = s.create(Cylinder, center=gs*vec3(0.5,0.1,0.5), radius=res*0.14, z=gs*vec3(0, 0.02, 0))\n";
	/*Wavelets noise field*/
	if (wavelets){
		ss << "xl_wltnoise = s.create(NoiseField, loadFromFile=True) \n";
		ss << "xl_wltnoise.posScale = vec3( int(1.0*gs.x) ) * 0.5 \n";
		/*scale according to lowres sim , smaller numbers mean larger vortices*/
		if (upres > 0){
			ss << "xl_wltnoise.posScale = xl_wltnoise.posScale * " << (1./(float)upres) << ((upres == 1)?". \n":"\n");
		}
		ss << "xl_wltnoise.timeAnim = 0.1 \n";
	}
	
	/*Flow solving steps, main loop*/
	ss << "def sim_step(t):\n";
	ss << "  forces.load('manta_forces.uni')\n";
	ss << "  addForceField(flags=flags, vel=vel,force=forces)\n";
	ss << "  addBuoyancy(density=density, vel=vel, gravity=vec3(0,0," << (-smd->domain->beta) << "), flags=flags) \n";
	
	manta_advect_SemiLagr(ss, 1, "flags", "vel", "density", 2);
	manta_advect_SemiLagr(ss, 1, "flags", "vel", "vel", 2);
	
	if(wavelets){
		ss << "  for i in range(uvs): \n";
		manta_advect_SemiLagr(ss, 2, "flags", "vel", "uv[i]", 2);
		ss << "    updateUvWeight( resetTime=16.5 , index=i, numUvs=uvs, uv=uv[i] )\n"; 
	}
	ss << "  applyInflow=False\n";
	ss << "  if (t>=0 and t<75):\n";
	if (noise_val_scale > 0.)
		ss << "    densityInflowMeshNoise( flags=flags, density=density, noise=noise, mesh=source, scale=3, sigma=0.5 )\n";
	else
		ss << "    densityInflowMesh(flags=flags, density=density, mesh=source, value=1)\n";	
	ss << "    applyInflow=True\n";
	
	ss << "  setWallBcs(flags=flags, vel=vel) \n";
//	ss << "  addBuoyancy(density=density, vel=vel, gravity=vec3(0,-6e-4,0), flags=flags) \n";
	ss << "  vorticityConfinement( vel=vel, flags=flags, strength=" << smd->domain->vorticity / 10. << " ) \n";
	
	manta_solve_pressure(ss,"flags", "vel", "pressure",true,smd->domain->border_collisions, smd->domain->manta_solver_res,1.0,0.01);
	ss << "  setWallBcs(flags=flags, vel=vel) \n";
	
	/*determine weighting*/
	ss << "  computeEnergy(flags=flags, vel=vel, energy=energy)\n";
	/* mark outer obstacle region by extrapolating flags for 2 layers */
	ss << "  tempFlag.copyFrom(flags)\n";
	ss << "  extrapolateSimpleFlags( flags=flags, val=tempFlag, distance=2, flagFrom=FlagObstacle, flagTo=FlagFluid )\n";
	/*now extrapolate energy weights into obstacles to fix boundary layer*/
	ss << "  extrapolateSimpleFlags( flags=tempFlag, val=energy, distance=6, flagFrom=FlagFluid, flagTo=FlagObstacle )\n";
	ss << "  computeWaveletCoeffs(energy)\n";
	/*Saving output*/
	ss << "  density.save('den%04d_temp.uni' % t) \n";
	ss << "  os.rename('den%04d_temp.uni' % t, 'den%04d.uni' % t) \n";
	ss << "  s.step()\n";
	ss << " \n";
	
	if (wavelets && upres > 0){
		ss << "  interpolateMACGrid( source=vel, target=xl_vel ) \n";
		/*add all necessary octaves*/
		ss << "  sStr = 1.0 * wltStrength  \n";
		ss << "  sPos = 2.0  \n";
		ss << "  for o in range(octaves): \n";
		ss << "    for i in range(uvs): \n";
		ss << "      uvWeight = getUvWeight(uv[i])  \n";
		ss << "      applyNoiseVec3( flags=xl_flags, target=xl_vel, noise=xl_wltnoise, scale=sStr * uvWeight, scaleSpatial=sPos , weight=energy, uv=uv[i] ) \n";
		ss << "    sStr *= 0.06 # magic kolmogorov factor \n";
		ss << "    sPos *= 2.0 \n";
		ss << "  for substep in range(upres):  \n";
		ss << "    advectSemiLagrange(flags=xl_flags, vel=xl_vel, grid=xl_density, order=2)  \n";
		ss << "  if (applyInflow): \n";
		if (noise_val_scale > 0.)
			ss << "    densityInflowMeshNoise( flags=xl_flags, density=xl_density, noise=xl_noise, mesh=xl_source, scale=3, sigma=0.5 ) \n";
		else
			ss << "    densityInflowMesh(flags=xl_flags, density=xl_density, mesh=source, value=1)\n";
		ss << "  xl_density.save('densityXl_%04d.uni' % t)\n";
		//ss << "    densityInflow( flags=xl_flags, density=xl_density, noise=xl_noise, shape=xl_source, scale=1, sigma=0.5 ) \n";
		ss << "  xl.step()   \n";
	}
//	manta_setup_file << ss.rdbuf();
//	manta_setup_file.close();
	/*constrcting final setup*/
	string smoke_script = smoke_setup_low + ((wavelets)?smoke_setup_high:"") + smoke_step_low + ((wavelets)?smoke_step_high:"");	
	manta_setup_file << smoke_script ;
	manta_setup_file.close();
	parseFile(smd, scene, smoke_script);
	vector<string> a;
	a.push_back("manta_scene.py");
	runMantaScript("",a);
}

std::string getRealValue(SmokeModifierData *smd, Scene *s, const std::string& varName)
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
	else if (varName == "XL_TIMESTEP")	/*add support*/
		ss << smd->domain->noise_time_anim * smd->domain->amplify+1;/*lowres time * upres*/
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
	else 
		cout<< "ERROR: Unknown option:"<< varName <<endl; 
	return ss.str();
}

std::string parseLine(SmokeModifierData *smd, Scene *s, const string& line)
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
			res 		+= getRealValue(smd,s,line.substr(start_del, currPos - start_del));
		}
		currPos ++;
	}
	res += line.substr(end_del+1, line.size()- end_del);
	return res;
}

void parseFile(SmokeModifierData *smd, Scene *s, const string & setup_string)
{
//	ifstream f (file);
std::istringstream f(setup_string);
	ofstream of;
	of.open("manta_scene.py", std::fstream::trunc);
	string line="";
//	if (f.is_open()){
		while(getline(f,line)){
			of << parseLine(smd,s,line) << "\n";
		}
//		f.close();
//	}
//	else{
//		printf ("Error: No scenario file found");
//	}
	of.close();
}


