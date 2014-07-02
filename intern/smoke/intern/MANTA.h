#ifndef MANTA_H
#define MANTA_H
#include "FLUID_3D.h"
#include "zlib.h"
#include "../../../source/blender/makesdna/DNA_scene_types.h"
#include "../../../source/blender/makesdna/DNA_modifier_types.h"
#include "../../../source/blender/makesdna/DNA_smoke_types.h"
#include <sstream>
#include <stdlib.h>
#include <fstream>
#include <pthread.h>
#include <Python.h>
#include "../../../extern/manta_pp/pwrapper/pymain.cpp"
extern "C" bool manta_check_grid_size(struct FLUID_3D *fluid, int dimX, int dimY, int dimZ)
{
	if (!(dimX == fluid->xRes() && dimY == fluid->yRes() && dimZ == fluid->zRes())) {
		for (int cnt(0); cnt < fluid->_totalCells; cnt++)
			fluid->_density[cnt] = 0.0f;
		return false;
	}
	return true;
}

extern "C" void read_mantaflow_sim(struct FLUID_3D *fluid, char *name)
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
		for (int cnt(0); cnt < fluid->_totalCells; cnt++)
			fluid->_density[cnt] = 0.0f;
		return;
	}
	
    char ID[5] = {0,0,0,0,0};
	gzread(gzf, ID, 4);
	
	/* legacy file format */
    if (!strcmp(ID, "DDF2")) {
        UniLegacyHeader head;
		gzread(gzf, &head, sizeof(UniLegacyHeader));
		if (!manta_check_grid_size(fluid, head.dimX, head.dimY, head.dimZ))	return;
        int numEl = head.dimX*head.dimY*head.dimZ;
        gzseek(gzf, numEl, SEEK_CUR);
        /* actual grid read */
        gzread(gzf, fluid->_density, sizeof(float)*numEl);
    } 
	/* legacy file format 2 */
    else if (!strcmp(ID, "MNT1")) {
        UniLegacyHeader2 head;
        gzread(gzf, &head, sizeof(UniLegacyHeader2));
		if (!manta_check_grid_size(fluid, head.dimX, head.dimY, head.dimZ))	return;
        /* actual grid read*/
        gzread(gzf, fluid->_density, sizeof(float)*head.dimX*head.dimY*head.dimZ);
    }
	/* current file format*/
    else if (!strcmp(ID, "MNT2")) {
        UniHeader head;
        gzread(gzf, &head, sizeof(UniHeader));
		if (!manta_check_grid_size(fluid, head.dimX, head.dimY, head.dimZ))	return;
		/* actual grid read */
        gzread(gzf,fluid->_density, sizeof(float)*head.dimX*head.dimY*head.dimZ);
    }
    gzclose(gzf);

#	endif	/*zlib*/
 }

static void indent_ss(stringstream& ss, int indent)
{
	/*two-spaces indent*/
	if (indent < 0) return;
	std::string indentation = ""; 
	for (size_t cnt(0); cnt < indent; ++cnt) {
		indentation += "  ";
	}
	ss << indentation;
}

static void manta_gen_noise(stringstream& ss, char* solver, int indent, char *noise, int seed, bool load, bool clamp, float clampNeg, float clampPos, float valScale, float valOffset, float timeAnim)
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

static void manta_solve_pressure(stringstream& ss, char *flags, char *vel, char *pressure, bool useResNorms, int openBound, int solver_res,float cgMaxIterFac=1.0, float cgAccuracy = 0.01)
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

static void manta_advect_SemiLagr(stringstream& ss, int indent, char *flags, char *vel, char *grid, int order)
{
	if((order <=1) || (flags == NULL) || (vel == NULL) || (grid == NULL)){return;}
	indent_ss(ss, indent);
	ss << "advectSemiLagrange(flags=" << flags << ", vel=" << vel \
	<< ", grid=" << grid << ", order=" << order << ") \n"; 
}

/*create solver, handle 2D case*/
static void manta_create_solver(stringstream& ss, char *name, char *nick, char *grid_size_name, int x_res, int y_res, int z_res, int dim)
{
	if ((dim != 2) && (dim != 3))
	{ return; }
	if (dim == 2)
	{ z_res = 1; }
	ss << grid_size_name << " = vec3(" << x_res << ", " << y_res << ", " << z_res << ")" << " \n";
	ss << name << " = Solver(name = '" << nick << "', gridSize = " << grid_size_name << ", dim = " << dim << ") \n";
}

inline bool file_exists (const std::string& name) {
    return ( access( name.c_str(), F_OK ) != -1 );
}

/*blender transforms obj coords to [-1,1]. This method transforms them back*/
static void add_mesh_transform_method(stringstream& ss)
{
	ss << "def transform_back(obj, res):\n" <<
	"  obj.scale(vec3(res/2, res/2, res/2))\n" <<
	"  obj.offset(vec3(res/2, res/2, res/2))\n\n";
}

//void *run_manta_scene(void *threadid)
void run_manta_scene()
{
	int a = Py_IsInitialized();
	//PyInterpreterState *st = PyThreadState_GET()->interp;
	//PyThreadState *ts = Py_NewInterpreter();
		
	vector<string> args;
	args.push_back("manta_scene.py");
	
	runScript(args);
	 
	//system("./manta manta_scene.py");
//	pthread_exit(NULL);
}

static void generate_manta_sim_file(Scene *scene, SmokeModifierData *smd)
{
	/*for now, simpleplume file creation
	*create python file with 2-spaces indentation*/
	
	bool wavelets = smd->domain->flags & MOD_SMOKE_HIGHRES;
	bool noise_clamp = smd->domain->flags & MOD_SMOKE_NOISE_CLAMP; 
	float noise_clamp_neg = smd->domain->noise_clamp_neg;
	float noise_clamp_pos = smd->domain->noise_clamp_pos;
	float noise_val_scale = smd->domain->noise_val_scale;
	float noise_val_offset = smd->domain->noise_val_offset;
	float noise_time_anim = smd->domain->noise_time_anim;
	
	FLUID_3D *fluid = smd->domain->fluid;
	ofstream manta_setup_file;
	manta_setup_file.open("manta_scene.py", std::fstream::trunc);
	stringstream ss; /*setup contents*/
	
	/*header*/
	ss << "import manta\n";//"from manta import * \n";
	ss << "import os, shutil, math, sys \n";
	if (!file_exists("manta_flow.obj")){
		return;
	}
	add_mesh_transform_method(ss);
/*Data Declaration*/
	/*Wavelets variables*/
	int upres = smd->domain->amplify;
	ss << "uvs = " << smd->domain->manta_uvs_num << "\n";
	ss << "velInflow = vec3(2, 0, 0)"<< "\n";	/*TODO:add UI*/
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
		/*Z axis in Blender = Y axis in Mantaflow*/
	manta_create_solver(ss, "s", "main", "gs", fluid->xRes(), fluid->zRes(), fluid->yRes(), smd->domain->manta_solver_res);
	ss << "s.timestep = " << smd->domain->time_scale << " \n";
		
/*Noise Field*/
	manta_gen_noise(ss, "s", 0, "noise", 256, true, noise_clamp, noise_clamp_neg, noise_clamp_pos, noise_val_scale, noise_val_offset, noise_time_anim);

/*Inflow source - for now, using mock sphere */
	ss << "source = s.create(Mesh)\n";
	ss << "source.load('manta_flow.obj')\n";
	ss << "transform_back(source, res)\n";
	ss << "sourceVel = s.create(Mesh)\n";
	ss << "sourceVel.load('manta_flow.obj')\n";
	ss << "transform_back(sourceVel, res)\n";
//	ss << "source    = s.create(Cylinder, center=gs*vec3(0.3,0.2,0.5), radius=res*0.081, z=gs*vec3(0.081, 0, 0))\n";
//	ss << "sourceVel = s.create(Cylinder, center=gs*vec3(0.3,0.2,0.5), radius=res*0.15 , z=gs*vec3(0.15 , 0, 0))\n";
	
/*Wavelets: larger solver*/
	if(wavelets && upres>0)
	{
		manta_create_solver(ss, "xl", "larger", "xl_gs", fluid->xRes() * upres, fluid->zRes()* upres, fluid->yRes() * upres, smd->domain->manta_solver_res);
		ss << "xl.timestep = " << smd->domain->time_scale << " \n";
		
		ss << "xl_vel = xl.create(MACGrid) \n";
		ss << "xl_density = xl.create(RealGrid) \n";/*smoke simulation*/
		ss << "xl_flags = xl.create(FlagGrid) \n";
		ss << "xl_flags.initDomain() \n";
		ss << "xl_flags.fillGrid() \n";
			
//		ss << "xl_source = xl.create(Cylinder, center=xl_gs*vec3(0.3,0.2,0.5), radius=xl_gs.x*0.081, z=xl_gs*vec3(0.081, 0, 0)) \n";
		ss << "xl_source = s.create(Mesh)\n";
		ss << "xl_source.load('manta_flow.obj')\n";
		ss << "transform_back(xl_source, res)\n";
//		ss << "xl_source.scale(vec3("<< upres <<", " << upres <<", " << upres << "))\n";

		/*Obstacle handling*/
		if (file_exists("manta_coll.obj"))
		{
			ss << "xl_obs = s.create(Mesh)\n";
			ss << "xl_obs.load('manta_coll.obj')\n";
			ss << "transform_back(xl_obs, res)\n";
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
		ss << "transform_back(obs, res)\n";
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
	ss << "source.meshSDF(source, sdf_flow, 1.1)\n";

	/*Wavelets noise field*/
	if (wavelets)
	{
		ss << "xl_wltnoise = s.create(NoiseField, loadFromFile=True) \n";
		ss << "xl_wltnoise.posScale = vec3( int(1.0*gs.x) ) * 0.5 \n";
		/*scale according to lowres sim , smaller numbers mean larger vortices*/
		if (upres > 0){
			ss << "xl_wltnoise.posScale = xl_wltnoise.posScale * " << (1./(float)upres) << ((upres == 1)?". \n":"\n");
		}
		ss << "xl_wltnoise.timeAnim = 0.1 \n";
	}
	
/*GUI for debugging purposes*/
	ss << "if (GUI):\n  gui = Gui()\n  gui.show() \n";

/*Flow solving stepsv, main loop*/
	ss << "for t in xrange(" << scene->r.sfra << ", " << scene->r.efra << "): \n";
	manta_advect_SemiLagr(ss, 1, "flags", "vel", "density", 2);
	manta_advect_SemiLagr(ss, 1, "flags", "vel", "vel", 2);
	
	if(wavelets){
		ss << "  for i in range(uvs): \n";
		manta_advect_SemiLagr(ss, 2, "flags", "vel", "uv[i]", 2);
		ss << "    updateUvWeight( resetTime=16.5 , index=i, numUvs=uvs, uv=uv[i] )\n"; 
	}
	ss << "  applyInflow=False\n";
	ss << "  if (t>=0 and t<75):\n";
		ss << "    densityInflowMesh( flags=flags, density=density, noise=noise, mesh=source, scale=1, sigma=0.5 )\n";
		//ss << "    densityInflow( flags=flags, density=density, noise=noise, shape=source, scale=1, sigma=0.5 )\n";
		ss << "    sourceVel.applyToGrid(grid=vel , value=velInflow,cutoff = 3)\n";
		//ss << "    sourceVel.applyToGrid( grid=vel , value=velInflow )\n";
	ss << "    applyInflow=True\n";
	
	ss << "  setWallBcs(flags=flags, vel=vel) \n";
	ss << "  addBuoyancy(density=density, vel=vel, gravity=vec3(0,-6e-4,0), flags=flags) \n";
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
	ss << "  density.save('den%04d.uni' % t) \n";
	ss << "  s.step()\n";
	ss << " \n";
	
	/**/
	if (wavelets && upres > 0)
	{
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
		ss << "    densityInflowMesh( flags=xl_flags, density=xl_density, noise=xl_noise, mesh=xl_source, scale=1, sigma=0.5 ) \n";
		//ss << "    densityInflow( flags=xl_flags, density=xl_density, noise=xl_noise, shape=xl_source, scale=1, sigma=0.5 ) \n";
		ss << "  xl.step()   \n";
	}
	manta_setup_file << ss.rdbuf();
	manta_setup_file.close();
	run_manta_scene();
//	pthread_t manta_thread;
//	int rc = pthread_create(&manta_thread, NULL, run_manta_scene, NULL);
//	pthread_detach(manta_thread);
}

#endif /* MANTA_H */

