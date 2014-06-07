#ifndef MANTA_H
#define MANTA_H
#include "FLUID_3D.h"
#include "zlib.h"
#include "../../../source/blender/makesdna/DNA_scene_types.h"
#include "../../../source/blender/makesdna/DNA_modifier_types.h"
#include "../../../source/blender/makesdna/DNA_smoke_types.h"
#include <sstream>
#include <fstream>

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

static void manta_gen_noise(stringstream& ss, int indent, char *noise, int seed, bool load, bool clamp, int clampNeg, int clampPos, float valScale, float valOffset, float timeAnim)
{
	if (ss == NULL)/*should never be here*/
	{
		return;
	}
	std::string indentation = ""; 
	for (size_t cnt(0); cnt < indent; ++cnt) {
		indentation += "  ";/*two-spaces indent*/
	}
	ss << indentation << noise << " = s.create(NoiseField, fixedSeed=" << seed << ", loadFromFile="<< (load?"True":"False") <<") \n";
	ss << indentation << noise << ".posScale = vec3(20) \n";
	ss << indentation << noise << ".clamp = " << ((clamp)?"True":"False") << " \n";
	ss << indentation << noise << ".clampNeg = " << clampNeg << " \n";
	ss << indentation << noise << ".clampPos = " << clampPos << " \n";
	ss << indentation << noise << ".valScale = " << valScale << " \n";
	ss << indentation << noise << ".valOffset = " << valOffset << " \n";
	ss << indentation << noise << ".timeAnim = " << timeAnim << " \n";
}

static void manta_solve_pressure(stringstream& ss, char *flags, char *vel, char *pressure, bool useResNorms, int openBound, int solver_res)
{
	/*open:0 ; vertical : 1; closed:2*/
	ss << "  solvePressure(flags=" << flags << ", vel=" << vel << ", pressure=" << pressure << ", useResNorm=" << (useResNorms?"True":"False") << ", openBound='";	
	
	if(openBound == 1) /*vertical*/
	{
		ss << "yY') \n";
	}
	else if (openBound == 0) /*open*/
	{
		if(solver_res == 2)
			ss << "xXyY') \n";
		else
			ss << "xXyYzZ') \n";
	}
	else	/*also for closed bounds*/ 
	{
			ss << "') \n";
	}
}

static void manta_advect_SemiLagr(stringstream& ss, char *indent, char *flags, char *vel, char *grid, int order)
{
	if((order <=1) || (indent == NULL) || (flags == NULL) || (vel == NULL) || (grid == NULL))
	{return;}
	ss << indent << "advectSemiLagrange(flags=" << flags << ", vel=" << vel \
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

static void generate_manta_sim_file(Scene *scene, SmokeModifierData *smd)
{
	/*for now, simpleplume file creation
	*create python file with 2-spaces indentation*/
	
	bool wavelets = smd->domain->flags & MOD_SMOKE_HIGHRES;
	FLUID_3D *fluid = smd->domain->fluid;
	
	ofstream manta_setup_file;
	manta_setup_file.open("manta_scene.py", std::fstream::trunc);
	stringstream ss; /*setup contents*/
	
	/*header*/
	ss << "from manta import * \n";
	ss << "import os, shutil, math, sys \n";
	
/*Data Declaration*/
	/*Wavelets variables*/
	int upres = smd->domain->amplify;
	if (wavelets) {
		ss << "upres = " << upres << "\n";
		ss << "wltStrength = " << smd->domain->strength << "\n";
		ss << "uvs = 1" << "\n";					/*TODO:add UI*/
		ss << "velInflow = vec3(2, 0, 0)"<< "\n";	/*TODO:add UI*/
		if(smd->domain->amplify > 0)/*TODO:add UI*/
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
	
/*Grids setup*/
/*For now, only one grid of each kind is needed*/
	ss << "vel = s.create(MACGrid) \n";
	ss << "density = s.create(RealGrid) \n";/*smoke simulation*/
	ss << "pressure = s.create(RealGrid) \n";/*must always be present*/
	if(wavelets){
		ss << "energy = s.create(RealGrid) \n";
		ss << "tempFlag  = s.create(FlagGrid)\n";
	}
/*Noise Field*/
	manta_gen_noise(ss, 0, "noise", 256, true, true, 0, 1, 1, 0.75, 0.2);

/*Wavelets: larger solver*/
	if(wavelets && upres>0)
	{
		manta_create_solver(ss, "xl", "larger", "xl_gs", fluid->xRes(), fluid->zRes(), fluid->yRes(), smd->domain->manta_solver_res);
		ss << "xl.timestep = " << smd->domain->time_scale << " \n";
		
		ss << "xl_vel = s.create(MACGrid) \n";
		ss << "xl_density = s.create(RealGrid) \n";/*smoke simulation*/
		ss << "xl_flags = s.create(FlagGrid) \n";
		ss << "xl_flags.initDomain() \n";
		ss << "xl_flags.fillGrid() \n";
			
		ss << "xl_source = xl.create(Cylinder, center=xl_gs*vec3(0.3,0.2,0.5), radius=xl_gs.x*0.081, z=xl_gs*vec3(0.081, 0, 0)) \n";
		ss << "xl_obs    = xl.create(Sphere,   center=xl_gs*vec3(0.5,0.5,0.5), radius=xl_gs.x*0.15) \n";
		ss << "xl_obs.applyToGrid(grid=xl_flags, value=FlagObstacle) \n";
		manta_gen_noise(ss, 0, "xl_noise", 256, true, true, 0, 2, 1, 0.075, 0.3 * upres);
	}
/*Flow setup*/
	ss << "flags = s.create(FlagGrid) \n";/*must always be present*/
	ss << "flags.initDomain() \n";
	ss << "flags.fillGrid() \n";

/*GUI for debugging purposes*/
	ss << "if (GUI):\n  gui = Gui()\n  gui.show() \n";

/*Inflow source - for now, using mock sphere */
	ss << "source    = s.create(Cylinder, center=gs*vec3(0.3,0.2,0.5), radius=res*0.081, z=gs*vec3(0.081, 0, 0))\n";
	ss << "sourceVel = s.create(Cylinder, center=gs*vec3(0.3,0.2,0.5), radius=res*0.15 , z=gs*vec3(0.15 , 0, 0))\n";
	ss << "obs       = s.create(Sphere,   center=gs*vec3(0.5,0.5,0.5), radius=res*0.15)\n";
	ss << "obs.applyToGrid(grid=flags, value=FlagObstacle)\n";
	
/*Flow solving stepsv, main loop*/
	ss << "for t in xrange(" << scene->r.sfra << ", " << scene->r.efra << "): \n";
	ss << "  densityInflow(flags=flags, density=density, noise=noise, shape=source, scale=1, sigma=0.5) \n"	;
	manta_advect_SemiLagr(ss, "  ", "flags", "vel", "density", 2);
	manta_advect_SemiLagr(ss, "  ", "flags", "vel", "vel", 2);
	ss << "  setWallBcs(flags=flags, vel=vel) \n";
	ss << "  addBuoyancy(density=density, vel=vel, gravity=vec3(0,-6e-4,0), flags=flags) \n";
	manta_solve_pressure(ss,"flags", "vel", "pressure",true,smd->domain->border_collisions, smd->domain->manta_solver_res);
	ss << "  setWallBcs(flags=flags, vel=vel) \n";

/*Saving output*/
	ss << "  density.save('den%04d.uni' % t) \n";
	ss << "  s.step()\n";
	ss << " \n";
	
	manta_setup_file << ss.rdbuf();
	manta_setup_file.close();		
}

#endif /* MANTA_H */

