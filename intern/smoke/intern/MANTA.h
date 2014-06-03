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

static void manta_gen_noise(stringstream& ss, bool clamp, int clampNeg, int clampPos, float valScale, float valOffset, float timeAnim)
{
	if (ss == NULL)/*should never be here*/
	{
		return;
	}
	ss << "noise = s.create(NoiseField) \n";
	ss << "noise.posScale = vec3(45) \n";
	ss << "noise.clamp = " << ((clamp)?"True":"False") << " \n";
	ss << "noise.clampNeg = " << clampNeg << " \n";
	ss << "noise.clampPos = " << clampPos << " \n";
	ss << "noise.valScale = " << valScale << " \n";
	ss << "noise.valOffset = " << valOffset << " \n";
	ss << "noise.timeAnim = " << timeAnim << " \n";
}

static void manta_advect_SemiLagr(stringstream& ss, char *indent, char *flags, char *vel, char *grid, int order)
{
	if((order <=1) || (indent == NULL) || (flags == NULL) || (vel == NULL) || (grid == NULL))
	{return;}
	ss << indent << "advectSemiLagrange(flags=" << flags << ", vel=" << vel \
	<< ", grid=" << grid << ", order=" << order << ") \n"; 
}

static void generate_manta_sim_file(Scene *scene, SmokeModifierData *smd)
{
	/*for now, simpleplume file creation
	*create python file with 2-spaces indentation*/
	
	FLUID_3D *fluid = smd->domain->fluid;
	
	ofstream manta_setup_file;
	manta_setup_file.open("manta_scene.py", std::fstream::trunc);
	stringstream ss; /*setup contents*/
	
	/*header*/
	ss << "from manta import * \n";

/*Data Declaration*/
	/*Solver Resolution*/
	ss << "res = " << smd->domain->maxres << " \n";
		/*Z axis in Blender = Y axis in Mantaflow*/
	ss << "gs = vec3(" << fluid->xRes() << ", " << fluid->zRes() << ", " << fluid->yRes() << ")" << " \n";
	ss << "s = Solver(name = 'main', gridSize = gs) \n";
	ss << "s.timestep = " << smd->domain->time_scale << " \n";
	
/*Grids setup*/
/*For now, only one grid of each kind is needed*/
	ss << "flags = s.create(FlagGrid) \n";/*must always be present*/
	ss << "vel = s.create(MACGrid) \n";
	ss << "density = s.create(RealGrid) \n";/*smoke simulation*/
	ss << "pressure = s.create(RealGrid) \n";/*must always be present*/

/*Noise Field*/
	manta_gen_noise(ss, true, 0, 1, 1, 0.75, 0.2);
/*	ss << "noise = s.create(NoiseField) \n");
	ss << "noise.posScale = vec3(45) \n");
	ss << "noise.clamp = True \n");
	ss << "noise.clampNeg = 0 \n");
	ss << "noise.clampPos = 1 \n");
	ss << "noise.valScale = 1 \n");
	ss << "noise.valOffset = 0.75 \n");
	ss << "noise.timeAnim = 0.2 \n");
*/
	
/*Flow setup*/
	ss << "flags.initDomain() \n";
	ss << "flags.fillGrid() \n";

/*GUI for debugging purposes*/
	ss << "if (GUI):\n  gui = Gui()\n  gui.show() \n";

/*Inflow source - for now, using mock sphere */
	ss << "source = s.create(Cylinder, center=gs*vec3(0.5,0.1,0.5), radius=res*0.14, z=gs*vec3(0, 0.02, 0)) \n";
	
/*Flow solving stepsv, main loop*/
	ss << "for t in xrange(" << scene->r.sfra << ", " << scene->r.efra << "): \n";
	ss << "  densityInflow(flags=flags, density=density, noise=noise, shape=source, scale=1, sigma=0.5) \n"	;
	manta_advect_SemiLagr(ss, "  ", "flags", "vel", "density", 2);
	manta_advect_SemiLagr(ss, "  ", "flags", "vel", "vel", 2);
/*	ss << "  advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2) \n");*/
/*	ss << "  advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2) \n");*/
	ss << "  setWallBcs(flags=flags, vel=vel) \n";
	ss << "  addBuoyancy(density=density, vel=vel, gravity=vec3(0,-6e-4,0), flags=flags) \n";
	ss << "  solvePressure(flags=flags, vel=vel, pressure=pressure, useResNorm=True, openBound='" << ((smd->domain->border_collisions == 2)?"N":"Y") << "') \n";/*2:closed border*/
	ss << "  setWallBcs(flags=flags, vel=vel) \n";

/*Saving output*/
	ss << "  density.save('den%04d.uni' % t) \n";
	ss << "  s.step()\n";
	ss << " \n";
	
	manta_setup_file << ss.rdbuf();
	manta_setup_file.close();		
}

#endif /* MANTA_H */

