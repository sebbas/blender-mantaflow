#ifndef MANTA_H
#define MANTA_H
#include "FLUID_3D.h"
#include "zlib.h"
#include "../../../source/blender/makesdna/DNA_scene_types.h"
#include "../../../source/blender/makesdna/DNA_modifier_types.h"
#include "../../../source/blender/makesdna/DNA_smoke_types.h"
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

static void manta_gen_noise(FILE *f, bool clamp, int clampNeg, int clampPos, float valScale, float valOffset, float timeAnim)
{
	if (f == NULL)/*should never be here*/
	{
		return;
	}
	fprintf(f, "  noise = s.create(NoiseField) \n");
	fprintf(f, "  noise.posScale = vec3(45) \n");
	fprintf(f, "  noise.clamp = %s \n", (clamp)?"True":"False");
	fprintf(f, "  noise.clampNeg = %d \n", clampNeg);
	fprintf(f, "  noise.clampPos = %d \n", clampPos);
	fprintf(f, "  noise.valScale = %f \n", valScale);
	fprintf(f, "  noise.valOffset = %f \n", valOffset);
	fprintf(f, "  noise.timeAnim = %f \n", timeAnim);
}

static void generate_manta_sim_file(Scene *scene, SmokeModifierData *smd)
{
	/*for now, simpleplume file creation
	*create python file with 2-spaces indentation*/
	
	FLUID_3D *fluid = smd->domain->fluid;
	FILE *f = fopen("manta_scene.py", "w");
	if (f == NULL)
	{
		exit(1);
	}
	/*header*/
	fprintf(f, "from manta import * \n");

/*Data Declaration*/
	/*Solver Resolution*/
	fprintf(f, "res = %d \n", smd->domain->maxres);
		/*Z axis in Blender = Y axis in Mantaflow*/
	fprintf(f, "gs = vec3(%d, %d, %d) \n", fluid->xRes(), fluid->zRes(), fluid->yRes());
	fprintf(f, "s = Solver(name = 'main', gridSize = gs) \n");
	fprintf(f, "s.timestep = %f \n", smd->domain->time_scale);
	
/*Grids setup*/
/*For now, only one grid of each kind is needed*/
	fprintf(f, "flags = s.create(FlagGrid) \n");/*must always be present*/
	fprintf(f, "vel = s.create(MACGrid) \n");
	fprintf(f, "density = s.create(RealGrid) \n");/*smoke simulation*/
	fprintf(f, "pressure = s.create(RealGrid) \n");/*must always be present*/

/*Noise Field*/
	manta_gen_noise(f, true, 0, 1, 1, 0.75, 0.2);
	
/*Flow setup*/
	fprintf(f, "flags.initDomain() \n");
	fprintf(f, "flags.fillGrid() \n");

/*GUI for debugging purposes*/
	fprintf(f, "if (GUI):\n  gui = Gui()\n  gui.show() \n");

/*Inflow source - for now, using mock sphere */
	fprintf(f, "source = s.create(Cylinder, center=gs*vec3(0.5,0.1,0.5), radius=res*0.14, z=gs*vec3(0, 0.02, 0)) \n");
	
/*Flow solving stepsv, main loop*/
	fprintf(f, "for t in xrange(%d, %d): \n", scene->r.sfra, scene->r.efra);
	fprintf(f, "  densityInflow(flags=flags, density=density, noise=noise, shape=source, scale=1, sigma=0.5) \n");
	fprintf(f, "  advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2) \n");
	fprintf(f, "  advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2) \n");
	fprintf(f, "  setWallBcs(flags=flags, vel=vel) \n");
	fprintf(f, "  addBuoyancy(density=density, vel=vel, gravity=vec3(0,-6e-4,0), flags=flags) \n");
	fprintf(f, "  solvePressure(flags=flags, vel=vel, pressure=pressure, useResNorm=True, openBound='%s') \n",(smd->domain->border_collisions == 2)?"N":"Y");/*2:closed border*/
	fprintf(f, "  setWallBcs(flags=flags, vel=vel) \n");

/*Saving output*/
	char format_str[] = "  density.save('den%04d.uni' % t) \n";
	fwrite(format_str, 1, sizeof(format_str)-1, f);
	fprintf(f, "  s.step()\n");
	fprintf(f, " \n");
	
	fclose(f);		
}

#endif /* MANTA_H */

