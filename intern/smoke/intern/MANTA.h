#ifndef MANTA_H
#define MANTA_H
#include "FLUID_3D.h"
#include "zlib.h"

extern "C" bool manta_check_grid_size(struct FLUID_3D *fluid, int dimX, int dimY, int dimZ)
{
	if (!(dimX == fluid->xRes() && dimY == fluid->yRes() && dimZ == fluid->zRes())) {
		for (int cnt(0); cnt < fluid->_totalCells; cnt++)
			fluid->_density[cnt] = 0.0f;
		return false;
	}
	return true;
}

extern "C" void read_mantaflow_sim(struct FLUID_3D *fluid, char* name)
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

#endif /* MANTA_H */

