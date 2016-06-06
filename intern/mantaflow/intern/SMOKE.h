/*
 * ***** BEGIN GPL LICENSE BLOCK *****
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version. 
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * The Original Code is Copyright (C) 2016 Blender Foundation.
 * All rights reserved.
 *
 * Contributor(s): Sebastian Barschkis (sebbas)
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file mantaflow/intern/SMOKE.h
 *  \ingroup mantaflow
 */

#ifndef SMOKE_A_H
#define SMOKE_A_H

#include <string>
#include <vector>

#include "Python.h"

struct SMOKE {
public:
	SMOKE(int *res, struct SmokeModifierData *smd);
	SMOKE() {};
	virtual ~SMOKE();
	
	// Manta step, handling everything
	void step(struct SmokeModifierData *smd);
	
	// Grid initialization functions
	void initHeat(struct SmokeModifierData *smd);
	void initFire(struct SmokeModifierData *smd);
	void initColors(struct SmokeModifierData *smd);
	void initFireHigh(struct SmokeModifierData *smd);
	void initColorsHigh(struct SmokeModifierData *smd);
	void initLiquid(SmokeModifierData *smd);
	
	// Pointer transfer Mantaflow -> Blender
	void updatePointers(struct SmokeModifierData *smd);
	void updatePointersHigh(struct SmokeModifierData *smd);

	// IO for Mantaflow scene script
	void exportScript(struct SmokeModifierData *smd);
	void exportGrids(struct SmokeModifierData *smd);
	
	// Getters
	inline size_t getTotalCells() { return mTotalCells; }
	inline size_t getTotalCellsHigh() { return mTotalCellsHigh; }
	inline bool usingHighRes() { return mUsingHighRes; }
	inline int getResX() { return mResX; }
	inline int getResY() { return mResY; }
	inline int getResZ() { return mResZ; }
	inline int getResXHigh() { return mResXHigh; }
	inline int getResYHigh() { return mResYHigh; }
	inline int getResZHigh() { return mResZHigh; }
	
	inline float* getDensity() { return mDensity; }
	inline float* getHeat() { return mHeat; }
	inline float* getVelocityX() { return mVelocityX; }
	inline float* getVelocityY() { return mVelocityY; }
	inline float* getVelocityZ() { return mVelocityZ; }
	inline float* getObVelocityX() { return mObVelocityX; }
	inline float* getObVelocityY() { return mObVelocityY; }
	inline float* getObVelocityZ() { return mObVelocityZ; }
	inline float* getForceX() { return mForceX; }
	inline float* getForceY() { return mForceY; }
	inline float* getForceZ() { return mForceZ; }
	inline unsigned char* getObstacles() { return mObstacles; }
	inline unsigned char* getObstaclesAnim() { return mObstaclesAnim; }
	inline float* getFlame() { return mFlame; }
	inline float* getFuel() { return mFuel; }
	inline float* getReact() { return mReact; }
	inline float* getColorR() { return mColorR; }
	inline float* getColorG() { return mColorG; }
	inline float* getColorB() { return mColorB; }
	inline float* getDensityInflow() { return mDensityInflow; }
	inline float* getFuelInflow() { return mFuelInflow; }
	inline int* getMantaFlags() { return mMantaFlags; }

	inline float* getDensityHigh() { return mDensityHigh; }
	inline float* getFlameHigh() { return mFlameHigh; }
	inline float* getFuelHigh() { return mFuelHigh; }
	inline float* getReactHigh() { return mReactHigh; }
	inline float* getColorRHigh() { return mColorRHigh; }
	inline float* getColorGHigh() { return mColorGHigh; }
	inline float* getColorBHigh() { return mColorBHigh; }
	inline float* getTextureU() { return mTextureU; }
	inline float* getTextureV() { return mTextureV; }
	inline float* getTextureW() { return mTextureW; }
	inline float* getTextureU2() { return mTextureU2; }
	inline float* getTextureV2() { return mTextureV2; }
	inline float* getTextureW2() { return mTextureW2; }
	
	inline float* getPhi() { return mPhi; }

	static bool mantaInitialized;

private:
	// simulation constants
	size_t mTotalCells;
	size_t mTotalCellsHigh;
	
	bool mUsingHeat;
	bool mUsingColors;
	bool mUsingFire;
	bool mUsingHighRes;
	bool mUsingLiquid;
	bool mUsingSmoke;
	
	int mResX;
	int mResY;
	int mResZ;
	
	int mResXHigh;
	int mResYHigh;
	int mResZHigh;
	
	float mTempAmb; /* ambient temperature */
	float mConstantScaling;
	std::vector<std::string> mCommands;

	// Smoke grids low res
	float* mDensity;
	float* mHeat;
	float* mVelocityX;
	float* mVelocityY;
	float* mVelocityZ;
	float* mObVelocityX;
	float* mObVelocityY;
	float* mObVelocityZ;
	float* mForceX;
	float* mForceY;
	float* mForceZ;
	unsigned char* mObstacles; /* only used (useful) for static obstacles like domain boundaries */
	unsigned char* mObstaclesAnim;
	float *mFlame;
	float *mFuel;
	float *mReact;
	float *mColorR;
	float *mColorG;
	float *mColorB;
	float* mDensityInflow;
	float* mFuelInflow;
	int* mMantaFlags;

	// Smoke grids high res
	float* mDensityHigh;
	float* mFlameHigh;
	float* mFuelHigh;
	float* mReactHigh;
	float* mColorRHigh;
	float* mColorGHigh;
	float* mColorBHigh;
	float* mTextureU;
	float* mTextureV;
	float* mTextureW;
	float* mTextureU2;
	float* mTextureV2;
	float* mTextureW2;
	
	// Liquids
	float* mPhi;
	
	void initDomain(struct SmokeModifierData *smd);
	void initDomainHigh(struct SmokeModifierData *smd);
	void initSmoke(struct SmokeModifierData *smd);
	void initSmokeHigh(struct SmokeModifierData *smd);
	void startMantaflow();
	void runPythonString(std::vector<std::string> commands);
	std::string getRealValue(const std::string& varName, SmokeModifierData *smd);
	std::string parseLine(const std::string& line, SmokeModifierData *smd);
	std::string parseScript(const std::string& setup_string, SmokeModifierData *smd);
	void* getGridPointer(std::string gridName, std::string solverName);
};

#endif