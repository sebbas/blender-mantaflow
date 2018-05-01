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

/** \file mantaflow/intern/FLUID.h
 *  \ingroup mantaflow
 */

#ifndef FLUID_A_H
#define FLUID_A_H

#include <string>
#include <vector>
#include <atomic>

struct FLUID {
public:
	FLUID(int *res, struct SmokeModifierData *smd);
	FLUID() {};
	virtual ~FLUID();
	
	// Mirroring Mantaflow structures for particle data
	typedef struct pData { float pos[3]; int flag; } pData;
	typedef struct pVel { float pos[3]; } pVel;

	// Manta step, handling everything
	void step(struct SmokeModifierData *smd, int startFrame);
	
	// Grid initialization functions
	void initHeat(struct SmokeModifierData *smd);
	void initFire(struct SmokeModifierData *smd);
	void initColors(struct SmokeModifierData *smd);
	void initFireHigh(struct SmokeModifierData *smd);
	void initColorsHigh(struct SmokeModifierData *smd);
	void initLiquid(SmokeModifierData *smd);
	void initLiquidMesh(SmokeModifierData *smd);
	void initObstacle(SmokeModifierData *smd);
	void initGuiding(SmokeModifierData *smd);
	void initInVelocity(SmokeModifierData *smd);
	void initSndParts(SmokeModifierData *smd);
	void initLiquidSndParts(SmokeModifierData *smd);

	// Pointer transfer Mantaflow -> Blender
	void updatePointers();
	void updatePointersHigh();

	// Bake
	int readCache(SmokeModifierData *smd, int framenr);
	int writeCache(SmokeModifierData *smd, int framenr);
	int bakeData(SmokeModifierData *smd, int framenr);
	int bakeNoise(SmokeModifierData *smd, int framenr);
	int bakeMesh(SmokeModifierData *smd, int framenr);
	int bakeParticles(SmokeModifierData *smd, int framenr);
	void updateVariablesLow(SmokeModifierData *smd);
	void updateVariablesHigh(SmokeModifierData *smd);

	// IO for Mantaflow scene script
	void exportSmokeScript(struct SmokeModifierData *smd);
	void exportSmokeData(struct SmokeModifierData *smd);
	void exportLiquidScript(struct SmokeModifierData *smd);
	void exportLiquidData(struct SmokeModifierData *smd);

	// Write files for fluid
	void saveFluidObstacleData(char *pathname);
	void saveFluidGuidingData(char *pathname);
	void saveFluidInvelData(char *pathname);
	void saveFluidSndPartsData(char *pathname);

	// Write files for smoke
	void saveSmokeData(char *pathname);
	void saveSmokeDataHigh(char *pathname);

	// Write files for liquids
	void saveMesh(char *filename, int framenr);
	void saveLiquidData(char *pathname);
	void saveLiquidDataHigh(char *pathname);

	// Write files for particles
	void saveParticles(char* filename);
	void saveParticleVelocities(char* filename);

	// Load files for liquids
	void loadLiquidData(char *pathname);
	void loadLiquidDataHigh(char *pathname);

	// Smoke getters
	inline size_t getTotalCells() { return mTotalCells; }
	inline size_t getTotalCellsHigh() { return mTotalCellsHigh; }
	inline bool usingHighRes() { return mUsingNoise; }
	inline int getResX() { return mResX; }
	inline int getResY() { return mResY; }
	inline int getResZ() { return mResZ; }
	inline int getParticleResX() { return mResXParticle; }
	inline int getParticleResY() { return mResYParticle; }
	inline int getParticleResZ() { return mResZParticle; }
	inline int getMeshResX() { return mResXMesh; }
	inline int getMeshResY() { return mResYMesh; }
	inline int getMeshResZ() { return mResZMesh; }
	inline int getResXHigh() { return mResXNoise; }
	inline int getResYHigh() { return mResYNoise; }
	inline int getResZHigh() { return mResZNoise; }
	inline int getMeshUpres() { return mUpresMesh; }
	inline int getParticleUpres() { return mUpresParticle; }
	
	inline float* getDensity() { return mDensity; }
	inline float* getDensityIn() { return mDensityIn; }
	inline float* getHeat() { return mHeat; }
	inline float* getHeatIn() { return mHeatIn; }
	inline float* getVelocityX() { return mVelocityX; }
	inline float* getVelocityY() { return mVelocityY; }
	inline float* getVelocityZ() { return mVelocityZ; }
	inline float* getObVelocityX() { return mObVelocityX; }
	inline float* getObVelocityY() { return mObVelocityY; }
	inline float* getObVelocityZ() { return mObVelocityZ; }
	inline float* getGuideVelocityX() { return mGuideVelocityX; }
	inline float* getGuideVelocityY() { return mGuideVelocityY; }
	inline float* getGuideVelocityZ() { return mGuideVelocityZ; }
	inline float* getInVelocityX() { return mInVelocityX; }
	inline float* getInVelocityY() { return mInVelocityY; }
	inline float* getInVelocityZ() { return mInVelocityZ; }
	inline float* getForceX() { return mForceX; }
	inline float* getForceY() { return mForceY; }
	inline float* getForceZ() { return mForceZ; }
	inline int* getObstacle() { return mObstacle; }
	inline int* getNumObstacle() { return mNumObstacle; }
	inline int* getNumGuide()    { return mNumGuide; }
	inline float* getFlame() { return mFlame; }
	inline float* getFuel()  { return mFuel; }
	inline float* getFuelIn()  { return mFuelIn; }
	inline float* getReact() { return mReact; }
	inline float* getColorR() { return mColorR; }
	inline float* getColorG() { return mColorG; }
	inline float* getColorB() { return mColorB; }
	inline float* getColorRIn() { return mColorRIn; }
	inline float* getColorGIn() { return mColorGIn; }
	inline float* getColorBIn() { return mColorBIn; }
	inline float* getEmissionIn() { return mEmissionIn; }
	inline float* getShadow() { return mShadow; }
	inline int* getFlowType() { return mFlowType; }
	inline int* getNumFlow()  { return mNumFlow; }

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
	
	inline float* getPhiIn()      { return mPhiIn; }
	inline float* getPhiObsIn()   { return mPhiObsIn; }
	inline float* getPhiGuideIn() { return mPhiGuideIn; }
	inline float* getPhiOutIn()   { return mPhiOutIn; }
	inline float* getPhi()        { return mPhi; }

	static std::atomic<bool> mantaInitialized;
	static std::atomic<int> solverID;
	static int with_debug; // on or off (1 or 0), also sets manta debug level
	
	// Liquid getters
	inline int getNumVertices()  { return mNumVertices; }
	inline int getNumNormals()   { return mNumNormals; }
	inline int getNumTriangles() { return mNumTriangles; }
	
	inline float getVertexXAt(int i) { return mVerticesX[i]; }
	inline float getVertexYAt(int i) { return mVerticesY[i]; }
	inline float getVertexZAt(int i) { return mVerticesZ[i]; }

	inline float getNormalXAt(int i) { return mNormalsX[i]; }
	inline float getNormalYAt(int i) { return mNormalsY[i]; }
	inline float getNormalZAt(int i) { return mNormalsZ[i]; }

	inline int getTriangleXAt(int i) { return mTrianglesX[i]; }
	inline int getTriangleYAt(int i) { return mTrianglesY[i]; }
	inline int getTriangleZAt(int i) { return mTrianglesZ[i]; }

	// Particle getters
	inline int getFlipParticleFlagAt(int i) { return (mFlipParticleData) ? ((std::vector<pData>*) mFlipParticleData)->at(i).flag : 0; }
	inline int getSndParticleFlagAt(int i) { return (mSndParticleData) ? ((std::vector<pData>*) mSndParticleData)->at(i).flag : 0; }

	inline float getFlipParticlePositionXAt(int i) { return (mFlipParticleData && !mFlipParticleData->empty()) ? mFlipParticleData->at(i).pos[0] : 0.f; }
	inline float getFlipParticlePositionYAt(int i) { return (mFlipParticleData && !mFlipParticleData->empty()) ? mFlipParticleData->at(i).pos[1] : 0.f; }
	inline float getFlipParticlePositionZAt(int i) { return (mFlipParticleData && !mFlipParticleData->empty()) ? mFlipParticleData->at(i).pos[2] : 0.f; }

	inline float getSndParticlePositionXAt(int i) { return (mSndParticleData && !mSndParticleData->empty()) ? mSndParticleData->at(i).pos[0] : 0.f; }
	inline float getSndParticlePositionYAt(int i) { return (mSndParticleData && !mSndParticleData->empty()) ? mSndParticleData->at(i).pos[1] : 0.f; }
	inline float getSndParticlePositionZAt(int i) { return (mSndParticleData && !mSndParticleData->empty()) ? mSndParticleData->at(i).pos[2] : 0.f; }

	inline float getFlipParticleVelocityXAt(int i) { return (mFlipParticleVelocity && !mFlipParticleVelocity->empty()) ? mFlipParticleVelocity->at(i).pos[0] : 0.f; }
	inline float getFlipParticleVelocityYAt(int i) { return (mFlipParticleVelocity && !mFlipParticleVelocity->empty()) ? mFlipParticleVelocity->at(i).pos[1] : 0.f; }
	inline float getFlipParticleVelocityZAt(int i) { return (mFlipParticleVelocity && !mFlipParticleVelocity->empty()) ? mFlipParticleVelocity->at(i).pos[2] : 0.f; }

	inline float getSndParticleVelocityXAt(int i) { return (mSndParticleVelocity && !mSndParticleVelocity->empty()) ? mSndParticleVelocity->at(i).pos[0] : 0.f; }
	inline float getSndParticleVelocityYAt(int i) { return (mSndParticleVelocity && !mSndParticleVelocity->empty()) ? mSndParticleVelocity->at(i).pos[1] : 0.f; }
	inline float getSndParticleVelocityZAt(int i) { return (mSndParticleVelocity && !mSndParticleVelocity->empty()) ? mSndParticleVelocity->at(i).pos[2] : 0.f; }

	inline float* getFlipParticleData() { return (mFlipParticleData && !mFlipParticleData->empty()) ? (float*) &mFlipParticleData->front() : NULL; }
	inline float* getSndParticleData()  { return (mSndParticleData && !mSndParticleData->empty()) ? (float*) &mSndParticleData->front() : NULL; }

	inline float* getFlipParticleVelocity() { return (mFlipParticleVelocity && !mFlipParticleVelocity->empty()) ? (float*) &mFlipParticleVelocity->front() : NULL; }
	inline float* getSndParticleVelocity()  { return (mSndParticleVelocity && !mSndParticleVelocity->empty()) ? (float*) &mSndParticleVelocity->front() : NULL; }
	inline float* getSndParticleLife()      { return (mSndParticleLife && !mSndParticleLife->empty()) ? (float*) &mSndParticleLife->front() : NULL; }

	inline int getNumFlipParticles() { return (mFlipParticleData && !mFlipParticleData->empty()) ? mFlipParticleData->size() : 0; }
	inline int getNumSndParticles() { return (mSndParticleData && !mSndParticleData->empty()) ? mSndParticleData->size() : 0; }

	void updateMeshData(const char* filename);
	void updateParticleData(const char* filename, bool isSecondary);

	void setFlipParticleData(float* buffer, int numParts);
	void setSndParticleData(float* buffer, int numParts);

	void setFlipParticleVelocity(float* buffer, int numParts);
	void setSndParticleVelocity(float* buffer, int numParts);
	void setSndParticleLife(float* buffer, int numParts);

	// Direct access to solver time attributes
	int getFrame();
	float getTimestep();
	void adaptTimestep();

private:
	// simulation constants
	size_t mTotalCells;
	size_t mTotalCellsHigh;
	size_t mTotalCellsMesh;
	size_t mTotalCellsParticles;

	int mCurrentID;
	
	bool mUsingHeat;
	bool mUsingColors;
	bool mUsingFire;
	bool mUsingObstacle;
	bool mUsingGuiding;
	bool mUsingInvel;
	bool mUsingNoise;
	bool mUsingMesh;
	bool mUsingLiquid;
	bool mUsingSmoke;
	bool mUsingDrops;
	bool mUsingBubbles;
	bool mUsingFloats;
	bool mUsingTracers;

	int mResX;
	int mResY;
	int mResZ;
	int mMaxRes;

	int mResXNoise;
	int mResYNoise;
	int mResZNoise;
	int mResXMesh;
	int mResYMesh;
	int mResZMesh;
	int mResXParticle;
	int mResYParticle;
	int mResZParticle;

	int mUpresMesh;
	int mUpresParticle;

	float mTempAmb; /* ambient temperature */
	float mConstantScaling;

	// Smoke grids
	float* mDensity;
	float* mDensityIn;
	float* mHeat;
	float* mHeatIn;
	float* mVelocityX;
	float* mVelocityY;
	float* mVelocityZ;
	float* mObVelocityX;
	float* mObVelocityY;
	float* mObVelocityZ;
	float* mGuideVelocityX;
	float* mGuideVelocityY;
	float* mGuideVelocityZ;
	float* mInVelocityX;
	float* mInVelocityY;
	float* mInVelocityZ;
	float* mForceX;
	float* mForceY;
	float* mForceZ;
	int* mObstacle;
	int* mNumObstacle;
	int* mNumGuide;
	float* mFlame;
	float* mFuel;
	float* mFuelIn;
	float* mReact;
	float* mColorR;
	float* mColorG;
	float* mColorB;
	float* mColorRIn;
	float* mColorGIn;
	float* mColorBIn;
	float* mEmissionIn;
	float* mShadow;
	int* mFlowType;
	int* mNumFlow;
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
	float* mPhiIn;
	float* mPhiObsIn;
	float* mPhiGuideIn;
	float* mPhiOutIn;
	float* mPhi;

	// Mesh fields for liquid surface
	int mNumVertices;
	int mNumNormals;
	int mNumTriangles;
	std::vector<float> mVerticesX;
	std::vector<float> mVerticesY;
	std::vector<float> mVerticesZ;
	std::vector<float> mNormalsX;
	std::vector<float> mNormalsY;
	std::vector<float> mNormalsZ;
	std::vector<int> mTrianglesX;
	std::vector<int> mTrianglesY;
	std::vector<int> mTrianglesZ;
	
	// Particle fields
	std::vector<pData>* mFlipParticleData;
	std::vector<pVel>* mFlipParticleVelocity;

	std::vector<pData>* mSndParticleData;
	std::vector<pVel>* mSndParticleVelocity;
	std::vector<float>* mSndParticleLife;

	void initDomain(struct SmokeModifierData *smd);
	void initNoise(struct SmokeModifierData *smd);
	void initMesh(struct SmokeModifierData *smd);
	void initSmoke(struct SmokeModifierData *smd);
	void initSmokeNoise(struct SmokeModifierData *smd);
	void initializeMantaflow();
	void terminateMantaflow();
	void runPythonString(std::vector<std::string> commands);
	std::string getRealValue(const std::string& varName, SmokeModifierData *smd);
	std::string parseLine(const std::string& line, SmokeModifierData *smd);
	std::string parseScript(const std::string& setup_string, SmokeModifierData *smd);
	void updateMeshDataFromBobj(const char* filename);
	void updateMeshDataFromObj(const char* filename);

};

#endif
