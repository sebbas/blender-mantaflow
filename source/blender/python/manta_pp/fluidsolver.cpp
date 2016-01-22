




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/Users/user/Developer/Xcode Projects/mantaflowDevelop/mantaflowgit/source/fluidsolver.cpp"
/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Main class for the fluid solver
 *
 ******************************************************************************/

#include "fluidsolver.h"
#include "grid.h"
#include <sstream>
#include <fstream>

using namespace std;
namespace Manta {

//******************************************************************************
// Gridstorage-related members

template<class T>
void FluidSolver::GridStorage<T>::free() {
	if (used != 0)
		errMsg("can't clean grid cache, some grids are still in use");
	for(size_t i = 0; i<grids.size(); i++)
		delete[] grids[i];
	grids.clear();
}
template<class T>
T* FluidSolver::GridStorage<T>::get(Vec3i size) {
	if ((int)grids.size() <= used) {
		grids.push_back(new T[size.x * size.y * size.z]);
	}
	if (used > 200)
		errMsg("too many temp grids used -- are they released properly ?");
	return grids[used++];
}
template<class T>
void FluidSolver::GridStorage<T>::release(T* ptr) {
	// rewrite pointer, as it may have changed due to swap operations
	used--;
	if (used < 0)
		errMsg("temp grid inconsistency");
	grids[used] = ptr;
}

template<> int* FluidSolver::getGridPointer<int>() {
	return mGridsInt.get(mGridSize);    
}
template<> Real* FluidSolver::getGridPointer<Real>() {
	return mGridsReal.get(mGridSize);    
}
template<> Vec3* FluidSolver::getGridPointer<Vec3>() {
	return mGridsVec.get(mGridSize);    
}
template<> void FluidSolver::freeGridPointer<int>(int *ptr) {
	mGridsInt.release(ptr);
}
template<> void FluidSolver::freeGridPointer<Real>(Real* ptr) {
	mGridsReal.release(ptr);
}
template<> void FluidSolver::freeGridPointer<Vec3>(Vec3* ptr) {
	mGridsVec.release(ptr);
}

//******************************************************************************
// FluidSolver members

FluidSolver::FluidSolver(Vec3i gridsize, int dim)
	: PbClass(this), mDt(1.0), mTimeTotal(0.), mFrame(0), 
	  mCflCond(1000), mDtMin(1.), mDtMax(1.), mFrameLength(1.),
	  mGridSize(gridsize), mDim(dim) , mTimePerFrame(0.), mLockDt(false), mAdaptDt(true)
{    
	assertMsg(dim==2 || dim==3, "Can only create 2D and 3D solvers");
	assertMsg(dim!=2 || gridsize.z == 1, "Trying to create 2D solver with size.z != 1");
}

FluidSolver::~FluidSolver() {
	mGridsInt.free();
	mGridsReal.free();
	mGridsVec.free();
}

PbClass* FluidSolver::create(PbType t, PbTypeVec T, const string& name) {        
	_args.add("nocheck",true);
	if (t.str() == "")
		errMsg("Need to specify object type. Use e.g. Solver.create(FlagGrid, ...) or Solver.create(type=FlagGrid, ...)");
	
	PbClass* ret = PbClass::createPyObject(t.str() + T.str(), name, _args, this);
	return ret;
}

void FluidSolver::step() {
	// update simulation time
	if(!mAdaptDt) {
		mTimeTotal += mDt;
		mFrame++;
	} else {
		// adaptive time stepping on (use eps to prevent roundoff errors)
		mTimePerFrame += mDt;
		if( (mTimePerFrame+VECTOR_EPSILON) >mFrameLength) {
			mFrame++;

			// re-calc total time, prevent drift...
			mTimeTotal = (double)mFrame * mFrameLength;
			mTimePerFrame = 0.;
			mLockDt = false;
		}
	}

	updateQtGui(true, mFrame,mTimeTotal, "FluidSolver::step");
}

void FluidSolver::printMemInfo() {
	std::ostringstream msg;
	msg << "Allocated grids: int " << mGridsInt.used  <<"/"<< mGridsInt.grids.size()  <<", ";
	msg <<                  "real "<< mGridsReal.used <<"/"<< mGridsReal.grids.size() <<", ";
	msg <<                  "vec3 "<< mGridsVec.used  <<"/"<< mGridsVec.grids.size()  <<". ";
	printf("%s\n", msg.str().c_str() );
}

std::string printBuildInfo() {
	string infoString = buildInfoString();
	debMsg( "Build info: "<<infoString.c_str()<<" ",1);
	return infoString;
} static PyObject* _W_0 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "printBuildInfo" ); PyObject *_retval = 0; { ArgLocker _lock;   _retval = toPy(printBuildInfo());  _args.check(); } pbFinalizePlugin(parent,"printBuildInfo" ); return _retval; } catch(std::exception& e) { pbSetError("printBuildInfo",e.what()); return 0; } } static const Pb::Register _RP_printBuildInfo ("","printBuildInfo",_W_0); 

void setDebugLevel(int level=1) {
	gDebugLevel = level; 
} static PyObject* _W_1 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "setDebugLevel" ); PyObject *_retval = 0; { ArgLocker _lock; int level = _args.getOpt<int >("level",0,1,&_lock);   _retval = getPyNone(); setDebugLevel(level);  _args.check(); } pbFinalizePlugin(parent,"setDebugLevel" ); return _retval; } catch(std::exception& e) { pbSetError("setDebugLevel",e.what()); return 0; } } static const Pb::Register _RP_setDebugLevel ("","setDebugLevel",_W_1); 

void FluidSolver::adaptTimestep(Real maxVel)
{
	if (!mLockDt) {
		// calculate current timestep from maxvel, clamp range
		mDt = std::max( std::min( (Real)(mCflCond/(maxVel+1e-05)), mDtMax) , mDtMin );
		if( (mTimePerFrame+mDt*1.05) > mFrameLength ) {
			// within 5% of full step? add epsilon to prevent roundoff errors...
			mDt = ( mFrameLength - mTimePerFrame ) + 1e-04;
		}
		else if ( (mTimePerFrame+mDt + mDtMin) > mFrameLength || (mTimePerFrame+(mDt*1.25)) > mFrameLength ) {
			// avoid tiny timesteps and strongly varying ones, do 2 medium size ones if necessary...
			mDt = (mFrameLength-mTimePerFrame+ 1e-04)*0.5;
			mLockDt = true;
		}
	}
	debMsg( "Frame "<<mFrame<<" current max vel: "<<maxVel<<" , dt: "<<mDt<<", "<<mTimePerFrame<<"/"<<mFrameLength<<" lock:"<<mLockDt , 1);
	mAdaptDt = true;

	// sanity check
	assertMsg( (mDt > (mDtMin/2.) ) , "Invalid dt encountered! Shouldnt happen..." );
}

} // manta



