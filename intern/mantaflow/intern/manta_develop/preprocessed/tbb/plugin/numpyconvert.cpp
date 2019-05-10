




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




/******************************************************************************
*
* MantaFlow fluid solver framework
* Copyright 2017 Steffen Wiewel, Moritz Becher 
*
* This program is free software, distributed under the terms of the
* Apache License, Version 2.0 
* http://www.apache.org/licenses/LICENSE-2.0
*
* Plugins to convert mantaflow grids to/from numpy arrays, also support pdata fields  
# (only compiled if NUMPY is enabled)
*
******************************************************************************/

#include "manta.h"
#include "kernel.h"
#include "grid.h"
#include "particle.h"
#include "levelset.h"

using namespace std;
namespace Manta
{

//====================================================================================================
// Grid numpy conversion
//----------------------------------------------------------------------------------------------------

template<typename T>
void copyArrayToGridScalar(const PyArrayContainer source, T& target)
{
	target.setConst(0.0f);
	unsigned int uGridSize = target.getSizeX() * target.getSizeY() * target.getSizeZ();
	assertMsg(source.TotalSize == uGridSize, "The size of the numpy array doesn't match the size of the Grid!");
	
	NumpyTypes eDataType  = source.DataType; 

	switch (eDataType)
	{
		case NumpyTypes::N_FLOAT:
			FOR_IDX(target) { target(idx) = (reinterpret_cast<float*>(source.pData))[idx]; }
			break;
		case NumpyTypes::N_DOUBLE:
			FOR_IDX(target) { target(idx) = (reinterpret_cast<double*>(source.pData))[idx]; } 
			break;
		default:
			errMsg("unknown/unsupported type of Numpy array");
			return;
	}
}

template<typename T>
void copyGridToArrayScalar(const T& source, PyArrayContainer target)
{
	unsigned int uGridsize = source.getSizeX() * source.getSizeY() * source.getSizeZ();
	assertMsg(target.TotalSize == uGridsize, "The size of the numpy array doesn't match the size of the grid!");
	
	NumpyTypes eDataType = target.DataType;

	switch (eDataType)
	{
		case NumpyTypes::N_FLOAT:
			FOR_IDX(source) { reinterpret_cast<float*>(target.pData)[idx] = source(idx); }
			break;
		case NumpyTypes::N_DOUBLE:
			FOR_IDX(source) { reinterpret_cast<double*>(target.pData)[idx] = source(idx); }
			break;
		default:
			errMsg("unknown/unsupported type of Numpy array");
			break;
	}
}

template<typename T>
void copyArrayToGridVector(const PyArrayContainer source, T& target)
{
	unsigned int uSizeX = target.getSizeX();
	unsigned int uSizeY = target.getSizeY();
	unsigned int uSizeZ = target.getSizeZ();
	unsigned int uSizeW = 3u;
	
	assertMsg(source.TotalSize == uSizeX * uSizeY * uSizeZ * uSizeW, "The size of the numpy array doesn't match the size of the grid!");
	
	NumpyTypes eDataType = source.DataType;

	switch (eDataType)
	{
		case NumpyTypes::N_FLOAT:
			FOR_IDX(target) { for(int w = 0; w < 3; ++w) { target(idx)[w] = (reinterpret_cast<float*>(source.pData))[idx*3+w]; } }
			break;
		case NumpyTypes::N_DOUBLE:
			FOR_IDX(target) { for(int w = 0; w < 3; ++w) { target(idx)[w] = (reinterpret_cast<double*>(source.pData))[idx*3+w]; } }
			break;
		default:
			errMsg("unknown/unsupported type of Vec3 Numpy array");
			break;
	}
}

template<typename T>
void copyGridToArrayVector(const T& source, PyArrayContainer target)
{
	unsigned int uSizeX = source.getSizeX();
	unsigned int uSizeY = source.getSizeY();
	unsigned int uSizeZ = source.getSizeZ();
	unsigned int uSizeW = 3u;

	assertMsg(target.TotalSize == uSizeX * uSizeY * uSizeZ * uSizeW, "The size of the numpy array doesn't match the size of the grid!");
	
	NumpyTypes eDataType = target.DataType;
	
	switch (eDataType)
	{
		case NumpyTypes::N_FLOAT:
			FOR_IDX(source) { for(int w = 0; w < 3; ++w) { (reinterpret_cast<float*>(target.pData))[idx*3+w] = source(idx)[w]; } }
			break;
		case NumpyTypes::N_DOUBLE:
			FOR_IDX(source) { for(int w = 0; w < 3; ++w) { (reinterpret_cast<double*>(target.pData))[idx*3+w] = source(idx)[w]; } }
			break;
		default:
			errMsg("unknown/unsupported type of Vec3 Numpy array");
			break;
	}
}

//====================================================================================================
// Python interface
//----------------------------------------------------------------------------------------------------

void copyArrayToGridReal(const PyArrayContainer source, Grid<Real>& target) {
	copyArrayToGridScalar<Grid<Real>>(source, target);
} static PyObject* _W_0 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) {
 try {
 PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "copyArrayToGridReal" , !noTiming ); PyObject *_retval = 0; {
 ArgLocker _lock; const PyArrayContainer source = _args.get<PyArrayContainer >("source",0,&_lock); Grid<Real>& target = *_args.getPtr<Grid<Real> >("target",1,&_lock);   _retval = getPyNone(); copyArrayToGridReal(source,target);  _args.check(); }
 pbFinalizePlugin(parent,"copyArrayToGridReal", !noTiming ); return _retval; }
 catch(std::exception& e) {
 pbSetError("copyArrayToGridReal",e.what()); return 0; }
 }
 static const Pb::Register _RP_copyArrayToGridReal ("","copyArrayToGridReal",_W_0); 
extern "C" {
 void PbRegister_copyArrayToGridReal() {
 KEEP_UNUSED(_RP_copyArrayToGridReal); }
 }



void copyGridToArrayReal(const Grid<Real>& source, PyArrayContainer target) {
	copyGridToArrayScalar<Grid<Real>>(source, target);
} static PyObject* _W_1 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) {
 try {
 PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "copyGridToArrayReal" , !noTiming ); PyObject *_retval = 0; {
 ArgLocker _lock; const Grid<Real>& source = *_args.getPtr<Grid<Real> >("source",0,&_lock); PyArrayContainer target = _args.get<PyArrayContainer >("target",1,&_lock);   _retval = getPyNone(); copyGridToArrayReal(source,target);  _args.check(); }
 pbFinalizePlugin(parent,"copyGridToArrayReal", !noTiming ); return _retval; }
 catch(std::exception& e) {
 pbSetError("copyGridToArrayReal",e.what()); return 0; }
 }
 static const Pb::Register _RP_copyGridToArrayReal ("","copyGridToArrayReal",_W_1); 
extern "C" {
 void PbRegister_copyGridToArrayReal() {
 KEEP_UNUSED(_RP_copyGridToArrayReal); }
 }



void copyArrayToGridLevelset(const PyArrayContainer source, LevelsetGrid& target) {
	copyArrayToGridScalar<LevelsetGrid>(source, target);
} static PyObject* _W_2 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) {
 try {
 PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "copyArrayToGridLevelset" , !noTiming ); PyObject *_retval = 0; {
 ArgLocker _lock; const PyArrayContainer source = _args.get<PyArrayContainer >("source",0,&_lock); LevelsetGrid& target = *_args.getPtr<LevelsetGrid >("target",1,&_lock);   _retval = getPyNone(); copyArrayToGridLevelset(source,target);  _args.check(); }
 pbFinalizePlugin(parent,"copyArrayToGridLevelset", !noTiming ); return _retval; }
 catch(std::exception& e) {
 pbSetError("copyArrayToGridLevelset",e.what()); return 0; }
 }
 static const Pb::Register _RP_copyArrayToGridLevelset ("","copyArrayToGridLevelset",_W_2); 
extern "C" {
 void PbRegister_copyArrayToGridLevelset() {
 KEEP_UNUSED(_RP_copyArrayToGridLevelset); }
 }



void copyGridToArrayLevelset(const LevelsetGrid& source, PyArrayContainer target) {
	copyGridToArrayScalar<LevelsetGrid>(source, target);
} static PyObject* _W_3 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) {
 try {
 PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "copyGridToArrayLevelset" , !noTiming ); PyObject *_retval = 0; {
 ArgLocker _lock; const LevelsetGrid& source = *_args.getPtr<LevelsetGrid >("source",0,&_lock); PyArrayContainer target = _args.get<PyArrayContainer >("target",1,&_lock);   _retval = getPyNone(); copyGridToArrayLevelset(source,target);  _args.check(); }
 pbFinalizePlugin(parent,"copyGridToArrayLevelset", !noTiming ); return _retval; }
 catch(std::exception& e) {
 pbSetError("copyGridToArrayLevelset",e.what()); return 0; }
 }
 static const Pb::Register _RP_copyGridToArrayLevelset ("","copyGridToArrayLevelset",_W_3); 
extern "C" {
 void PbRegister_copyGridToArrayLevelset() {
 KEEP_UNUSED(_RP_copyGridToArrayLevelset); }
 }



void copyArrayToGridVec3(const PyArrayContainer source, Grid<Vec3>& target) {
	copyArrayToGridVector<Grid<Vec3>>(source, target);
} static PyObject* _W_4 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) {
 try {
 PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "copyArrayToGridVec3" , !noTiming ); PyObject *_retval = 0; {
 ArgLocker _lock; const PyArrayContainer source = _args.get<PyArrayContainer >("source",0,&_lock); Grid<Vec3>& target = *_args.getPtr<Grid<Vec3> >("target",1,&_lock);   _retval = getPyNone(); copyArrayToGridVec3(source,target);  _args.check(); }
 pbFinalizePlugin(parent,"copyArrayToGridVec3", !noTiming ); return _retval; }
 catch(std::exception& e) {
 pbSetError("copyArrayToGridVec3",e.what()); return 0; }
 }
 static const Pb::Register _RP_copyArrayToGridVec3 ("","copyArrayToGridVec3",_W_4); 
extern "C" {
 void PbRegister_copyArrayToGridVec3() {
 KEEP_UNUSED(_RP_copyArrayToGridVec3); }
 }



void copyGridToArrayVec3(const Grid<Vec3>& source, PyArrayContainer target) {
	copyGridToArrayVector<Grid<Vec3>>(source, target);
} static PyObject* _W_5 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) {
 try {
 PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "copyGridToArrayVec3" , !noTiming ); PyObject *_retval = 0; {
 ArgLocker _lock; const Grid<Vec3>& source = *_args.getPtr<Grid<Vec3> >("source",0,&_lock); PyArrayContainer target = _args.get<PyArrayContainer >("target",1,&_lock);   _retval = getPyNone(); copyGridToArrayVec3(source,target);  _args.check(); }
 pbFinalizePlugin(parent,"copyGridToArrayVec3", !noTiming ); return _retval; }
 catch(std::exception& e) {
 pbSetError("copyGridToArrayVec3",e.what()); return 0; }
 }
 static const Pb::Register _RP_copyGridToArrayVec3 ("","copyGridToArrayVec3",_W_5); 
extern "C" {
 void PbRegister_copyGridToArrayVec3() {
 KEEP_UNUSED(_RP_copyGridToArrayVec3); }
 }



void copyArrayToGridMAC(const PyArrayContainer source, MACGrid& target) {
	copyArrayToGridVector<MACGrid>(source, target);
} static PyObject* _W_6 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) {
 try {
 PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "copyArrayToGridMAC" , !noTiming ); PyObject *_retval = 0; {
 ArgLocker _lock; const PyArrayContainer source = _args.get<PyArrayContainer >("source",0,&_lock); MACGrid& target = *_args.getPtr<MACGrid >("target",1,&_lock);   _retval = getPyNone(); copyArrayToGridMAC(source,target);  _args.check(); }
 pbFinalizePlugin(parent,"copyArrayToGridMAC", !noTiming ); return _retval; }
 catch(std::exception& e) {
 pbSetError("copyArrayToGridMAC",e.what()); return 0; }
 }
 static const Pb::Register _RP_copyArrayToGridMAC ("","copyArrayToGridMAC",_W_6); 
extern "C" {
 void PbRegister_copyArrayToGridMAC() {
 KEEP_UNUSED(_RP_copyArrayToGridMAC); }
 }



void copyGridToArrayMAC(const MACGrid& source, PyArrayContainer target) {
	copyGridToArrayVector<MACGrid>(source, target);
} static PyObject* _W_7 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) {
 try {
 PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "copyGridToArrayMAC" , !noTiming ); PyObject *_retval = 0; {
 ArgLocker _lock; const MACGrid& source = *_args.getPtr<MACGrid >("source",0,&_lock); PyArrayContainer target = _args.get<PyArrayContainer >("target",1,&_lock);   _retval = getPyNone(); copyGridToArrayMAC(source,target);  _args.check(); }
 pbFinalizePlugin(parent,"copyGridToArrayMAC", !noTiming ); return _retval; }
 catch(std::exception& e) {
 pbSetError("copyGridToArrayMAC",e.what()); return 0; }
 }
 static const Pb::Register _RP_copyGridToArrayMAC ("","copyGridToArrayMAC",_W_7); 
extern "C" {
 void PbRegister_copyGridToArrayMAC() {
 KEEP_UNUSED(_RP_copyGridToArrayMAC); }
 }



//====================================================================================================
// pdata conversion functions
//----------------------------------------------------------------------------------------------------

template<typename T>
void numpyToParticleDataImpl(const PyArrayContainer source, ParticleDataImpl<T> &target) {
	assertMsg(source.TotalSize == target.size(), "The size of the numpy array doesn't match the size of the pdata field!");
	std::copy(reinterpret_cast<const T*>(source.pData), reinterpret_cast<const T*>(source.pData)+source.TotalSize,  &(target[0]));
}
template<typename T>
void particleDataImplToNumpy(const ParticleDataImpl<T> &source, PyArrayContainer target) {
	assertMsg(target.TotalSize == source.size(), "The size of the numpy array doesn't match the size of the pdata field!");
	std::copy(&(source[0]), &(source[0])+target.TotalSize, reinterpret_cast<T*>(target.pData));
}

// python interface

void copyArrayToPdataInt(const PyArrayContainer source, ParticleDataImpl<int> &target) { 
	numpyToParticleDataImpl<int>(source, target); 
} static PyObject* _W_8 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) {
 try {
 PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "copyArrayToPdataInt" , !noTiming ); PyObject *_retval = 0; {
 ArgLocker _lock; const PyArrayContainer source = _args.get<PyArrayContainer >("source",0,&_lock); ParticleDataImpl<int> & target = *_args.getPtr<ParticleDataImpl<int>  >("target",1,&_lock);   _retval = getPyNone(); copyArrayToPdataInt(source,target);  _args.check(); }
 pbFinalizePlugin(parent,"copyArrayToPdataInt", !noTiming ); return _retval; }
 catch(std::exception& e) {
 pbSetError("copyArrayToPdataInt",e.what()); return 0; }
 }
 static const Pb::Register _RP_copyArrayToPdataInt ("","copyArrayToPdataInt",_W_8); 
extern "C" {
 void PbRegister_copyArrayToPdataInt() {
 KEEP_UNUSED(_RP_copyArrayToPdataInt); }
 }


void copyPdataToArrayInt(const ParticleDataImpl<int> &source, PyArrayContainer target) { 
	particleDataImplToNumpy<int>(source, target); 
} static PyObject* _W_9 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) {
 try {
 PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "copyPdataToArrayInt" , !noTiming ); PyObject *_retval = 0; {
 ArgLocker _lock; const ParticleDataImpl<int> & source = *_args.getPtr<ParticleDataImpl<int>  >("source",0,&_lock); PyArrayContainer target = _args.get<PyArrayContainer >("target",1,&_lock);   _retval = getPyNone(); copyPdataToArrayInt(source,target);  _args.check(); }
 pbFinalizePlugin(parent,"copyPdataToArrayInt", !noTiming ); return _retval; }
 catch(std::exception& e) {
 pbSetError("copyPdataToArrayInt",e.what()); return 0; }
 }
 static const Pb::Register _RP_copyPdataToArrayInt ("","copyPdataToArrayInt",_W_9); 
extern "C" {
 void PbRegister_copyPdataToArrayInt() {
 KEEP_UNUSED(_RP_copyPdataToArrayInt); }
 }



void copyArrayToPdataReal(const PyArrayContainer source, ParticleDataImpl<Real> &target) { 
	numpyToParticleDataImpl<Real>(source, target); 
} static PyObject* _W_10 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) {
 try {
 PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "copyArrayToPdataReal" , !noTiming ); PyObject *_retval = 0; {
 ArgLocker _lock; const PyArrayContainer source = _args.get<PyArrayContainer >("source",0,&_lock); ParticleDataImpl<Real> & target = *_args.getPtr<ParticleDataImpl<Real>  >("target",1,&_lock);   _retval = getPyNone(); copyArrayToPdataReal(source,target);  _args.check(); }
 pbFinalizePlugin(parent,"copyArrayToPdataReal", !noTiming ); return _retval; }
 catch(std::exception& e) {
 pbSetError("copyArrayToPdataReal",e.what()); return 0; }
 }
 static const Pb::Register _RP_copyArrayToPdataReal ("","copyArrayToPdataReal",_W_10); 
extern "C" {
 void PbRegister_copyArrayToPdataReal() {
 KEEP_UNUSED(_RP_copyArrayToPdataReal); }
 }


void copyPdataToArrayReal(const ParticleDataImpl<Real> &source, PyArrayContainer target) { 
	particleDataImplToNumpy<Real>(source, target); 
} static PyObject* _W_11 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) {
 try {
 PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "copyPdataToArrayReal" , !noTiming ); PyObject *_retval = 0; {
 ArgLocker _lock; const ParticleDataImpl<Real> & source = *_args.getPtr<ParticleDataImpl<Real>  >("source",0,&_lock); PyArrayContainer target = _args.get<PyArrayContainer >("target",1,&_lock);   _retval = getPyNone(); copyPdataToArrayReal(source,target);  _args.check(); }
 pbFinalizePlugin(parent,"copyPdataToArrayReal", !noTiming ); return _retval; }
 catch(std::exception& e) {
 pbSetError("copyPdataToArrayReal",e.what()); return 0; }
 }
 static const Pb::Register _RP_copyPdataToArrayReal ("","copyPdataToArrayReal",_W_11); 
extern "C" {
 void PbRegister_copyPdataToArrayReal() {
 KEEP_UNUSED(_RP_copyPdataToArrayReal); }
 }



void copyArrayToPdataVec3(const PyArrayContainer source, ParticleDataImpl<Vec3> &target) {
	assertMsg(source.TotalSize == target.size()*3, "The size of the numpy array doesn't match the size of the pdata field!");
	std::copy(reinterpret_cast<const Real*>(source.pData), reinterpret_cast<const Real*>(source.pData)+source.TotalSize,  &(target[0][0]));
} static PyObject* _W_12 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) {
 try {
 PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "copyArrayToPdataVec3" , !noTiming ); PyObject *_retval = 0; {
 ArgLocker _lock; const PyArrayContainer source = _args.get<PyArrayContainer >("source",0,&_lock); ParticleDataImpl<Vec3> & target = *_args.getPtr<ParticleDataImpl<Vec3>  >("target",1,&_lock);   _retval = getPyNone(); copyArrayToPdataVec3(source,target);  _args.check(); }
 pbFinalizePlugin(parent,"copyArrayToPdataVec3", !noTiming ); return _retval; }
 catch(std::exception& e) {
 pbSetError("copyArrayToPdataVec3",e.what()); return 0; }
 }
 static const Pb::Register _RP_copyArrayToPdataVec3 ("","copyArrayToPdataVec3",_W_12); 
extern "C" {
 void PbRegister_copyArrayToPdataVec3() {
 KEEP_UNUSED(_RP_copyArrayToPdataVec3); }
 }


void copyPdataToArrayVec3(const ParticleDataImpl<Vec3> &source, PyArrayContainer target) {
	assertMsg(target.TotalSize == source.size()*3, "The size of the numpy array doesn't match the size of the pdata field!");
	std::copy(&(source[0][0]), &(source[0][0])+target.TotalSize, reinterpret_cast<Real*>(target.pData));
} static PyObject* _W_13 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) {
 try {
 PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "copyPdataToArrayVec3" , !noTiming ); PyObject *_retval = 0; {
 ArgLocker _lock; const ParticleDataImpl<Vec3> & source = *_args.getPtr<ParticleDataImpl<Vec3>  >("source",0,&_lock); PyArrayContainer target = _args.get<PyArrayContainer >("target",1,&_lock);   _retval = getPyNone(); copyPdataToArrayVec3(source,target);  _args.check(); }
 pbFinalizePlugin(parent,"copyPdataToArrayVec3", !noTiming ); return _retval; }
 catch(std::exception& e) {
 pbSetError("copyPdataToArrayVec3",e.what()); return 0; }
 }
 static const Pb::Register _RP_copyPdataToArrayVec3 ("","copyPdataToArrayVec3",_W_13); 
extern "C" {
 void PbRegister_copyPdataToArrayVec3() {
 KEEP_UNUSED(_RP_copyPdataToArrayVec3); }
 }



} // manta



