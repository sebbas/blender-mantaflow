




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/Users/user/Developer/Xcode Projects/blenderFireIntegration/mantaflowgit/source/grid.cpp"
/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Grid representation
 *
 ******************************************************************************/

#include "grid.h"
#include "levelset.h"
#include "kernel.h"
#include <limits>
#include <sstream>
#include <cstring>
#include "fileio.h"

using namespace std;
namespace Manta {

//******************************************************************************
// GridBase members

GridBase::GridBase (FluidSolver* parent) 
	: PbClass(parent), mType(TypeNone)
{
	checkParent();
	m3D = getParent()->is3D();
}

//******************************************************************************
// Grid<T> members

// helpers to set type
template<class T> inline GridBase::GridType typeList() { return GridBase::TypeNone; }
template<> inline GridBase::GridType typeList<Real>()  { return GridBase::TypeReal; }
template<> inline GridBase::GridType typeList<int>()   { return GridBase::TypeInt;  }
template<> inline GridBase::GridType typeList<Vec3>()  { return GridBase::TypeVec3; }

template<class T>
Grid<T>::Grid(FluidSolver* parent, bool show)
	: GridBase(parent)
{
	mType = typeList<T>();
	mSize = parent->getGridSize();
	mData = parent->getGridPointer<T>();
	
	mStrideZ = parent->is2D() ? 0 : (mSize.x * mSize.y);
	mDx = 1.0 / mSize.max();
	clear();
	setHidden(!show);
}

template<class T>
Grid<T>::Grid(const Grid<T>& a) : GridBase(a.getParent()) {
	mSize = a.mSize;
	mType = a.mType;
	mStrideZ = a.mStrideZ;
	mDx = a.mDx;
	FluidSolver *gp = a.getParent();
	mData = gp->getGridPointer<T>();
	memcpy(mData, a.mData, sizeof(T) * a.mSize.x * a.mSize.y * a.mSize.z);
}

template<class T>
Grid<T>::~Grid() {
	mParent->freeGridPointer<T>(mData);
}

template<class T>
void Grid<T>::clear() {
	memset(mData, 0, sizeof(T) * mSize.x * mSize.y * mSize.z);    
}

template<class T>
void Grid<T>::swap(Grid<T>& other) {
	if (other.getSizeX() != getSizeX() || other.getSizeY() != getSizeY() || other.getSizeZ() != getSizeZ())
		errMsg("Grid::swap(): Grid dimensions mismatch.");
	
	T* dswap = other.mData;
	other.mData = mData;
	mData = dswap;
}

template<class T>
void Grid<T>::load(string name) {
	if (name.find_last_of('.') == string::npos)
		errMsg("file '" + name + "' does not have an extension");
	string ext = name.substr(name.find_last_of('.'));
	if (ext == ".raw")
		readGridRaw(name, this);
	else if (ext == ".uni")
		readGridUni(name, this);
	else
		errMsg("file '" + name +"' filetype not supported");
}

template<class T>
void Grid<T>::save(string name) {
	if (name.find_last_of('.') == string::npos)
		errMsg("file '" + name + "' does not have an extension");
	string ext = name.substr(name.find_last_of('.'));
	if (ext == ".raw")
		writeGridRaw(name, this);
	else if (ext == ".uni")
		writeGridUni(name, this);
	else if (ext == ".vol")
		writeGridVol(name, this);
	else if (ext == ".txt")
		writeGridTxt(name, this);
	else
		errMsg("file '" + name +"' filetype not supported");
}

//******************************************************************************
// Grid<T> operators

//! Kernel: Compute min value of Real grid

 struct CompMinReal : public KernelBase { CompMinReal(Grid<Real>& val) :  KernelBase(&val,0) ,val(val) ,minVal(std::numeric_limits<Real>::max())  { run(); }  inline void op(int idx, Grid<Real>& val ,Real& minVal)  {
	if (val[idx] < minVal)
		minVal = val[idx];
}   inline operator Real () { return minVal; } inline Real  & getRet() { return minVal; }  inline Grid<Real>& getArg0() { return val; } typedef Grid<Real> type0; void run() {  const int _sz = size; for (int i=0; i < _sz; i++) op(i, val,minVal);  } Grid<Real>& val;  Real minVal;  };

//! Kernel: Compute max value of Real grid

 struct CompMaxReal : public KernelBase { CompMaxReal(Grid<Real>& val) :  KernelBase(&val,0) ,val(val) ,maxVal(-std::numeric_limits<Real>::max())  { run(); }  inline void op(int idx, Grid<Real>& val ,Real& maxVal)  {
	if (val[idx] > maxVal)
		maxVal = val[idx];
}   inline operator Real () { return maxVal; } inline Real  & getRet() { return maxVal; }  inline Grid<Real>& getArg0() { return val; } typedef Grid<Real> type0; void run() {  const int _sz = size; for (int i=0; i < _sz; i++) op(i, val,maxVal);  } Grid<Real>& val;  Real maxVal;  };

//! Kernel: Compute min value of int grid

 struct CompMinInt : public KernelBase { CompMinInt(Grid<int>& val) :  KernelBase(&val,0) ,val(val) ,minVal(std::numeric_limits<int>::max())  { run(); }  inline void op(int idx, Grid<int>& val ,int& minVal)  {
	if (val[idx] < minVal)
		minVal = val[idx];
}   inline operator int () { return minVal; } inline int  & getRet() { return minVal; }  inline Grid<int>& getArg0() { return val; } typedef Grid<int> type0; void run() {  const int _sz = size; for (int i=0; i < _sz; i++) op(i, val,minVal);  } Grid<int>& val;  int minVal;  };

//! Kernel: Compute max value of int grid

 struct CompMaxInt : public KernelBase { CompMaxInt(Grid<int>& val) :  KernelBase(&val,0) ,val(val) ,maxVal(-std::numeric_limits<int>::max())  { run(); }  inline void op(int idx, Grid<int>& val ,int& maxVal)  {
	if (val[idx] > maxVal)
		maxVal = val[idx];
}   inline operator int () { return maxVal; } inline int  & getRet() { return maxVal; }  inline Grid<int>& getArg0() { return val; } typedef Grid<int> type0; void run() {  const int _sz = size; for (int i=0; i < _sz; i++) op(i, val,maxVal);  } Grid<int>& val;  int maxVal;  };

//! Kernel: Compute min norm of vec grid

 struct CompMinVec : public KernelBase { CompMinVec(Grid<Vec3>& val) :  KernelBase(&val,0) ,val(val) ,minVal(std::numeric_limits<Real>::max())  { run(); }  inline void op(int idx, Grid<Vec3>& val ,Real& minVal)  {
	const Real s = normSquare(val[idx]);
	if (s < minVal)
		minVal = s;
}   inline operator Real () { return minVal; } inline Real  & getRet() { return minVal; }  inline Grid<Vec3>& getArg0() { return val; } typedef Grid<Vec3> type0; void run() {  const int _sz = size; for (int i=0; i < _sz; i++) op(i, val,minVal);  } Grid<Vec3>& val;  Real minVal;  };

//! Kernel: Compute max norm of vec grid

 struct CompMaxVec : public KernelBase { CompMaxVec(Grid<Vec3>& val) :  KernelBase(&val,0) ,val(val) ,maxVal(-std::numeric_limits<Real>::max())  { run(); }  inline void op(int idx, Grid<Vec3>& val ,Real& maxVal)  {
	const Real s = normSquare(val[idx]);
	if (s > maxVal)
		maxVal = s;
}   inline operator Real () { return maxVal; } inline Real  & getRet() { return maxVal; }  inline Grid<Vec3>& getArg0() { return val; } typedef Grid<Vec3> type0; void run() {  const int _sz = size; for (int i=0; i < _sz; i++) op(i, val,maxVal);  } Grid<Vec3>& val;  Real maxVal;  };


template<class T> Grid<T>& Grid<T>::safeDivide (const Grid<T>& a) {
	gridSafeDiv<T> (*this, a);
	return *this;
}
template<class T> Grid<T>& Grid<T>::copyFrom (const Grid<T>& a) {
	assertMsg (a.mSize.x == mSize.x && a.mSize.y == mSize.y && a.mSize.z == mSize.z, "different grid resolutions "<<a.mSize<<" vs "<<this->mSize );
	memcpy(mData, a.mData, sizeof(T) * mSize.x * mSize.y * mSize.z);
	mType = a.mType; // copy type marker
	return *this;
}
/*template<class T> Grid<T>& Grid<T>::operator= (const Grid<T>& a) {
	note: do not use , use copyFrom instead
}*/

template <class T>  struct knGridSetConstReal : public KernelBase { knGridSetConstReal(Grid<T>& me, T val) :  KernelBase(&me,0) ,me(me),val(val)   { run(); }  inline void op(int idx, Grid<T>& me, T val )  { me[idx]  = val; }   inline Grid<T>& getArg0() { return me; } typedef Grid<T> type0;inline T& getArg1() { return val; } typedef T type1; void run() {  const int _sz = size; for (int i=0; i < _sz; i++) op(i, me,val);  } Grid<T>& me; T val;   };
template <class T>  struct knGridAddConstReal : public KernelBase { knGridAddConstReal(Grid<T>& me, T val) :  KernelBase(&me,0) ,me(me),val(val)   { run(); }  inline void op(int idx, Grid<T>& me, T val )  { me[idx] += val; }   inline Grid<T>& getArg0() { return me; } typedef Grid<T> type0;inline T& getArg1() { return val; } typedef T type1; void run() {  const int _sz = size; for (int i=0; i < _sz; i++) op(i, me,val);  } Grid<T>& me; T val;   };
template <class T>  struct knGridMultConst : public KernelBase { knGridMultConst(Grid<T>& me, T val) :  KernelBase(&me,0) ,me(me),val(val)   { run(); }  inline void op(int idx, Grid<T>& me, T val )  { me[idx] *= val; }   inline Grid<T>& getArg0() { return me; } typedef Grid<T> type0;inline T& getArg1() { return val; } typedef T type1; void run() {  const int _sz = size; for (int i=0; i < _sz; i++) op(i, me,val);  } Grid<T>& me; T val;   };
template <class T>  struct knGridClamp : public KernelBase { knGridClamp(Grid<T>& me, T min, T max) :  KernelBase(&me,0) ,me(me),min(min),max(max)   { run(); }  inline void op(int idx, Grid<T>& me, T min, T max )  { me[idx] = clamp( me[idx], min, max); }   inline Grid<T>& getArg0() { return me; } typedef Grid<T> type0;inline T& getArg1() { return min; } typedef T type1;inline T& getArg2() { return max; } typedef T type2; void run() {  const int _sz = size; for (int i=0; i < _sz; i++) op(i, me,min,max);  } Grid<T>& me; T min; T max;   };

template<class T> void Grid<T>::add(const Grid<T>& a) {
	gridAdd<T,T>(*this, a);
}
template<class T> void Grid<T>::sub(const Grid<T>& a) {
	gridSub<T,T>(*this, a);
}
template<class T> void Grid<T>::addScaled(const Grid<T>& a, const T& factor) { 
	gridScaledAdd<T,T> (*this, a, factor); 
}
template<class T> void Grid<T>::setConst(T a) {
	knGridSetConstReal<T>( *this, T(a) );
}
template<class T> void Grid<T>::addConst(T a) {
	knGridAddConstReal<T>( *this, T(a) );
}
template<class T> void Grid<T>::multConst(T a) {
	knGridMultConst<T>( *this, a );
}

template<class T> void Grid<T>::mult(const Grid<T>& a) {
	gridMult<T,T> (*this, a);
}

template<class T> void Grid<T>::clamp(Real min, Real max) {
	knGridClamp<T> (*this, T(min), T(max) );
}

template<> Real Grid<Real>::getMax() {
	return CompMaxReal (*this);
}
template<> Real Grid<Real>::getMin() {
	return CompMinReal (*this);
}
template<> Real Grid<Real>::getMaxAbs() {
	Real amin = CompMinReal (*this);
	Real amax = CompMaxReal (*this);
	return max( fabs(amin), fabs(amax));
}
template<> Real Grid<Vec3>::getMax() {
	return sqrt(CompMaxVec (*this));
}
template<> Real Grid<Vec3>::getMin() { 
	return sqrt(CompMinVec (*this));
}
template<> Real Grid<Vec3>::getMaxAbs() {
	return sqrt(CompMaxVec (*this));
}
template<> Real Grid<int>::getMax() {
	return (Real) CompMaxInt (*this);
}
template<> Real Grid<int>::getMin() {
	return (Real) CompMinInt (*this);
}
template<> Real Grid<int>::getMaxAbs() {
	int amin = CompMinInt (*this);
	int amax = CompMaxInt (*this);
	return max( fabs((Real)amin), fabs((Real)amax));
}
template<class T> void Grid<T>::writeGridToMemory(const std::string& memLoc, const std::string& sizeAllowed) {
	if (memLoc == "" ||memLoc == "0" ){
		debMsg("Cant write grid to NULL pointer",1);
		return;
	}
	istringstream iss(sizeAllowed);
    size_t sizeAllowed_num;
    iss >> sizeAllowed_num;
	if (sizeof(T) * mSize.x * mSize.y * mSize.z != sizeAllowed_num){
		debMsg("Cant write grid with incompatible size",1);
		return;
	}
	stringstream ss(memLoc);
	void *gridPointer = NULL;
	ss >> gridPointer;
	memcpy(gridPointer, mData, sizeAllowed_num);
}
template<class T> void Grid<T>::readGridFromMemory(const std::string& memLoc, int x, int y, int z) {
	debMsg("Reading grid from " + memLoc,1);
	if (memLoc == "" ||memLoc == "0" ){
		debMsg("Can not write grid to NULL pointer",1);
		return;
	}
	if (x != mSize.x || y != mSize.y || z != mSize.z)
	{
		debMsg("Can not write grid with different domain size",1);
		return;
	}
	stringstream ss(memLoc);
	void *gridPointer = NULL;
	ss >> gridPointer;
	memcpy(mData, gridPointer, sizeof(T) * x * y * z);
}
template<class T> void Grid<T>::readAdaptiveGridFromMemory(const std::string& memLoc, Vec3i minSize, Vec3i maxSize) {
	if (memLoc == "" ||memLoc == "0" ){
		debMsg("Can not write grid to NULL pointer",1);
		return;
	}
	if (minSize.x < 0 || minSize.y < 0 || minSize.z < 0){
		debMsg("Adaptive grid smaller than 0",1);
		return;
	}
	if (maxSize.x > mSize.x || maxSize.y > mSize.y || maxSize.z > mSize.z){
		debMsg("Adaptive grid larger than current",1);
		return;
	}
	Vec3i adaptiveSize = maxSize - minSize;
	stringstream ss(memLoc);
	void *gridPointer = NULL;
	ss >> gridPointer;
	float *data_Array = (float* )gridPointer;
	for (int x = 0; x < adaptiveSize.x; ++x){
		for (int y = 0; y < adaptiveSize.y; ++y){
			for (int z = 0; z < adaptiveSize.z; ++z){
				get(x + minSize.x, y + minSize.y, z + minSize.z) = data_Array[x  + adaptiveSize.x * y + adaptiveSize.x * adaptiveSize.y * z];
			}
		}
	}
}
template<class T> std::string Grid<T>::getDataPointer() {
	ostringstream ss;
	ss << mData ;
	return ss.str();
}

// compute maximal diference of two cells in the grid
// used for testing system

Real gridMaxDiff(Grid<Real>& g1, Grid<Real>& g2 ) {
	double maxVal = 0.;
	FOR_IJK(g1) {
		maxVal = std::max(maxVal, (double)fabs( g1(i,j,k)-g2(i,j,k) ));
	}
	return maxVal; 
} static PyObject* _W_0 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "gridMaxDiff" ); PyObject *_retval = 0; { ArgLocker _lock; Grid<Real>& g1 = *_args.getPtr<Grid<Real> >("g1",0,&_lock); Grid<Real>& g2 = *_args.getPtr<Grid<Real> >("g2",1,&_lock);   _retval = toPy(gridMaxDiff(g1,g2));  _args.check(); } pbFinalizePlugin(parent,"gridMaxDiff" ); return _retval; } catch(std::exception& e) { pbSetError("gridMaxDiff",e.what()); return 0; } } static const Pb::Register _RP_gridMaxDiff ("","gridMaxDiff",_W_0); 

Real gridMaxDiffInt(Grid<int>& g1, Grid<int>& g2 ) {
	double maxVal = 0.;
	FOR_IJK(g1) {
		maxVal = std::max(maxVal, (double)fabs( (double)g1(i,j,k)-g2(i,j,k) ));
	}
	return maxVal; 
} static PyObject* _W_1 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "gridMaxDiffInt" ); PyObject *_retval = 0; { ArgLocker _lock; Grid<int>& g1 = *_args.getPtr<Grid<int> >("g1",0,&_lock); Grid<int>& g2 = *_args.getPtr<Grid<int> >("g2",1,&_lock);   _retval = toPy(gridMaxDiffInt(g1,g2));  _args.check(); } pbFinalizePlugin(parent,"gridMaxDiffInt" ); return _retval; } catch(std::exception& e) { pbSetError("gridMaxDiffInt",e.what()); return 0; } } static const Pb::Register _RP_gridMaxDiffInt ("","gridMaxDiffInt",_W_1); 

Real gridMaxDiffVec3(Grid<Vec3>& g1, Grid<Vec3>& g2 ) {
	double maxVal = 0.;
	FOR_IJK(g1) {
		// accumulate differences with double precision
		// note - don't use norm here! should be as precise as possible...
		double d = 0.;
		for(int c=0; c<3; ++c) { 
			d += fabs( (double)g1(i,j,k)[c] - (double)g2(i,j,k)[c] );
		}
		maxVal = std::max(maxVal, d );
		//maxVal = std::max(maxVal, (double)fabs( norm(g1(i,j,k)-g2(i,j,k)) ));
	}
	return maxVal; 
} static PyObject* _W_2 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "gridMaxDiffVec3" ); PyObject *_retval = 0; { ArgLocker _lock; Grid<Vec3>& g1 = *_args.getPtr<Grid<Vec3> >("g1",0,&_lock); Grid<Vec3>& g2 = *_args.getPtr<Grid<Vec3> >("g2",1,&_lock);   _retval = toPy(gridMaxDiffVec3(g1,g2));  _args.check(); } pbFinalizePlugin(parent,"gridMaxDiffVec3" ); return _retval; } catch(std::exception& e) { pbSetError("gridMaxDiffVec3",e.what()); return 0; } } static const Pb::Register _RP_gridMaxDiffVec3 ("","gridMaxDiffVec3",_W_2); 

// simple helper functions to copy (convert) mac to vec3 , and levelset to real grids
// (are assumed to be the same for running the test cases - in general they're not!)

void copyMacToVec3(MACGrid &source, Grid<Vec3>& target) {
	FOR_IJK(target) {
		target(i,j,k) = source(i,j,k);
	}
} static PyObject* _W_3 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "copyMacToVec3" ); PyObject *_retval = 0; { ArgLocker _lock; MACGrid& source = *_args.getPtr<MACGrid >("source",0,&_lock); Grid<Vec3>& target = *_args.getPtr<Grid<Vec3> >("target",1,&_lock);   _retval = getPyNone(); copyMacToVec3(source,target);  _args.check(); } pbFinalizePlugin(parent,"copyMacToVec3" ); return _retval; } catch(std::exception& e) { pbSetError("copyMacToVec3",e.what()); return 0; } } static const Pb::Register _RP_copyMacToVec3 ("","copyMacToVec3",_W_3); 
void convertMacToVec3(MACGrid &source , Grid<Vec3> &target) { debMsg("Deprecated - do not use convertMacToVec3... use copyMacToVec3 instead",1); copyMacToVec3(source,target); } static PyObject* _W_4 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "convertMacToVec3" ); PyObject *_retval = 0; { ArgLocker _lock; MACGrid& source = *_args.getPtr<MACGrid >("source",0,&_lock); Grid<Vec3> & target = *_args.getPtr<Grid<Vec3>  >("target",1,&_lock);   _retval = getPyNone(); convertMacToVec3(source,target);  _args.check(); } pbFinalizePlugin(parent,"convertMacToVec3" ); return _retval; } catch(std::exception& e) { pbSetError("convertMacToVec3",e.what()); return 0; } } static const Pb::Register _RP_convertMacToVec3 ("","convertMacToVec3",_W_4); 


void copyLevelsetToReal(LevelsetGrid &source , Grid<Real> &target) {
	FOR_IJK(target) {
		target(i,j,k) = source(i,j,k);
	}
} static PyObject* _W_5 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "copyLevelsetToReal" ); PyObject *_retval = 0; { ArgLocker _lock; LevelsetGrid& source = *_args.getPtr<LevelsetGrid >("source",0,&_lock); Grid<Real> & target = *_args.getPtr<Grid<Real>  >("target",1,&_lock);   _retval = getPyNone(); copyLevelsetToReal(source,target);  _args.check(); } pbFinalizePlugin(parent,"copyLevelsetToReal" ); return _retval; } catch(std::exception& e) { pbSetError("copyLevelsetToReal",e.what()); return 0; } } static const Pb::Register _RP_copyLevelsetToReal ("","copyLevelsetToReal",_W_5); 
void convertLevelsetToReal(LevelsetGrid &source , Grid<Real> &target) { debMsg("Deprecated - do not use convertLevelsetToReal... use copyLevelsetToReal instead",1); copyLevelsetToReal(source,target); } static PyObject* _W_6 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "convertLevelsetToReal" ); PyObject *_retval = 0; { ArgLocker _lock; LevelsetGrid& source = *_args.getPtr<LevelsetGrid >("source",0,&_lock); Grid<Real> & target = *_args.getPtr<Grid<Real>  >("target",1,&_lock);   _retval = getPyNone(); convertLevelsetToReal(source,target);  _args.check(); } pbFinalizePlugin(parent,"convertLevelsetToReal" ); return _retval; } catch(std::exception& e) { pbSetError("convertLevelsetToReal",e.what()); return 0; } } static const Pb::Register _RP_convertLevelsetToReal ("","convertLevelsetToReal",_W_6); 


void copyMacToReal(MACGrid &source, Grid<Real> &targetX, Grid<Real> &targetY, Grid<Real> &targetZ) {
	FOR_IJK(source) {
		targetX(i,j,k) = source(i,j,k).x;
		targetY(i,j,k) = source(i,j,k).y;
		targetZ(i,j,k) = source(i,j,k).z;
	}
} static PyObject* _W_7 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "copyMacToReal" ); PyObject *_retval = 0; { ArgLocker _lock; MACGrid& source = *_args.getPtr<MACGrid >("source",0,&_lock); Grid<Real> & targetX = *_args.getPtr<Grid<Real>  >("targetX",1,&_lock); Grid<Real> & targetY = *_args.getPtr<Grid<Real>  >("targetY",2,&_lock); Grid<Real> & targetZ = *_args.getPtr<Grid<Real>  >("targetZ",3,&_lock);   _retval = getPyNone(); copyMacToReal(source,targetX,targetY,targetZ);  _args.check(); } pbFinalizePlugin(parent,"copyMacToReal" ); return _retval; } catch(std::exception& e) { pbSetError("copyMacToReal",e.what()); return 0; } } static const Pb::Register _RP_copyMacToReal ("","copyMacToReal",_W_7); 


void copyRealToMac(Grid<Real> &sourceX, Grid<Real> &sourceY, Grid<Real> &sourceZ, MACGrid &target) {
	FOR_IJK(target) {
		target(i,j,k).x = sourceX(i,j,k);
		target(i,j,k).y = sourceY(i,j,k);
		target(i,j,k).z = sourceZ(i,j,k);
	}
} static PyObject* _W_8 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "copyRealToMac" ); PyObject *_retval = 0; { ArgLocker _lock; Grid<Real> & sourceX = *_args.getPtr<Grid<Real>  >("sourceX",0,&_lock); Grid<Real> & sourceY = *_args.getPtr<Grid<Real>  >("sourceY",1,&_lock); Grid<Real> & sourceZ = *_args.getPtr<Grid<Real>  >("sourceZ",2,&_lock); MACGrid& target = *_args.getPtr<MACGrid >("target",3,&_lock);   _retval = getPyNone(); copyRealToMac(sourceX,sourceY,sourceZ,target);  _args.check(); } pbFinalizePlugin(parent,"copyRealToMac" ); return _retval; } catch(std::exception& e) { pbSetError("copyRealToMac",e.what()); return 0; } } static const Pb::Register _RP_copyRealToMac ("","copyRealToMac",_W_8); 

template<class T> void Grid<T>::printGrid(int zSlice, bool printIndex) {
	std::ostringstream out;
	out << std::endl;
	const int bnd = 1;
	FOR_IJK_BND(*this,bnd) {
		int idx = (*this).index(i,j,k);
		if(zSlice>=0 && k==zSlice) { 
			out << " ";
			if(printIndex &&  this->is3D()) out << "  "<<i<<","<<j<<","<<k <<":";
			if(printIndex && !this->is3D()) out << "  "<<i<<","<<j<<":";
			out << (*this)[idx]; 
			if(i==(*this).getSizeX()-1 -bnd) out << std::endl; 
		}
	}
	out << endl; debMsg("Printing " << this->getName() << out.str().c_str() , 1);
}

//! helper to swap components of a grid (eg for data import)
void swapComponents(Grid<Vec3>& vel, int c1=0, int c2=1, int c3=2) {
	FOR_IJK(vel) {
		Vec3 v = vel(i,j,k);
		vel(i,j,k)[0] = v[c1];
		vel(i,j,k)[1] = v[c2];
		vel(i,j,k)[2] = v[c3];
	}
} static PyObject* _W_9 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "swapComponents" ); PyObject *_retval = 0; { ArgLocker _lock; Grid<Vec3>& vel = *_args.getPtr<Grid<Vec3> >("vel",0,&_lock); int c1 = _args.getOpt<int >("c1",1,0,&_lock); int c2 = _args.getOpt<int >("c2",2,1,&_lock); int c3 = _args.getOpt<int >("c3",3,2,&_lock);   _retval = getPyNone(); swapComponents(vel,c1,c2,c3);  _args.check(); } pbFinalizePlugin(parent,"swapComponents" ); return _retval; } catch(std::exception& e) { pbSetError("swapComponents",e.what()); return 0; } } static const Pb::Register _RP_swapComponents ("","swapComponents",_W_9); 

// helper functions for UV grid data (stored grid coordinates as Vec3 values, and uv weight in entry zero)

// make uv weight accesible in python
Real getUvWeight(Grid<Vec3> &uv) { return uv[0][0]; } static PyObject* _W_10 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "getUvWeight" ); PyObject *_retval = 0; { ArgLocker _lock; Grid<Vec3> & uv = *_args.getPtr<Grid<Vec3>  >("uv",0,&_lock);   _retval = toPy(getUvWeight(uv));  _args.check(); } pbFinalizePlugin(parent,"getUvWeight" ); return _retval; } catch(std::exception& e) { pbSetError("getUvWeight",e.what()); return 0; } } static const Pb::Register _RP_getUvWeight ("","getUvWeight",_W_10); 

// note - right now the UV grids have 0 values at the border after advection... could be fixed with an extrapolation step...

// compute normalized modulo interval
static inline Real computeUvGridTime(Real t, Real resetTime) {
	return fmod( (t / resetTime), (Real)1. );
}
// create ramp function in 0..1 range with half frequency
static inline Real computeUvRamp(Real t) {
	Real uvWeight = 2. * t; 
	if (uvWeight>1.) uvWeight=2.-uvWeight;
	return uvWeight;
}

 struct knResetUvGrid : public KernelBase { knResetUvGrid(Grid<Vec3>& target) :  KernelBase(&target,0) ,target(target)   { run(); }  inline void op(int i, int j, int k, Grid<Vec3>& target )  { target(i,j,k) = Vec3((Real)i,(Real)j,(Real)k); }   inline Grid<Vec3>& getArg0() { return target; } typedef Grid<Vec3> type0; void run() {  const int _maxX = maxX; const int _maxY = maxY; for (int k=minZ; k< maxZ; k++) for (int j=0; j< _maxY; j++) for (int i=0; i< _maxX; i++) op(i,j,k, target);  } Grid<Vec3>& target;   };


void resetUvGrid(Grid<Vec3> &target) {
	knResetUvGrid reset(target); // note, llvm complains about anonymous declaration here... ?
} static PyObject* _W_11 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "resetUvGrid" ); PyObject *_retval = 0; { ArgLocker _lock; Grid<Vec3> & target = *_args.getPtr<Grid<Vec3>  >("target",0,&_lock);   _retval = getPyNone(); resetUvGrid(target);  _args.check(); } pbFinalizePlugin(parent,"resetUvGrid" ); return _retval; } catch(std::exception& e) { pbSetError("resetUvGrid",e.what()); return 0; } } static const Pb::Register _RP_resetUvGrid ("","resetUvGrid",_W_11); 

void updateUvWeight(Real resetTime, int index, int numUvs, Grid<Vec3> &uv , bool info=false) {
	const Real t   = uv.getParent()->getTime();
	Real  timeOff  = resetTime/(Real)numUvs;

	Real lastt = computeUvGridTime(t +(Real)index*timeOff - uv.getParent()->getDt(), resetTime);
	Real currt = computeUvGridTime(t +(Real)index*timeOff                  , resetTime);
	Real uvWeight = computeUvRamp(currt);

	// normalize the uvw weights , note: this is a bit wasteful...
	Real uvWTotal = 0.;
	for(int i=0; i<numUvs; ++i) {
		uvWTotal += computeUvRamp( computeUvGridTime(t +(Real)i*timeOff , resetTime) );
	}
	if(uvWTotal<=VECTOR_EPSILON) { uvWeight =  uvWTotal = 1.; }
	else                           uvWeight /= uvWTotal;

	// check for reset
	if( currt < lastt ) 
		knResetUvGrid reset( uv );

	// write new weight value to grid
	uv[0] = Vec3( uvWeight, 0.,0.);

	// print info about uv weights?
	if(info) debMsg("Uv grid "<<index<<"/"<<numUvs<< " t="<<currt<<" w="<<uvWeight<<", reset:"<<(int)(currt<lastt) , 1);
} static PyObject* _W_12 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "updateUvWeight" ); PyObject *_retval = 0; { ArgLocker _lock; Real resetTime = _args.get<Real >("resetTime",0,&_lock); int index = _args.get<int >("index",1,&_lock); int numUvs = _args.get<int >("numUvs",2,&_lock); Grid<Vec3> & uv = *_args.getPtr<Grid<Vec3>  >("uv",3,&_lock); bool info = _args.getOpt<bool >("info",4,false,&_lock);   _retval = getPyNone(); updateUvWeight(resetTime,index,numUvs,uv,info);  _args.check(); } pbFinalizePlugin(parent,"updateUvWeight" ); return _retval; } catch(std::exception& e) { pbSetError("updateUvWeight",e.what()); return 0; } } static const Pb::Register _RP_updateUvWeight ("","updateUvWeight",_W_12); 

template <class T>  struct knSetBoundary : public KernelBase { knSetBoundary(Grid<T>& grid, T value, int w) :  KernelBase(&grid,0) ,grid(grid),value(value),w(w)   { run(); }  inline void op(int i, int j, int k, Grid<T>& grid, T value, int w )  { 
	bool bnd = (i<=w || i>=grid.getSizeX()-1-w || j<=w || j>=grid.getSizeY()-1-w || (grid.is3D() && (k<=w || k>=grid.getSizeZ()-1-w)));
	if (bnd) 
		grid(i,j,k) = value;
}   inline Grid<T>& getArg0() { return grid; } typedef Grid<T> type0;inline T& getArg1() { return value; } typedef T type1;inline int& getArg2() { return w; } typedef int type2; void run() {  const int _maxX = maxX; const int _maxY = maxY; for (int k=minZ; k< maxZ; k++) for (int j=0; j< _maxY; j++) for (int i=0; i< _maxX; i++) op(i,j,k, grid,value,w);  } Grid<T>& grid; T value; int w;   };

template<class T> void Grid<T>::setBound(T value, int boundaryWidth) {
	knSetBoundary<T>( *this, value, boundaryWidth );
}


template <class T>  struct knSetBoundaryNeumann : public KernelBase { knSetBoundaryNeumann(Grid<T>& grid, int w) :  KernelBase(&grid,0) ,grid(grid),w(w)   { run(); }  inline void op(int i, int j, int k, Grid<T>& grid, int w )  { 
	bool set = false;
	int  si=i, sj=j, sk=k;
	if( i<=w) {
		si = w+1; set=true;
	}
	if( i>=grid.getSizeX()-1-w){
		si = grid.getSizeX()-1-w-1; set=true;
	}
	if( j<=w){
		sj = w+1; set=true;
	}
	if( j>=grid.getSizeY()-1-w){
		sj = grid.getSizeY()-1-w-1; set=true;
	}
	if( grid.is3D() ){
		 if( k<=w ) {
			sk = w+1; set=true;
		 }
		 if( k>=grid.getSizeZ()-1-w ) {
			sk = grid.getSizeZ()-1-w-1; set=true;
		 }
	}
	if(set)
		grid(i,j,k) = grid(si, sj, sk);
}   inline Grid<T>& getArg0() { return grid; } typedef Grid<T> type0;inline int& getArg1() { return w; } typedef int type1; void run() {  const int _maxX = maxX; const int _maxY = maxY; for (int k=minZ; k< maxZ; k++) for (int j=0; j< _maxY; j++) for (int i=0; i< _maxX; i++) op(i,j,k, grid,w);  } Grid<T>& grid; int w;   };

template<class T> void Grid<T>::setBoundNeumann(int boundaryWidth) {
	knSetBoundaryNeumann<T>( *this, boundaryWidth );
}

//! helper kernels for getGridAvg

 struct knGridTotalSum : public KernelBase { knGridTotalSum(const Grid<Real>& a, FlagGrid* flags) :  KernelBase(&a,0) ,a(a),flags(flags) ,result(0.0)  { run(); }  inline void op(int idx, const Grid<Real>& a, FlagGrid* flags ,double& result)  {
	if(flags) {	if(flags->isFluid(idx)) result += a[idx]; } 
	else      {	result += a[idx]; } 
}   inline operator double () { return result; } inline double  & getRet() { return result; }  inline const Grid<Real>& getArg0() { return a; } typedef Grid<Real> type0;inline FlagGrid* getArg1() { return flags; } typedef FlagGrid type1; void run() {  const int _sz = size; for (int i=0; i < _sz; i++) op(i, a,flags,result);  } const Grid<Real>& a; FlagGrid* flags;  double result;  };


 struct knCountFluidCells : public KernelBase { knCountFluidCells(FlagGrid& flags) :  KernelBase(&flags,0) ,flags(flags) ,numEmpty(0)  { run(); }  inline void op(int idx, FlagGrid& flags ,int& numEmpty)  { if (flags.isFluid(idx) ) numEmpty++; }   inline operator int () { return numEmpty; } inline int  & getRet() { return numEmpty; }  inline FlagGrid& getArg0() { return flags; } typedef FlagGrid type0; void run() {  const int _sz = size; for (int i=0; i < _sz; i++) op(i, flags,numEmpty);  } FlagGrid& flags;  int numEmpty;  };

//! averaged value for all cells (if flags are given, only for fluid cells)

Real getGridAvg(Grid<Real>& source, FlagGrid* flags=NULL) {
	double sum = knGridTotalSum(source, flags);

	double cells;
	if(flags) { cells = knCountFluidCells(*flags); }
	else      { cells = source.getSizeX()*source.getSizeY()*source.getSizeZ(); }

	if(cells>0.) sum *= 1./cells;
	else         sum = -1.;
	return sum;
} static PyObject* _W_13 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "getGridAvg" ); PyObject *_retval = 0; { ArgLocker _lock; Grid<Real>& source = *_args.getPtr<Grid<Real> >("source",0,&_lock); FlagGrid* flags = _args.getPtrOpt<FlagGrid >("flags",1,NULL,&_lock);   _retval = toPy(getGridAvg(source,flags));  _args.check(); } pbFinalizePlugin(parent,"getGridAvg" ); return _retval; } catch(std::exception& e) { pbSetError("getGridAvg",e.what()); return 0; } } static const Pb::Register _RP_getGridAvg ("","getGridAvg",_W_13); 

//! transfer data between real and vec3 grids

 struct knGetComponent : public KernelBase { knGetComponent(Grid<Vec3>& source, Grid<Real>& target, int component) :  KernelBase(&source,0) ,source(source),target(target),component(component)   { run(); }  inline void op(int idx, Grid<Vec3>& source, Grid<Real>& target, int component )  { 
	target[idx] = source[idx][component]; 
}   inline Grid<Vec3>& getArg0() { return source; } typedef Grid<Vec3> type0;inline Grid<Real>& getArg1() { return target; } typedef Grid<Real> type1;inline int& getArg2() { return component; } typedef int type2; void run() {  const int _sz = size; for (int i=0; i < _sz; i++) op(i, source,target,component);  } Grid<Vec3>& source; Grid<Real>& target; int component;   };
void getComponent(Grid<Vec3>& source, Grid<Real>& target, int component) { knGetComponent(source, target, component); } static PyObject* _W_14 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "getComponent" ); PyObject *_retval = 0; { ArgLocker _lock; Grid<Vec3>& source = *_args.getPtr<Grid<Vec3> >("source",0,&_lock); Grid<Real>& target = *_args.getPtr<Grid<Real> >("target",1,&_lock); int component = _args.get<int >("component",2,&_lock);   _retval = getPyNone(); getComponent(source,target,component);  _args.check(); } pbFinalizePlugin(parent,"getComponent" ); return _retval; } catch(std::exception& e) { pbSetError("getComponent",e.what()); return 0; } } static const Pb::Register _RP_getComponent ("","getComponent",_W_14); 

 struct knSetComponent : public KernelBase { knSetComponent(Grid<Real>& source, Grid<Vec3>& target, int component) :  KernelBase(&source,0) ,source(source),target(target),component(component)   { run(); }  inline void op(int idx, Grid<Real>& source, Grid<Vec3>& target, int component )  { 
	target[idx][component] = source[idx]; 
}   inline Grid<Real>& getArg0() { return source; } typedef Grid<Real> type0;inline Grid<Vec3>& getArg1() { return target; } typedef Grid<Vec3> type1;inline int& getArg2() { return component; } typedef int type2; void run() {  const int _sz = size; for (int i=0; i < _sz; i++) op(i, source,target,component);  } Grid<Real>& source; Grid<Vec3>& target; int component;   };
void setComponent(Grid<Real>& source, Grid<Vec3>& target, int component) { knSetComponent(source, target, component); } static PyObject* _W_15 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "setComponent" ); PyObject *_retval = 0; { ArgLocker _lock; Grid<Real>& source = *_args.getPtr<Grid<Real> >("source",0,&_lock); Grid<Vec3>& target = *_args.getPtr<Grid<Vec3> >("target",1,&_lock); int component = _args.get<int >("component",2,&_lock);   _retval = getPyNone(); setComponent(source,target,component);  _args.check(); } pbFinalizePlugin(parent,"setComponent" ); return _retval; } catch(std::exception& e) { pbSetError("setComponent",e.what()); return 0; } } static const Pb::Register _RP_setComponent ("","setComponent",_W_15); 

//******************************************************************************
// Specialization classes

void FlagGrid::InitMinXWall(const int &boundaryWidth, Grid<Real>& phiWalls) {
	const int w = boundaryWidth;
	FOR_IJK(phiWalls) {
		phiWalls(i,j,k) = std::min(i - w - .5, (double)phiWalls(i,j,k));
	}
}

void FlagGrid::InitMaxXWall(const int &boundaryWidth, Grid<Real>& phiWalls) {
	const int w = boundaryWidth;
	FOR_IJK(phiWalls) {
		phiWalls(i,j,k) = std::min(mSize.x-i-1.5-w, (double)phiWalls(i,j,k));
	}
}

void FlagGrid::InitMinYWall(const int &boundaryWidth, Grid<Real>& phiWalls) {
	const int w = boundaryWidth;
	FOR_IJK(phiWalls) {
		phiWalls(i,j,k) = std::min(j - w - .5, (double)phiWalls(i,j,k));
	}
}

void FlagGrid::InitMaxYWall(const int &boundaryWidth, Grid<Real>& phiWalls) {
	const int w = boundaryWidth;
	FOR_IJK(phiWalls) {
		phiWalls(i,j,k) = std::min(mSize.y-j-1.5-w, (double)phiWalls(i,j,k));
	}
}

void FlagGrid::InitMinZWall(const int &boundaryWidth, Grid<Real>& phiWalls) {
	const int w = boundaryWidth;
	FOR_IJK(phiWalls) {
		phiWalls(i,j,k) = std::min(k - w - .5, (double)phiWalls(i,j,k));
	}
}

void FlagGrid::InitMaxZWall(const int &boundaryWidth, Grid<Real>& phiWalls) {
	const int w = boundaryWidth;
	FOR_IJK(phiWalls) {
		phiWalls(i,j,k) = std::min(mSize.z-k-1.5-w, (double)phiWalls(i,j,k));
	}
}

void FlagGrid::initDomain( const int &boundaryWidth
	                     , const string &wall    
						 , const string &open    
						 , const string &inflow  
						 , const string &outflow
						 , Grid<Real>* phiWalls ) {
	
	int  types[6] = {0};
	bool set  [6] = {false};

	if(phiWalls) phiWalls->setConst(1000000000);

	for (char i = 0; i<6; ++i) {
		//min x-direction
		if(!set[0]) {
			if(open[i]=='x')         {types[0] = TypeOpen;set[0] = true;}
			else if(inflow[i]=='x')  {types[0] = TypeInflow;set[0] = true;}
			else if(outflow[i]=='x') {types[0] = TypeOutflow;set[0] = true;}
			else if(wall[i]=='x') {
				types[0]    = TypeObstacle;
				if(phiWalls) InitMinXWall(boundaryWidth, *phiWalls);
				set[0] = true;
			}			
		}
		//max x-direction
		if(!set[1]) {
			if(open[i]=='X')         {types[1] = TypeOpen;set[1] = true;}
			else if(inflow[i]=='X')  {types[1] = TypeInflow;set[1] = true;}
			else if(outflow[i]=='X') {types[1] = TypeOutflow;set[1] = true;}
			else if(wall[i]=='X')  {
				types[1]    = TypeObstacle;
				if(phiWalls) InitMaxXWall(boundaryWidth, *phiWalls);
				set[1] = true;
			}			
		}
		//min y-direction
		if(!set[2]) {
			if(open[i]=='y')         {types[2] = TypeOpen;set[2] = true;}
			else if(inflow[i]=='y')  {types[2] = TypeInflow;set[2] = true;}
			else if(outflow[i]=='y') {types[2] = TypeOutflow;set[2] = true;}
			else if(wall[i]=='y') {
				types[2]    = TypeObstacle;
				if(phiWalls) InitMinYWall(boundaryWidth, *phiWalls);
				set[2] = true;
			}			
		}
		//max y-direction
		if(!set[3]) {
			if(open[i]=='Y')         {types[3] = TypeOpen;set[3] = true;}
			else if(inflow[i]=='Y')  {types[3] = TypeInflow;set[3] = true;}
			else if(outflow[i]=='Y') {types[3] = TypeOutflow;set[3] = true;}
			else if(wall[i]=='Y') {
				types[3]    = TypeObstacle;
				if(phiWalls) InitMaxYWall(boundaryWidth, *phiWalls);
				set[3] = true;
			}			
		}
		if(this->is3D()) {
		//min z-direction
			if(!set[4]) {
				if(open[i]=='z')         {types[4] = TypeOpen;set[4] = true;}
				else if(inflow[i]=='z')  {types[4] = TypeInflow;set[4] = true;}
				else if(outflow[i]=='z') {types[4] = TypeOutflow;set[4] = true;}
				else if(wall[i]=='z') {
					types[4]    = TypeObstacle;
					if(phiWalls) InitMinZWall(boundaryWidth, *phiWalls);
					set[4] = true;
				}				
			}
			//max z-direction
			if(!set[5]) {
				if(open[i]=='Z')         {types[5] = TypeOpen;set[5] = true;}
				else if(inflow[i]=='Z')  {types[5] = TypeInflow;set[5] = true;}
				else if(outflow[i]=='Z') {types[5] = TypeOutflow;set[5] = true;}
				else if(wall[i]=='Z') {
					types[5]    = TypeObstacle;
					if(phiWalls) InitMaxZWall(boundaryWidth, *phiWalls);
					set[5] = true;
				}				
			}
		}
	}

	FOR_IDX(*this)
		mData[idx] = TypeEmpty;
		initBoundaries(boundaryWidth, types);
	
}

void FlagGrid::initBoundaries(const int &boundaryWidth, const int *types) {
	const int w = boundaryWidth;
	FOR_IJK(*this) {
		bool bnd = (i <= w);
		if (bnd) mData[index(i,j,k)] = types[0];
		bnd = (i >= mSize.x-1-w);
		if (bnd) mData[index(i,j,k)] = types[1];
		bnd = (j <= w);
		if (bnd) mData[index(i,j,k)] = types[2];
		bnd = (j >= mSize.y-1-w);
		if (bnd) mData[index(i,j,k)] = types[3];
		if(is3D()) {
			bnd = (k <= w);
			if (bnd) mData[index(i,j,k)] = types[4];
			bnd = (k >= mSize.z-1-w);
			if (bnd) mData[index(i,j,k)] = types[5];
		}
	}
}

void FlagGrid::updateFromLevelset(LevelsetGrid& levelset) {
	FOR_IDX(*this) {
		if (!isObstacle(idx)) {
			const Real phi = levelset[idx];
			if (phi <= levelset.invalidTimeValue()) continue;
			
			mData[idx] &= ~(TypeEmpty | TypeFluid); // clear empty/fluid flags
			mData[idx] |= (phi <= 0) ? TypeFluid : TypeEmpty; // set resepctive flag
		}
	}
}   

void FlagGrid::fillGrid(int type) {
	FOR_IDX(*this) {
		if ((mData[idx] & TypeObstacle)==0 && (mData[idx] & TypeInflow)==0&& (mData[idx] & TypeOutflow)==0&& (mData[idx] & TypeOpen)==0)
			mData[idx] = (mData[idx] & ~(TypeEmpty | TypeFluid)) | type;
	}
}

// explicit instantiation
template class Grid<int>;
template class Grid<Real>;
template class Grid<Vec3>;


//******************************************************************************
// enable compilation of a more complicated test data type
// enable in grid.h

#if ENABLE_GRID_TEST_DATATYPE==1
// todo fix, missing:  template<> const char* Namify<nbVector>::S = "TestDatatype";

template<> Real Grid<nbVector>::getMin() { return 0.; }
template<> Real Grid<nbVector>::getMax() { return 0.; }
template<> Real Grid<nbVector>::getMaxAbs()      { return 0.; }

 struct knNbvecTestKernel : public KernelBase { knNbvecTestKernel(Grid<nbVector>& target) :  KernelBase(&target,0) ,target(target)   { run(); }  inline void op(int i, int j, int k, Grid<nbVector>& target )  { target(i,j,k).push_back(i+j+k); }   inline Grid<nbVector>& getArg0() { return target; } typedef Grid<nbVector> type0; void run() {  const int _maxX = maxX; const int _maxY = maxY; for (int k=minZ; k< maxZ; k++) for (int j=0; j< _maxY; j++) for (int i=0; i< _maxX; i++) op(i,j,k, target);  } Grid<nbVector>& target;   };

void nbvecTestOp(Grid<nbVector> &target) {
	knNbvecTestKernel nbvecTest(target); 
} static PyObject* _W_16 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "nbvecTestOp" ); PyObject *_retval = 0; { ArgLocker _lock; Grid<nbVector> & target = *_args.getPtr<Grid<nbVector>  >("target",0,&_lock);   _retval = getPyNone(); nbvecTestOp(target);  _args.check(); } pbFinalizePlugin(parent,"nbvecTestOp" ); return _retval; } catch(std::exception& e) { pbSetError("nbvecTestOp",e.what()); return 0; } } static const Pb::Register _RP_nbvecTestOp ("","nbvecTestOp",_W_16); 

// instantiate test datatype , not really required for simulations, mostly here for demonstration purposes
template class Grid<nbVector>;
#endif // ENABLE_GRID_TEST_DATATYPE


} //namespace


