




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/Users/sbarschkis/Developer/Mantaflow/blenderIntegration/mantaflowgit/source/grid.h"
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

#ifndef _GRID_H
#define _GRID_H

#include "manta.h"
#include "vectorbase.h"
#include "interpol.h"
#include "interpolHigh.h"
#include "kernel.h"

namespace Manta {
class LevelsetGrid;
	
//! Base class for all grids
class GridBase : public PbClass {public:
	enum GridType { TypeNone = 0, TypeReal = 1, TypeInt = 2, TypeVec3 = 4, TypeMAC = 8, TypeLevelset = 16, TypeFlags = 32 };
		
	GridBase(FluidSolver* parent); static int _W_0 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = Pb::objFromPy(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(0, "GridBase::GridBase" , !noTiming ); { ArgLocker _lock; FluidSolver* parent = _args.getPtr<FluidSolver >("parent",0,&_lock);  obj = new GridBase(parent); obj->registerObject(_self, &_args); _args.check(); } pbFinalizePlugin(obj->getParent(),"GridBase::GridBase" , !noTiming ); return 0; } catch(std::exception& e) { pbSetError("GridBase::GridBase",e.what()); return -1; } }
	
	//! Get the grids X dimension
	inline int getSizeX() const { return mSize.x; }
	//! Get the grids Y dimension
	inline int getSizeY() const { return mSize.y; }
	//! Get the grids Z dimension
	inline int getSizeZ() const { return mSize.z; }
	//! Get the grids dimensions
	inline Vec3i getSize() const { return mSize; }
	
	//! Get Stride in X dimension
	inline IndexInt getStrideX() const { return 1; }
	//! Get Stride in Y dimension
	inline IndexInt getStrideY() const { return mSize.x; }
	//! Get Stride in Z dimension
	inline IndexInt getStrideZ() const { return mStrideZ; }
	
	inline Real getDx() { return mDx; }
	
	//! Check if indices are within bounds, otherwise error (should only be called when debugging)
	inline void checkIndex(int i, int j, int k) const;
	//! Check if indices are within bounds, otherwise error (should only be called when debugging)
	inline void checkIndex(IndexInt idx) const;
	//! Check if index is within given boundaries
	inline bool isInBounds(const Vec3i& p, int bnd) const;
	//! Check if index is within given boundaries
	inline bool isInBounds(const Vec3i& p) const;
	//! Check if index is within given boundaries
	inline bool isInBounds(const Vec3& p, int bnd = 0) const { return isInBounds(toVec3i(p), bnd); }
	//! Check if linear index is in the range of the array
	inline bool isInBounds(IndexInt idx) const;
	
	//! Get the type of grid
	inline GridType getType() const { return mType; }
	//! Check dimensionality
	inline bool is3D() const { return m3D; }
	
	//! Get index into the data
	inline IndexInt index(int i, int j, int k) const { DEBUG_ONLY(checkIndex(i,j,k)); return (IndexInt)i + (IndexInt)mSize.x * j + (IndexInt)mStrideZ * k; }
	//! Get index into the data
	inline IndexInt index(const Vec3i& pos) const    { DEBUG_ONLY(checkIndex(pos.x,pos.y,pos.z)); return (IndexInt)pos.x + (IndexInt)mSize.x * pos.y + (IndexInt)mStrideZ * pos.z; }

	//! grid4d compatibility functions 
	inline bool is4D() const { return false; }
	inline int getSizeT() const { return 1; }
	inline int getStrideT() const { return 0; }
	inline int index(int i, int j, int k, int unused) const { return index(i,j,k); }
	inline bool isInBounds(int i,int j, int k, int t, int bnd) const { if(t!=0) return false; return isInBounds( Vec3i(i,j,k), bnd ); }

protected:
	
	GridType mType;
	Vec3i mSize;
	Real mDx;
	bool m3D; 	// precomputed Z shift: to ensure 2D compatibility, always use this instead of sx*sy !
	IndexInt mStrideZ;  public: PbArgs _args;}
#define _C_GridBase
;

//! Grid class

template<class T> class Grid : public GridBase {public:
	//! init new grid, values are set to zero
	Grid(FluidSolver* parent, bool show = true); static int _W_1 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = Pb::objFromPy(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(0, "Grid::Grid" , !noTiming ); { ArgLocker _lock; FluidSolver* parent = _args.getPtr<FluidSolver >("parent",0,&_lock); bool show = _args.getOpt<bool >("show",1,true,&_lock);  obj = new Grid(parent,show); obj->registerObject(_self, &_args); _args.check(); } pbFinalizePlugin(obj->getParent(),"Grid::Grid" , !noTiming ); return 0; } catch(std::exception& e) { pbSetError("Grid::Grid",e.what()); return -1; } }
	//! create new & copy content from another grid
	Grid(const Grid<T>& a);
	//! return memory to solver
	virtual ~Grid();
	
	typedef T BASETYPE;
	typedef GridBase BASETYPE_GRID;
	
	void save(std::string name); static PyObject* _W_2 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Grid* pbo = dynamic_cast<Grid*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "Grid::save" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock; std::string name = _args.get<std::string >("name",0,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->save(name);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Grid::save" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("Grid::save",e.what()); return 0; } }
	void load(std::string name); static PyObject* _W_3 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Grid* pbo = dynamic_cast<Grid*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "Grid::load" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock; std::string name = _args.get<std::string >("name",0,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->load(name);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Grid::load" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("Grid::load",e.what()); return 0; } }
	
	//! set all cells to zero
	void clear(); static PyObject* _W_4 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Grid* pbo = dynamic_cast<Grid*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "Grid::clear" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock;  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->clear();  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Grid::clear" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("Grid::clear",e.what()); return 0; } }
	
	//! all kinds of access functions, use grid(), grid[] or grid.get()
	//! access data
	inline T get(int i,int j, int k) const         { return mData[index(i,j,k)]; }
	//! access data
	inline T& get(int i,int j, int k)              { return mData[index(i,j,k)]; }
	//! access data
	inline T get(IndexInt idx) const               { DEBUG_ONLY(checkIndex(idx)); return mData[idx]; }
	//! access data
	inline T get(const Vec3i& pos) const           { return mData[index(pos)]; }
	//! access data
	inline T& operator()(int i, int j, int k)      { return mData[index(i, j, k)]; }
	//! access data
	inline T operator()(int i, int j, int k) const { return mData[index(i, j, k)]; }
	//! access data
	inline T& operator()(IndexInt idx)             { DEBUG_ONLY(checkIndex(idx)); return mData[idx]; }
	//! access data
	inline T operator()(IndexInt idx) const        { DEBUG_ONLY(checkIndex(idx)); return mData[idx]; }
	//! access data
	inline T& operator()(const Vec3i& pos)         { return mData[index(pos)]; }
	//! access data
	inline T operator()(const Vec3i& pos) const    { return mData[index(pos)]; }
	//! access data
	inline T& operator[](IndexInt idx)             { DEBUG_ONLY(checkIndex(idx)); return mData[idx]; }
	//! access data
	inline const T operator[](IndexInt idx) const  { DEBUG_ONLY(checkIndex(idx)); return mData[idx]; }
	
	// interpolated access
	inline T    getInterpolated(const Vec3& pos) const { return interpol<T>(mData, mSize, mStrideZ, pos); }
	inline void setInterpolated(const Vec3& pos, const T& val, Grid<Real>& sumBuffer) const { setInterpol<T>(mData, mSize, mStrideZ, pos, val, &sumBuffer[0]); }
	// higher order interpolation (1=linear, 2=cubic)
	inline T getInterpolatedHi(const Vec3& pos, int order) const { 
		switch(order) {
		case 1:  return interpol     <T>(mData, mSize, mStrideZ, pos); 
		case 2:  return interpolCubic<T>(mData, mSize, mStrideZ, pos); 
		default: assertMsg(false, "Unknown interpolation order "<<order); }
	}
	
	// assignment / copy

	//! warning - do not use "=" for grids in python, this copies the reference! not the grid content...
	//Grid<T>& operator=(const Grid<T>& a);
	//! copy content from other grid (use this one instead of operator= !)
	Grid<T>& copyFrom(const Grid<T>& a, bool copyType=true ); static PyObject* _W_5 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Grid* pbo = dynamic_cast<Grid*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "Grid::copyFrom" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock; const Grid<T>& a = *_args.getPtr<Grid<T> >("a",0,&_lock); bool copyType = _args.getOpt<bool >("copyType",1,true ,&_lock);  pbo->_args.copy(_args);  _retval = toPy(pbo->copyFrom(a,copyType));  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Grid::copyFrom" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("Grid::copyFrom",e.what()); return 0; } } // old: { *this = a; }

	// helper functions to work with grids in scene files 

	//! add/subtract other grid
	void add(const Grid<T>& a); static PyObject* _W_6 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Grid* pbo = dynamic_cast<Grid*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "Grid::add" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock; const Grid<T>& a = *_args.getPtr<Grid<T> >("a",0,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->add(a);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Grid::add" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("Grid::add",e.what()); return 0; } }
	void sub(const Grid<T>& a); static PyObject* _W_7 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Grid* pbo = dynamic_cast<Grid*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "Grid::sub" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock; const Grid<T>& a = *_args.getPtr<Grid<T> >("a",0,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->sub(a);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Grid::sub" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("Grid::sub",e.what()); return 0; } }
	//! set all cells to constant value
	void setConst(T s); static PyObject* _W_8 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Grid* pbo = dynamic_cast<Grid*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "Grid::setConst" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock; T s = _args.get<T >("s",0,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->setConst(s);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Grid::setConst" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("Grid::setConst",e.what()); return 0; } }
	//! add constant to all grid cells
	void addConst(T s); static PyObject* _W_9 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Grid* pbo = dynamic_cast<Grid*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "Grid::addConst" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock; T s = _args.get<T >("s",0,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->addConst(s);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Grid::addConst" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("Grid::addConst",e.what()); return 0; } }
	//! add scaled other grid to current one (note, only "Real" factor, "T" type not supported here!)
	void addScaled(const Grid<T>& a, const T& factor); static PyObject* _W_10 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Grid* pbo = dynamic_cast<Grid*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "Grid::addScaled" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock; const Grid<T>& a = *_args.getPtr<Grid<T> >("a",0,&_lock); const T& factor = *_args.getPtr<T >("factor",1,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->addScaled(a,factor);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Grid::addScaled" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("Grid::addScaled",e.what()); return 0; } } 
	//! multiply contents of grid
	void mult( const Grid<T>& a); static PyObject* _W_11 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Grid* pbo = dynamic_cast<Grid*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "Grid::mult" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock; const Grid<T>& a = *_args.getPtr<Grid<T> >("a",0,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->mult(a);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Grid::mult" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("Grid::mult",e.what()); return 0; } }
	//! multiply each cell by a constant scalar value
	void multConst(T s); static PyObject* _W_12 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Grid* pbo = dynamic_cast<Grid*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "Grid::multConst" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock; T s = _args.get<T >("s",0,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->multConst(s);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Grid::multConst" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("Grid::multConst",e.what()); return 0; } }
	//! clamp content to range (for vec3, clamps each component separately)
	void clamp(Real min, Real max); static PyObject* _W_13 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Grid* pbo = dynamic_cast<Grid*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "Grid::clamp" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock; Real min = _args.get<Real >("min",0,&_lock); Real max = _args.get<Real >("max",1,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->clamp(min,max);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Grid::clamp" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("Grid::clamp",e.what()); return 0; } }
	
	// common compound operators
	//! get absolute max value in grid 
	Real getMaxAbs(); static PyObject* _W_14 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Grid* pbo = dynamic_cast<Grid*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "Grid::getMaxAbs" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock;  pbo->_args.copy(_args);  _retval = toPy(pbo->getMaxAbs());  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Grid::getMaxAbs" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("Grid::getMaxAbs",e.what()); return 0; } }
	//! get max value in grid 
	Real getMax(); static PyObject* _W_15 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Grid* pbo = dynamic_cast<Grid*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "Grid::getMax" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock;  pbo->_args.copy(_args);  _retval = toPy(pbo->getMax());  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Grid::getMax" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("Grid::getMax",e.what()); return 0; } }
	//! get min value in grid 
	Real getMin(); static PyObject* _W_16 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Grid* pbo = dynamic_cast<Grid*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "Grid::getMin" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock;  pbo->_args.copy(_args);  _retval = toPy(pbo->getMin());  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Grid::getMin" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("Grid::getMin",e.what()); return 0; } }    
	//! set all boundary cells to constant value (Dirichlet)
	void setBound(T value, int boundaryWidth=1); static PyObject* _W_17 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Grid* pbo = dynamic_cast<Grid*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "Grid::setBound" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock; T value = _args.get<T >("value",0,&_lock); int boundaryWidth = _args.getOpt<int >("boundaryWidth",1,1,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->setBound(value,boundaryWidth);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Grid::setBound" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("Grid::setBound",e.what()); return 0; } }
	//! set all boundary cells to last inner value (Neumann)
	void setBoundNeumann(int boundaryWidth=1); static PyObject* _W_18 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Grid* pbo = dynamic_cast<Grid*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "Grid::setBoundNeumann" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock; int boundaryWidth = _args.getOpt<int >("boundaryWidth",0,1,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->setBoundNeumann(boundaryWidth);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Grid::setBoundNeumann" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("Grid::setBoundNeumann",e.what()); return 0; } }

	//! for compatibility, old names:
	Real getMaxAbsValue() { return getMaxAbs(); } static PyObject* _W_19 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Grid* pbo = dynamic_cast<Grid*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "Grid::getMaxAbsValue" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock;  pbo->_args.copy(_args);  _retval = toPy(pbo->getMaxAbsValue());  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Grid::getMaxAbsValue" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("Grid::getMaxAbsValue",e.what()); return 0; } }
	Real getMaxValue() { return getMax(); } static PyObject* _W_20 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Grid* pbo = dynamic_cast<Grid*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "Grid::getMaxValue" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock;  pbo->_args.copy(_args);  _retval = toPy(pbo->getMaxValue());  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Grid::getMaxValue" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("Grid::getMaxValue",e.what()); return 0; } }
	Real getMinValue() { return getMin(); } static PyObject* _W_21 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Grid* pbo = dynamic_cast<Grid*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "Grid::getMinValue" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock;  pbo->_args.copy(_args);  _retval = toPy(pbo->getMinValue());  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Grid::getMinValue" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("Grid::getMinValue",e.what()); return 0; } }

	//! get data pointer of grid
	std::string getDataPointer(); static PyObject* _W_22 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Grid* pbo = dynamic_cast<Grid*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "Grid::getDataPointer" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock;  pbo->_args.copy(_args);  _retval = toPy(pbo->getDataPointer());  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Grid::getDataPointer" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("Grid::getDataPointer",e.what()); return 0; } }
	
	//! debugging helper, print grid from python
	void printGrid(int zSlice=-1, bool printIndex=false); static PyObject* _W_23 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Grid* pbo = dynamic_cast<Grid*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "Grid::printGrid" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock; int zSlice = _args.getOpt<int >("zSlice",0,-1,&_lock); bool printIndex = _args.getOpt<bool >("printIndex",1,false,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->printGrid(zSlice,printIndex);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Grid::printGrid" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("Grid::printGrid",e.what()); return 0; } } 

	// c++ only operators
	template<class S> Grid<T>& operator+=(const Grid<S>& a);
	template<class S> Grid<T>& operator+=(const S& a);
	template<class S> Grid<T>& operator-=(const Grid<S>& a);
	template<class S> Grid<T>& operator-=(const S& a);
	template<class S> Grid<T>& operator*=(const Grid<S>& a);
	template<class S> Grid<T>& operator*=(const S& a);
	template<class S> Grid<T>& operator/=(const Grid<S>& a);
	template<class S> Grid<T>& operator/=(const S& a);
	Grid<T>& safeDivide(const Grid<T>& a);    
	
	//! Swap data with another grid (no actual data is moved)
	void swap(Grid<T>& other);

	//! grid4d compatibility functions 
	inline T& operator()(int i, int j, int k, int unused)       { return mData[index(i, j, k)]; }
	inline T operator() (int i, int j, int k, int unused) const { return mData[index(i, j, k)]; }

protected: 	T* mData; public: PbArgs _args;}
#define _C_Grid
;

// Python doesn't know about templates: explicit aliases needed




//! Special function for staggered grids
class MACGrid : public Grid<Vec3> {public:
	MACGrid(FluidSolver* parent, bool show=true) :Grid<Vec3>(parent,show){ 
		mType = (GridType)(TypeMAC | TypeVec3); } static int _W_24 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = Pb::objFromPy(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(0, "MACGrid::MACGrid" , !noTiming ); { ArgLocker _lock; FluidSolver* parent = _args.getPtr<FluidSolver >("parent",0,&_lock); bool show = _args.getOpt<bool >("show",1,true,&_lock);  obj = new MACGrid(parent,show); obj->registerObject(_self, &_args); _args.check(); } pbFinalizePlugin(obj->getParent(),"MACGrid::MACGrid" , !noTiming ); return 0; } catch(std::exception& e) { pbSetError("MACGrid::MACGrid",e.what()); return -1; } }
	
	// specialized functions for interpolating MAC information
	inline Vec3 getCentered(int i, int j, int k) const;
	inline Vec3 getCentered(const Vec3i& pos) const { return getCentered(pos.x, pos.y, pos.z); }
	inline Vec3 getAtMACX(int i, int j, int k) const;
	inline Vec3 getAtMACY(int i, int j, int k) const;
	inline Vec3 getAtMACZ(int i, int j, int k) const;
	// interpolation
	inline Vec3 getInterpolated(const Vec3& pos) const { return interpolMAC(mData, mSize, mStrideZ, pos); }
	inline void setInterpolated(const Vec3& pos, const Vec3& val, Vec3* tmp) { return setInterpolMAC(mData, mSize, mStrideZ, pos, val, tmp); }
	inline Vec3 getInterpolatedHi(const Vec3& pos, int order) const { 
		switch(order) {
		case 1:  return interpolMAC     (mData, mSize, mStrideZ, pos); 
		case 2:  return interpolCubicMAC(mData, mSize, mStrideZ, pos); 
		default: assertMsg(false, "Unknown interpolation order "<<order); }
	}
	// specials for mac grid:
	template<int comp> inline Real getInterpolatedComponent(Vec3 pos) const { return interpolComponent<comp>(mData, mSize, mStrideZ, pos); }
	template<int comp> inline Real getInterpolatedComponentHi(const Vec3& pos, int order) const { 
		switch(order) {
		case 1:  return interpolComponent<comp>(mData, mSize, mStrideZ, pos); 
		case 2:  return interpolCubicMAC(mData, mSize, mStrideZ, pos)[comp];  // warning - not yet optimized
		default: assertMsg(false, "Unknown interpolation order "<<order); }
	}
	 protected: public: PbArgs _args;}
#define _C_MACGrid
;

//! Special functions for FlagGrid
class FlagGrid : public Grid<int> {public:
	FlagGrid(FluidSolver* parent, int dim=3, bool show=true) :Grid<int>(parent,show){ 
		mType = (GridType)(TypeFlags | TypeInt); } static int _W_25 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = Pb::objFromPy(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(0, "FlagGrid::FlagGrid" , !noTiming ); { ArgLocker _lock; FluidSolver* parent = _args.getPtr<FluidSolver >("parent",0,&_lock); int dim = _args.getOpt<int >("dim",1,3,&_lock); bool show = _args.getOpt<bool >("show",2,true,&_lock);  obj = new FlagGrid(parent,dim,show); obj->registerObject(_self, &_args); _args.check(); } pbFinalizePlugin(obj->getParent(),"FlagGrid::FlagGrid" , !noTiming ); return 0; } catch(std::exception& e) { pbSetError("FlagGrid::FlagGrid",e.what()); return -1; } }
	
	//! types of cells, in/outflow can be combined, e.g., TypeFluid|TypeInflow
	enum CellType { 
		TypeNone     = 0,
		TypeFluid    = 1,
		TypeObstacle = 2,
		TypeEmpty    = 4,
		TypeInflow   = 8,
		TypeOutflow  = 16,
		TypeOpen     = 32,
		TypeStick    = 128,
		TypeReserved = 256,
		// 2^10 - 2^14 reserved for moving obstacles
		TypeZeroPressure = (1<<15) 
	};
		
	//! access for particles
	inline int getAt(const Vec3& pos) const { return mData[index((int)pos.x, (int)pos.y, (int)pos.z)]; }
			
	//! check for different flag types
	inline bool isObstacle(IndexInt idx) const { return get(idx) & TypeObstacle; }
	inline bool isObstacle(int i, int j, int k) const { return get(i,j,k) & TypeObstacle; }
	inline bool isObstacle(const Vec3i& pos) const { return get(pos) & TypeObstacle; }
	inline bool isObstacle(const Vec3& pos) const { return getAt(pos) & TypeObstacle; }
	inline bool isFluid(IndexInt idx) const { return get(idx) & TypeFluid; }
	inline bool isFluid(int i, int j, int k) const { return get(i,j,k) & TypeFluid; }
	inline bool isFluid(const Vec3i& pos) const { return get(pos) & TypeFluid; }
	inline bool isFluid(const Vec3& pos) const { return getAt(pos) & TypeFluid; }
	inline bool isInflow(IndexInt idx) const { return get(idx) & TypeInflow; }
	inline bool isInflow(int i, int j, int k) const { return get(i,j,k) & TypeInflow; }
	inline bool isInflow(const Vec3i& pos) const { return get(pos) & TypeInflow; }
	inline bool isInflow(const Vec3& pos) const { return getAt(pos) & TypeInflow; }
	inline bool isEmpty(IndexInt idx) const { return get(idx) & TypeEmpty; }
	inline bool isEmpty(int i, int j, int k) const { return get(i,j,k) & TypeEmpty; }
	inline bool isEmpty(const Vec3i& pos) const { return get(pos) & TypeEmpty; }
	inline bool isEmpty(const Vec3& pos) const { return getAt(pos) & TypeEmpty; }
	inline bool isOutflow(IndexInt idx) const { return get(idx) & TypeOutflow; }
	inline bool isOutflow(int i, int j, int k) const { return get(i, j, k) & TypeOutflow; }
	inline bool isOutflow(const Vec3i& pos) const { return get(pos) & TypeOutflow; }
	inline bool isOutflow(const Vec3& pos) const { return getAt(pos) & TypeOutflow; }
	inline bool isOpen(IndexInt idx) const { return get(idx) & TypeOpen; }
	inline bool isOpen(int i, int j, int k) const { return get(i, j, k) & TypeOpen; }
	inline bool isOpen(const Vec3i& pos) const { return get(pos) & TypeOpen; }
	inline bool isOpen(const Vec3& pos) const { return getAt(pos) & TypeOpen; }
	inline bool isStick(IndexInt idx) const { return get(idx) & TypeStick; }
	inline bool isStick(int i, int j, int k) const { return get(i,j,k) & TypeStick; }
	inline bool isStick(const Vec3i& pos) const { return get(pos) & TypeStick; }
	inline bool isStick(const Vec3& pos) const { return getAt(pos) & TypeStick; }
	
	void InitMinXWall(const int &boundaryWidth, Grid<Real>& phiWalls);
	void InitMaxXWall(const int &boundaryWidth, Grid<Real>& phiWalls);
	void InitMinYWall(const int &boundaryWidth, Grid<Real>& phiWalls);
	void InitMaxYWall(const int &boundaryWidth, Grid<Real>& phiWalls);
	void InitMinZWall(const int &boundaryWidth, Grid<Real>& phiWalls);
	void InitMaxZWall(const int &boundaryWidth, Grid<Real>& phiWalls);
	// Python callables
	




void initDomain( const int &boundaryWidth = 0 , const std::string &wall = "xXyYzZ" , const std::string &open = "      " , const std::string &inflow = "      " , const std::string &outflow = "      " , Grid<Real>* phiWalls = 0x00 ); static PyObject* _W_26 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FlagGrid* pbo = dynamic_cast<FlagGrid*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "FlagGrid::initDomain" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock; const int& boundaryWidth = _args.getOpt<int >("boundaryWidth",0,0 ,&_lock); const std::string& wall = _args.getOpt<std::string >("wall",1,"xXyYzZ" ,&_lock); const std::string& open = _args.getOpt<std::string >("open",2,"      " ,&_lock); const std::string& inflow = _args.getOpt<std::string >("inflow",3,"      " ,&_lock); const std::string& outflow = _args.getOpt<std::string >("outflow",4,"      " ,&_lock); Grid<Real>* phiWalls = _args.getPtrOpt<Grid<Real> >("phiWalls",5,0x00 ,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->initDomain(boundaryWidth,wall,open,inflow,outflow,phiWalls);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"FlagGrid::initDomain" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("FlagGrid::initDomain",e.what()); return 0; } }

	void initBoundaries( const int &boundaryWidth, const int *types );

	void updateFromLevelset(LevelsetGrid& levelset); static PyObject* _W_27 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FlagGrid* pbo = dynamic_cast<FlagGrid*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "FlagGrid::updateFromLevelset" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock; LevelsetGrid& levelset = *_args.getPtr<LevelsetGrid >("levelset",0,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->updateFromLevelset(levelset);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"FlagGrid::updateFromLevelset" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("FlagGrid::updateFromLevelset",e.what()); return 0; } }    
	void fillGrid(int type=TypeFluid); static PyObject* _W_28 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FlagGrid* pbo = dynamic_cast<FlagGrid*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "FlagGrid::fillGrid" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock; int type = _args.getOpt<int >("type",0,TypeFluid,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->fillGrid(type);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"FlagGrid::fillGrid" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("FlagGrid::fillGrid",e.what()); return 0; } }  public: PbArgs _args;}
#define _C_FlagGrid
;

//! helper to compute grid conversion factor between local coordinates of two grids
inline Vec3 calcGridSizeFactor(Vec3i s1, Vec3i s2) {
	return Vec3( Real(s1[0])/s2[0], Real(s1[1])/s2[1], Real(s1[2])/s2[2] );
}

// prototypes for grid plugins
void copyMacToVec3(MACGrid &source, Grid<Vec3>& target);
void convertMacToVec3(MACGrid &source, Grid<Vec3>& target);
void resampleVec3ToMac(Grid<Vec3>& source, MACGrid &target);
void resampleMacToVec3 (MACGrid &source, Grid<Vec3>& target );
void getComponent(const Grid<Vec3>& src, Grid<Real>& dst, int c);
void setComponent(const Grid<Real>& src, Grid<Vec3>& dst, int c);

//******************************************************************************
// enable compilation of a more complicated test data type
// for grids... note - this also enables code parts in fileio.cpp!
// the code below is meant only as an example for a grid with a more complex data type
// and illustrates which functions need to be implemented; it's not needed
// to run any simulations in mantaflow!

#define ENABLE_GRID_TEST_DATATYPE 0

#if ENABLE_GRID_TEST_DATATYPE==1

typedef std::vector<int> nbVectorBaseType;
class nbVector : public nbVectorBaseType {
	public: 
		inline nbVector() : nbVectorBaseType() {};
		inline ~nbVector() {};

		// grid operators require certain functions
		inline nbVector(Real v) : nbVectorBaseType() { this->push_back( (int)v ); };

		inline const nbVector& operator+= ( const nbVector &v1 ) {
			assertMsg(false,"Never call!"); return *this; 
		}
		inline const nbVector& operator-= ( const nbVector &v1 ) {
			assertMsg(false,"Never call!"); return *this; 
		}
		inline const nbVector& operator*= ( const nbVector &v1 ) {
			assertMsg(false,"Never call!"); return *this; 
		}
};

template<> inline nbVector* FluidSolver::getGridPointer<nbVector>() {
	return new nbVector[mGridSize.x * mGridSize.y * mGridSize.z];
}
template<> inline void FluidSolver::freeGridPointer<nbVector>(nbVector* ptr) {
	return delete[] ptr;
}

inline nbVector operator+ ( const nbVector &v1, const nbVector &v2 ) {
	assertMsg(false,"Never call!"); return nbVector(); 
}
inline nbVector operator- ( const nbVector &v1, const nbVector &v2 ) {
	assertMsg(false,"Never call!"); return nbVector(); 
}
inline nbVector operator* ( const nbVector &v1, const nbVector &v2 ) {
	assertMsg(false,"Never call!"); return nbVector(); 
}
template<class S>
inline nbVector operator* ( const nbVector& v, S s ) {
	assertMsg(false,"Never call!"); return nbVector(); 
} 
template<class S> 
inline nbVector operator* ( S s, const nbVector& v ) {
	assertMsg(false,"Never call!"); return nbVector(); 
}

template<> inline nbVector safeDivide<nbVector>(const nbVector &a, const nbVector& b) { 
	assertMsg(false,"Never call!"); return nbVector(); 
}

std::ostream& operator<< ( std::ostream& os, const nbVectorBaseType& i ) {
	os << " nbVectorBaseType NYI ";
	return os;
}

// make data type known to python
// (python keyword changed here, because the preprocessor does not yet parse #ifdefs correctly)
PY THON alias Grid<nbVector> TestDataGrid;
// ? PY THON alias nbVector TestDatatype;

#endif // end ENABLE_GRID_TEST_DATATYPE



//******************************************************************************
// Implementation of inline functions

inline void GridBase::checkIndex(int i, int j, int k) const {
	//if (i<0 || j<0  || i>=mSize.x || j>=mSize.y || (is3D() && (k<0|| k>= mSize.z))) {
	if (i<0 || j<0  || i>=mSize.x || j>=mSize.y || k<0|| k>= mSize.z ) {
		std::ostringstream s;
		s << "Grid " << mName << " dim " << mSize << " : index " << i << "," << j << "," << k << " out of bound ";
		errMsg(s.str());
	}
}

inline void GridBase::checkIndex(IndexInt idx) const {
	if (idx<0 || idx >= mSize.x * mSize.y * mSize.z) {
		std::ostringstream s;
		s << "Grid " << mName << " dim " << mSize << " : index " << idx << " out of bound ";
		errMsg(s.str());
	}
}

bool GridBase::isInBounds(const Vec3i& p) const { 
	return (p.x >= 0 && p.y >= 0 && p.z >= 0 && p.x < mSize.x && p.y < mSize.y && p.z < mSize.z); 
}

bool GridBase::isInBounds(const Vec3i& p, int bnd) const { 
	bool ret = (p.x >= bnd && p.y >= bnd && p.x < mSize.x-bnd && p.y < mSize.y-bnd);
	if(this->is3D()) {
		ret &= (p.z >= bnd && p.z < mSize.z-bnd); 
	} else {
		ret &= (p.z == 0);
	}
	return ret;
}
//! Check if linear index is in the range of the array
bool GridBase::isInBounds(IndexInt idx) const {
	if (idx<0 || idx >= mSize.x * mSize.y * mSize.z) {
		return false;
	}
	return true;
}

inline Vec3 MACGrid::getCentered(int i, int j, int k) const {
	DEBUG_ONLY(checkIndex(i+1,j+1,k));
	const IndexInt idx = index(i,j,k);
	Vec3 v = Vec3(0.5* (mData[idx].x + mData[idx+1].x),
				  0.5* (mData[idx].y + mData[idx+mSize.x].y),
				  0.);
	if( this->is3D() ) {
		DEBUG_ONLY(checkIndex(idx+mStrideZ));
		v[2] =    0.5* (mData[idx].z + mData[idx+mStrideZ].z);
	}
	return v;
}

inline Vec3 MACGrid::getAtMACX(int i, int j, int k) const {
	DEBUG_ONLY(checkIndex(i-1,j+1,k));
	const IndexInt idx = index(i,j,k);
	Vec3 v =  Vec3(   (mData[idx].x),
				0.25* (mData[idx].y + mData[idx-1].y + mData[idx+mSize.x].y + mData[idx+mSize.x-1].y),
				0.);
	if( this->is3D() ) {
		DEBUG_ONLY(checkIndex(idx+mStrideZ-1));
		v[2] = 0.25* (mData[idx].z + mData[idx-1].z + mData[idx+mStrideZ].z + mData[idx+mStrideZ-1].z);
	}
	return v;
}

inline Vec3 MACGrid::getAtMACY(int i, int j, int k) const {
	DEBUG_ONLY(checkIndex(i+1,j-1,k));
	const IndexInt idx = index(i,j,k);
	Vec3 v =  Vec3(0.25* (mData[idx].x + mData[idx-mSize.x].x + mData[idx+1].x + mData[idx+1-mSize.x].x),
						 (mData[idx].y),   0. );
	if( this->is3D() ) {
		DEBUG_ONLY(checkIndex(idx+mStrideZ-mSize.x));
		v[2] = 0.25* (mData[idx].z + mData[idx-mSize.x].z + mData[idx+mStrideZ].z + mData[idx+mStrideZ-mSize.x].z);
	}
	return v;
}

inline Vec3 MACGrid::getAtMACZ(int i, int j, int k) const {
	const IndexInt idx = index(i,j,k);
	DEBUG_ONLY(checkIndex(idx-mStrideZ));
	DEBUG_ONLY(checkIndex(idx+mSize.x-mStrideZ));
	Vec3 v =  Vec3(0.25* (mData[idx].x + mData[idx-mStrideZ].x + mData[idx+1].x + mData[idx+1-mStrideZ].x),
				   0.25* (mData[idx].y + mData[idx-mStrideZ].y + mData[idx+mSize.x].y + mData[idx+mSize.x-mStrideZ].y),
						 (mData[idx].z) );
	return v;
}

template <class T, class S>  struct gridAdd : public KernelBase { gridAdd(Grid<T>& me, const Grid<S>& other) :  KernelBase(&me,0) ,me(me),other(other)   { runMessage(); run(); }   inline void op(IndexInt idx, Grid<T>& me, const Grid<S>& other )  { me[idx] += other[idx]; }    inline Grid<T>& getArg0() { return me; } typedef Grid<T> type0;inline const Grid<S>& getArg1() { return other; } typedef Grid<S> type1; void runMessage() { debMsg("Executing kernel gridAdd ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (IndexInt i = 0; i < _sz; i++) op(i,me,other);  }   } Grid<T>& me; const Grid<S>& other;   };
#line 507 "grid.h"


template <class T, class S>  struct gridSub : public KernelBase { gridSub(Grid<T>& me, const Grid<S>& other) :  KernelBase(&me,0) ,me(me),other(other)   { runMessage(); run(); }   inline void op(IndexInt idx, Grid<T>& me, const Grid<S>& other )  { me[idx] -= other[idx]; }    inline Grid<T>& getArg0() { return me; } typedef Grid<T> type0;inline const Grid<S>& getArg1() { return other; } typedef Grid<S> type1; void runMessage() { debMsg("Executing kernel gridSub ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (IndexInt i = 0; i < _sz; i++) op(i,me,other);  }   } Grid<T>& me; const Grid<S>& other;   };
#line 508 "grid.h"


template <class T, class S>  struct gridMult : public KernelBase { gridMult(Grid<T>& me, const Grid<S>& other) :  KernelBase(&me,0) ,me(me),other(other)   { runMessage(); run(); }   inline void op(IndexInt idx, Grid<T>& me, const Grid<S>& other )  { me[idx] *= other[idx]; }    inline Grid<T>& getArg0() { return me; } typedef Grid<T> type0;inline const Grid<S>& getArg1() { return other; } typedef Grid<S> type1; void runMessage() { debMsg("Executing kernel gridMult ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (IndexInt i = 0; i < _sz; i++) op(i,me,other);  }   } Grid<T>& me; const Grid<S>& other;   };
#line 509 "grid.h"


template <class T, class S>  struct gridDiv : public KernelBase { gridDiv(Grid<T>& me, const Grid<S>& other) :  KernelBase(&me,0) ,me(me),other(other)   { runMessage(); run(); }   inline void op(IndexInt idx, Grid<T>& me, const Grid<S>& other )  { me[idx] /= other[idx]; }    inline Grid<T>& getArg0() { return me; } typedef Grid<T> type0;inline const Grid<S>& getArg1() { return other; } typedef Grid<S> type1; void runMessage() { debMsg("Executing kernel gridDiv ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (IndexInt i = 0; i < _sz; i++) op(i,me,other);  }   } Grid<T>& me; const Grid<S>& other;   };
#line 510 "grid.h"


template <class T, class S>  struct gridAddScalar : public KernelBase { gridAddScalar(Grid<T>& me, const S& other) :  KernelBase(&me,0) ,me(me),other(other)   { runMessage(); run(); }   inline void op(IndexInt idx, Grid<T>& me, const S& other )  { me[idx] += other; }    inline Grid<T>& getArg0() { return me; } typedef Grid<T> type0;inline const S& getArg1() { return other; } typedef S type1; void runMessage() { debMsg("Executing kernel gridAddScalar ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (IndexInt i = 0; i < _sz; i++) op(i,me,other);  }   } Grid<T>& me; const S& other;   };
#line 511 "grid.h"


template <class T, class S>  struct gridMultScalar : public KernelBase { gridMultScalar(Grid<T>& me, const S& other) :  KernelBase(&me,0) ,me(me),other(other)   { runMessage(); run(); }   inline void op(IndexInt idx, Grid<T>& me, const S& other )  { me[idx] *= other; }    inline Grid<T>& getArg0() { return me; } typedef Grid<T> type0;inline const S& getArg1() { return other; } typedef S type1; void runMessage() { debMsg("Executing kernel gridMultScalar ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (IndexInt i = 0; i < _sz; i++) op(i,me,other);  }   } Grid<T>& me; const S& other;   };
#line 512 "grid.h"


template <class T, class S>  struct gridScaledAdd : public KernelBase { gridScaledAdd(Grid<T>& me, const Grid<T>& other, const S& factor) :  KernelBase(&me,0) ,me(me),other(other),factor(factor)   { runMessage(); run(); }   inline void op(IndexInt idx, Grid<T>& me, const Grid<T>& other, const S& factor )  { me[idx] += factor * other[idx]; }    inline Grid<T>& getArg0() { return me; } typedef Grid<T> type0;inline const Grid<T>& getArg1() { return other; } typedef Grid<T> type1;inline const S& getArg2() { return factor; } typedef S type2; void runMessage() { debMsg("Executing kernel gridScaledAdd ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (IndexInt i = 0; i < _sz; i++) op(i,me,other,factor);  }   } Grid<T>& me; const Grid<T>& other; const S& factor;   };
#line 513 "grid.h"



template <class T>  struct gridSafeDiv : public KernelBase { gridSafeDiv(Grid<T>& me, const Grid<T>& other) :  KernelBase(&me,0) ,me(me),other(other)   { runMessage(); run(); }   inline void op(IndexInt idx, Grid<T>& me, const Grid<T>& other )  { me[idx] = safeDivide(me[idx], other[idx]); }    inline Grid<T>& getArg0() { return me; } typedef Grid<T> type0;inline const Grid<T>& getArg1() { return other; } typedef Grid<T> type1; void runMessage() { debMsg("Executing kernel gridSafeDiv ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (IndexInt i = 0; i < _sz; i++) op(i,me,other);  }   } Grid<T>& me; const Grid<T>& other;   };
#line 515 "grid.h"


template <class T>  struct gridSetConst : public KernelBase { gridSetConst(Grid<T>& grid, T value) :  KernelBase(&grid,0) ,grid(grid),value(value)   { runMessage(); run(); }   inline void op(IndexInt idx, Grid<T>& grid, T value )  { grid[idx] = value; }    inline Grid<T>& getArg0() { return grid; } typedef Grid<T> type0;inline T& getArg1() { return value; } typedef T type1; void runMessage() { debMsg("Executing kernel gridSetConst ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (IndexInt i = 0; i < _sz; i++) op(i,grid,value);  }   } Grid<T>& grid; T value;   };
#line 516 "grid.h"



template<class T> template<class S> Grid<T>& Grid<T>::operator+= (const Grid<S>& a) {
	gridAdd<T,S> (*this, a);
	return *this;
}
template<class T> template<class S> Grid<T>& Grid<T>::operator+= (const S& a) {
	gridAddScalar<T,S> (*this, a);
	return *this;
}
template<class T> template<class S> Grid<T>& Grid<T>::operator-= (const Grid<S>& a) {
	gridSub<T,S> (*this, a);
	return *this;
}
template<class T> template<class S> Grid<T>& Grid<T>::operator-= (const S& a) {
	gridAddScalar<T,S> (*this, -a);
	return *this;
}
template<class T> template<class S> Grid<T>& Grid<T>::operator*= (const Grid<S>& a) {
	gridMult<T,S> (*this, a);
	return *this;
}
template<class T> template<class S> Grid<T>& Grid<T>::operator*= (const S& a) {
	gridMultScalar<T,S> (*this, a);
	return *this;
}
template<class T> template<class S> Grid<T>& Grid<T>::operator/= (const Grid<S>& a) {
	gridDiv<T,S> (*this, a);
	return *this;
}
template<class T> template<class S> Grid<T>& Grid<T>::operator/= (const S& a) {
	S rez((S)1.0 / a);
	gridMultScalar<T,S> (*this, rez);
	return *this;
}


//******************************************************************************
// Other helper functions

// compute gradient of a scalar grid
inline Vec3 getGradient(const Grid<Real>& data, int i, int j, int k) {
	Vec3 v;

	if (i > data.getSizeX()-2) i= data.getSizeX()-2;
	if (j > data.getSizeY()-2) j= data.getSizeY()-2;
	if (i < 1) i = 1;
	if (j < 1) j = 1;
	v = Vec3( data(i+1,j  ,k  ) - data(i-1,j  ,k  ) ,
			  data(i  ,j+1,k  ) - data(i  ,j-1,k  ) , 0. );

	if(data.is3D()) {
		if (k > data.getSizeZ()-2) k= data.getSizeZ()-2;
		if (k < 1) k = 1;
		v[2]= data(i  ,j  ,k+1) - data(i  ,j  ,k-1);
	} 

	return v;
}

// interpolate grid from one size to another size

template <class S>  struct knInterpolateGridTempl : public KernelBase { knInterpolateGridTempl(Grid<S>& target, Grid<S>& source, const Vec3& sourceFactor , Vec3 offset, int orderSpace=1 ) :  KernelBase(&target,0) ,target(target),source(source),sourceFactor(sourceFactor),offset(offset),orderSpace(orderSpace)   { runMessage(); run(); }  inline void op(int i, int j, int k, Grid<S>& target, Grid<S>& source, const Vec3& sourceFactor , Vec3 offset, int orderSpace=1  )  {
	Vec3 pos = Vec3(i,j,k) * sourceFactor + offset;
	if(!source.is3D()) pos[2] = 0; // allow 2d -> 3d
	target(i,j,k) = source.getInterpolatedHi(pos, orderSpace);
}   inline Grid<S>& getArg0() { return target; } typedef Grid<S> type0;inline Grid<S>& getArg1() { return source; } typedef Grid<S> type1;inline const Vec3& getArg2() { return sourceFactor; } typedef Vec3 type2;inline Vec3& getArg3() { return offset; } typedef Vec3 type3;inline int& getArg4() { return orderSpace; } typedef int type4; void runMessage() { debMsg("Executing kernel knInterpolateGridTempl ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int k=minZ; k < maxZ; k++) for (int j=0; j < _maxY; j++) for (int i=0; i < _maxX; i++) op(i,j,k,target,source,sourceFactor,offset,orderSpace);  } } else { const int k=0; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int j=0; j < _maxY; j++) for (int i=0; i < _maxX; i++) op(i,j,k,target,source,sourceFactor,offset,orderSpace);  } }  } Grid<S>& target; Grid<S>& source; const Vec3& sourceFactor; Vec3 offset; int orderSpace;   };
#line 578 "grid.h"

 
// template glue code - choose interpolation based on template arguments
template<class GRID>
void interpolGridTempl( GRID& target, GRID& source ) {
	errMsg("interpolGridTempl - Only valid for specific instantiations");
}


} //namespace
#endif


