




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * shapes classes
 *
 ******************************************************************************/

#ifndef _SHAPES_H
#define _SHAPES_H

#include "pwrapper/manta.h"
#include "util/vectorbase.h"
#include "levelset.h"

namespace Manta {

// forward declaration
class Mesh;
	
//! Base class for all shapes
class Shape : public PbClass {public:
	enum GridType { TypeNone = 0, TypeBox = 1, TypeSphere = 2, TypeCylinder };
	
	Shape(FluidSolver* parent); static int _W_0 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = Pb::objFromPy(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "Shape::Shape" ); { ArgLocker _lock; FluidSolver* parent = _args.getPtr<FluidSolver >("parent",0,&_lock);  obj = new Shape(parent); obj->registerObject(_self, &_args); _args.check(); } pbFinalizePlugin(obj->getParent(),"Shape::Shape" ); return 0; } catch(std::exception& e) { pbSetError("Shape::Shape",e.what()); return -1; } }
	
	//! Get the type of grid
	inline GridType getType() const { return mType; }
	
	//! Apply shape to flag grid, set inside cells to <value>
	void applyToGrid(GridBase* grid, FlagGrid* respectFlags=0); static PyObject* _W_1 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Shape* pbo = dynamic_cast<Shape*>(Pb::objFromPy(_self)); pbPreparePlugin(pbo->getParent(), "Shape::applyToGrid"); PyObject *_retval = 0; { ArgLocker _lock; GridBase* grid = _args.getPtr<GridBase >("grid",0,&_lock); FlagGrid* respectFlags = _args.getPtrOpt<FlagGrid >("respectFlags",1,0,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->applyToGrid(grid,respectFlags);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Shape::applyToGrid"); return _retval; } catch(std::exception& e) { pbSetError("Shape::applyToGrid",e.what()); return 0; } }
	void applyToGridSmooth(GridBase* grid, Real sigma=1.0, Real shift=0, FlagGrid* respectFlags=0); static PyObject* _W_2 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Shape* pbo = dynamic_cast<Shape*>(Pb::objFromPy(_self)); pbPreparePlugin(pbo->getParent(), "Shape::applyToGridSmooth"); PyObject *_retval = 0; { ArgLocker _lock; GridBase* grid = _args.getPtr<GridBase >("grid",0,&_lock); Real sigma = _args.getOpt<Real >("sigma",1,1.0,&_lock); Real shift = _args.getOpt<Real >("shift",2,0,&_lock); FlagGrid* respectFlags = _args.getPtrOpt<FlagGrid >("respectFlags",3,0,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->applyToGridSmooth(grid,sigma,shift,respectFlags);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Shape::applyToGridSmooth"); return _retval; } catch(std::exception& e) { pbSetError("Shape::applyToGridSmooth",e.what()); return 0; } }
	LevelsetGrid computeLevelset(); static PyObject* _W_3 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Shape* pbo = dynamic_cast<Shape*>(Pb::objFromPy(_self)); pbPreparePlugin(pbo->getParent(), "Shape::computeLevelset"); PyObject *_retval = 0; { ArgLocker _lock;  pbo->_args.copy(_args);  _retval = toPy(pbo->computeLevelset());  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Shape::computeLevelset"); return _retval; } catch(std::exception& e) { pbSetError("Shape::computeLevelset",e.what()); return 0; } }
	void collideMesh(Mesh& mesh); static PyObject* _W_4 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Shape* pbo = dynamic_cast<Shape*>(Pb::objFromPy(_self)); pbPreparePlugin(pbo->getParent(), "Shape::collideMesh"); PyObject *_retval = 0; { ArgLocker _lock; Mesh& mesh = *_args.getPtr<Mesh >("mesh",0,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->collideMesh(mesh);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Shape::collideMesh"); return _retval; } catch(std::exception& e) { pbSetError("Shape::collideMesh",e.what()); return 0; } }
	virtual Vec3 getCenter() const { return Vec3::Zero; } static PyObject* _W_5 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Shape* pbo = dynamic_cast<Shape*>(Pb::objFromPy(_self)); pbPreparePlugin(pbo->getParent(), "Shape::getCenter"); PyObject *_retval = 0; { ArgLocker _lock;  pbo->_args.copy(_args);  _retval = toPy(pbo->getCenter());  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Shape::getCenter"); return _retval; } catch(std::exception& e) { pbSetError("Shape::getCenter",e.what()); return 0; } }
	virtual void setCenter(const Vec3& center) {} static PyObject* _W_6 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Shape* pbo = dynamic_cast<Shape*>(Pb::objFromPy(_self)); pbPreparePlugin(pbo->getParent(), "Shape::setCenter"); PyObject *_retval = 0; { ArgLocker _lock; const Vec3& center = _args.get<Vec3 >("center",0,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->setCenter(center);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Shape::setCenter"); return _retval; } catch(std::exception& e) { pbSetError("Shape::setCenter",e.what()); return 0; } }
	virtual Vec3 getExtent() const { return Vec3::Zero; } static PyObject* _W_7 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Shape* pbo = dynamic_cast<Shape*>(Pb::objFromPy(_self)); pbPreparePlugin(pbo->getParent(), "Shape::getExtent"); PyObject *_retval = 0; { ArgLocker _lock;  pbo->_args.copy(_args);  _retval = toPy(pbo->getExtent());  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Shape::getExtent"); return _retval; } catch(std::exception& e) { pbSetError("Shape::getExtent",e.what()); return 0; } }
	
	//! Inside test of the shape
	virtual bool isInside(const Vec3& pos) const;
	inline bool isInsideGrid(int i, int j, int k) const { return isInside(Vec3(i+0.5,j+0.5,k+0.5)); };
	
	virtual void generateMesh(Mesh* mesh) {} ;    
	virtual void generateLevelset(Grid<Real>& phi) {};    
	
protected: 	GridType mType; public: PbArgs _args;}
#define _C_Shape
;

//! Dummy shape
class NullShape : public Shape {   
public:
	NullShape(FluidSolver* parent) :Shape(parent){} static int _W_8 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = Pb::objFromPy(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "NullShape::NullShape" ); { ArgLocker _lock; FluidSolver* parent = _args.getPtr<FluidSolver >("parent",0,&_lock);  obj = new NullShape(parent); obj->registerObject(_self, &_args); _args.check(); } pbFinalizePlugin(obj->getParent(),"NullShape::NullShape" ); return 0; } catch(std::exception& e) { pbSetError("NullShape::NullShape",e.what()); return -1; } }
	
	virtual bool isInside(const Vec3& pos) const { return false; }
	virtual void generateMesh(Mesh* mesh) {}
	
protected: 	virtual void generateLevelset(Grid<Real>& phi) { gridSetConst<Real>( phi , 1000.0f ); } public: PbArgs _args;}
#define _C_NullShape
;

//! Box shape
class Box : public Shape {   
public:
	Box(FluidSolver* parent, Vec3 center = Vec3::Invalid, Vec3 p0 = Vec3::Invalid, Vec3 p1 = Vec3::Invalid, Vec3 size = Vec3::Invalid); static int _W_9 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = Pb::objFromPy(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "Box::Box" ); { ArgLocker _lock; FluidSolver* parent = _args.getPtr<FluidSolver >("parent",0,&_lock); Vec3 center = _args.getOpt<Vec3 >("center",1,Vec3::Invalid,&_lock); Vec3 p0 = _args.getOpt<Vec3 >("p0",2,Vec3::Invalid,&_lock); Vec3 p1 = _args.getOpt<Vec3 >("p1",3,Vec3::Invalid,&_lock); Vec3 size = _args.getOpt<Vec3 >("size",4,Vec3::Invalid,&_lock);  obj = new Box(parent,center,p0,p1,size); obj->registerObject(_self, &_args); _args.check(); } pbFinalizePlugin(obj->getParent(),"Box::Box" ); return 0; } catch(std::exception& e) { pbSetError("Box::Box",e.what()); return -1; } }
	
	inline Vec3 getSize() const { return mP1-mP0; }
	inline Vec3 getP0() const { return mP0; }
	inline Vec3 getP1() const { return mP1; }
	virtual void setCenter(const Vec3& center) { Vec3 dh=0.5*(mP1-mP0); mP0 = center-dh; mP1 = center+dh;}
	virtual Vec3 getCenter() const { return 0.5*(mP1+mP0); }
	virtual Vec3 getExtent() const { return getSize(); }
	virtual bool isInside(const Vec3& pos) const;
	virtual void generateMesh(Mesh* mesh);
	virtual void generateLevelset(Grid<Real>& phi);
	
protected: 	Vec3 mP0, mP1; public: PbArgs _args;}
#define _C_Box
;

//! Spherical shape
class Sphere : public Shape {   
public:
	Sphere(FluidSolver* parent, Vec3 center, Real radius, Vec3 scale=Vec3(1,1,1)); static int _W_10 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = Pb::objFromPy(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "Sphere::Sphere" ); { ArgLocker _lock; FluidSolver* parent = _args.getPtr<FluidSolver >("parent",0,&_lock); Vec3 center = _args.get<Vec3 >("center",1,&_lock); Real radius = _args.get<Real >("radius",2,&_lock); Vec3 scale = _args.getOpt<Vec3 >("scale",3,Vec3(1,1,1),&_lock);  obj = new Sphere(parent,center,radius,scale); obj->registerObject(_self, &_args); _args.check(); } pbFinalizePlugin(obj->getParent(),"Sphere::Sphere" ); return 0; } catch(std::exception& e) { pbSetError("Sphere::Sphere",e.what()); return -1; } }
	
	virtual void setCenter(const Vec3& center) { mCenter = center; }
	virtual Vec3 getCenter() const { return mCenter; }
	inline Real getRadius() const { return mRadius; }
	virtual Vec3 getExtent() const { return Vec3(2.0*mRadius); }    
	virtual bool isInside(const Vec3& pos) const;
	virtual void generateMesh(Mesh* mesh);
	virtual void generateLevelset(Grid<Real>& phi);
	
protected:
	Vec3 mCenter, mScale; 	Real mRadius; public: PbArgs _args;}
#define _C_Sphere
;

//! Cylindrical shape
class Cylinder : public Shape {   
public:
	Cylinder(FluidSolver* parent, Vec3 center, Real radius, Vec3 z); static int _W_11 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = Pb::objFromPy(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "Cylinder::Cylinder" ); { ArgLocker _lock; FluidSolver* parent = _args.getPtr<FluidSolver >("parent",0,&_lock); Vec3 center = _args.get<Vec3 >("center",1,&_lock); Real radius = _args.get<Real >("radius",2,&_lock); Vec3 z = _args.get<Vec3 >("z",3,&_lock);  obj = new Cylinder(parent,center,radius,z); obj->registerObject(_self, &_args); _args.check(); } pbFinalizePlugin(obj->getParent(),"Cylinder::Cylinder" ); return 0; } catch(std::exception& e) { pbSetError("Cylinder::Cylinder",e.what()); return -1; } }
	
	void setRadius(Real r) { mRadius = r; } static PyObject* _W_12 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Cylinder* pbo = dynamic_cast<Cylinder*>(Pb::objFromPy(_self)); pbPreparePlugin(pbo->getParent(), "Cylinder::setRadius"); PyObject *_retval = 0; { ArgLocker _lock; Real r = _args.get<Real >("r",0,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->setRadius(r);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Cylinder::setRadius"); return _retval; } catch(std::exception& e) { pbSetError("Cylinder::setRadius",e.what()); return 0; } }
	void setZ(Vec3 z) { mZDir=z; mZ=normalize(mZDir); } static PyObject* _W_13 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Cylinder* pbo = dynamic_cast<Cylinder*>(Pb::objFromPy(_self)); pbPreparePlugin(pbo->getParent(), "Cylinder::setZ"); PyObject *_retval = 0; { ArgLocker _lock; Vec3 z = _args.get<Vec3 >("z",0,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->setZ(z);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Cylinder::setZ"); return _retval; } catch(std::exception& e) { pbSetError("Cylinder::setZ",e.what()); return 0; } }
	
	virtual void setCenter(const Vec3& center) { mCenter=center; }
	virtual Vec3 getCenter() const { return mCenter; }
	inline Real getRadius() const { return mRadius; }
	inline Vec3 getZ() const { return mZ*mZDir; }
	virtual Vec3 getExtent() const { return Vec3(2.0*sqrt(square(mZ)+square(mRadius))); }    
	virtual bool isInside(const Vec3& pos) const;
	virtual void generateMesh(Mesh* mesh);
	virtual void generateLevelset(Grid<Real>& phi);

protected:
	Vec3 mCenter, mZDir; 	Real mRadius, mZ; public: PbArgs _args;}
#define _C_Cylinder
;

	

} //namespace
#endif


