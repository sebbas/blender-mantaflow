




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
 * Set boundary conditions, gravity
 *
 ******************************************************************************/

#include "vectorbase.h"
#include "grid.h"
#include "commonkernels.h"

using namespace std;

namespace Manta { 

// MLE 2014-07-05 copy from pressure.cpp
inline void convertDescToVec(const string& desc, Vector3D<bool>& lo, Vector3D<bool>& up) {
    for(size_t i=0; i<desc.size(); i++) {
        if (desc[i] == 'x') lo.x = true;
        else if (desc[i] == 'y') lo.y = true;
        else if (desc[i] == 'z') lo.z = true;
        else if (desc[i] == 'X') up.x = true;
        else if (desc[i] == 'Y') up.y = true;
        else if (desc[i] == 'Z') up.z = true;
        else errMsg("invalid character in boundary description string. Only [xyzXYZ] allowed.");
    }
}

//! add Forces between fl/fl and fl/em cells
 struct KnAddForceField : public KernelBase { KnAddForceField(FlagGrid& flags, MACGrid& vel, Grid<Vec3>& force) :  KernelBase(&flags,1) ,flags(flags),vel(vel),force(force)   { run(); }  inline void op(int i, int j, int k, FlagGrid& flags, MACGrid& vel, Grid<Vec3>& force )  {
	bool curFluid = flags.isFluid(i,j,k);
	bool curEmpty = flags.isEmpty(i,j,k);
	if (!curFluid && !curEmpty) return;
	
	if (flags.isFluid(i-1,j,k) || (curFluid && flags.isEmpty(i-1,j,k))) 
		vel(i,j,k).x += 0.5*(force(i-1,j,k).x + force(i,j,k).x);
	if (flags.isFluid(i,j-1,k) || (curFluid && flags.isEmpty(i,j-1,k))) 
		vel(i,j,k).y += 0.5*(force(i,j-1,k).y + force(i,j,k).y);
	if (vel.is3D() && (flags.isFluid(i,j,k-1) || (curFluid && flags.isEmpty(i,j,k-1))))
		vel(i,j,k).z += 0.5*(force(i,j,k-1).z + force(i,j,k).z);
}   inline FlagGrid& getArg0() { return flags; } typedef FlagGrid type0;inline MACGrid& getArg1() { return vel; } typedef MACGrid type1;inline Grid<Vec3>& getArg2() { return force; } typedef Grid<Vec3> type2; void run() {  const int _maxX = maxX; const int _maxY = maxY; for (int k=minZ; k< maxZ; k++) for (int j=1; j< _maxY; j++) for (int i=1; i< _maxX; i++) op(i,j,k, flags,vel,force);  } FlagGrid& flags; MACGrid& vel; Grid<Vec3>& force;   };

//! add Forces between fl/fl and fl/em cells
 struct KnAddForce : public KernelBase { KnAddForce(FlagGrid& flags, MACGrid& vel, Vec3 force) :  KernelBase(&flags,1) ,flags(flags),vel(vel),force(force)   { run(); }  inline void op(int i, int j, int k, FlagGrid& flags, MACGrid& vel, Vec3 force )  {
	bool curFluid = flags.isFluid(i,j,k);
	bool curEmpty = flags.isEmpty(i,j,k);
	if (!curFluid && !curEmpty) return;
	
	if (flags.isFluid(i-1,j,k) || (curFluid && flags.isEmpty(i-1,j,k))) 
		vel(i,j,k).x += force.x;
	if (flags.isFluid(i,j-1,k) || (curFluid && flags.isEmpty(i,j-1,k))) 
		vel(i,j,k).y += force.y;
	if (vel.is3D() && (flags.isFluid(i,j,k-1) || (curFluid && flags.isEmpty(i,j,k-1))))
		vel(i,j,k).z += force.z;
}   inline FlagGrid& getArg0() { return flags; } typedef FlagGrid type0;inline MACGrid& getArg1() { return vel; } typedef MACGrid type1;inline Vec3& getArg2() { return force; } typedef Vec3 type2; void run() {  const int _maxX = maxX; const int _maxY = maxY; for (int k=minZ; k< maxZ; k++) for (int j=1; j< _maxY; j++) for (int i=1; i< _maxX; i++) op(i,j,k, flags,vel,force);  } FlagGrid& flags; MACGrid& vel; Vec3 force;   };

//! add external force fields to all fluid cells
void addForceField(FlagGrid& flags, MACGrid& vel, Grid<Vec3>& force){    
	KnAddForceField(flags, vel, force);
} static PyObject* _W_0 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "addForceField" ); PyObject *_retval = 0; { ArgLocker _lock; FlagGrid& flags = *_args.getPtr<FlagGrid >("flags",0,&_lock); MACGrid& vel = *_args.getPtr<MACGrid >("vel",1,&_lock); Grid<Vec3>& force = *_args.getPtr<Grid<Vec3> >("force",2,&_lock);   _retval = getPyNone(); addForceField(flags,vel,force);  _args.check(); } pbFinalizePlugin(parent,"addForceField" ); return _retval; } catch(std::exception& e) { pbSetError("addForceField",e.what()); return 0; } } static const Pb::Register _RP_addForceField ("","addForceField",_W_0); 
	
//! add gravity forces to all fluid cells
void addGravity(FlagGrid& flags, MACGrid& vel, Vec3 gravity) {    
	Vec3 f = gravity * flags.getParent()->getDt() / flags.getDx();
	KnAddForce(flags, vel, f);
} static PyObject* _W_1 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "addGravity" ); PyObject *_retval = 0; { ArgLocker _lock; FlagGrid& flags = *_args.getPtr<FlagGrid >("flags",0,&_lock); MACGrid& vel = *_args.getPtr<MACGrid >("vel",1,&_lock); Vec3 gravity = _args.get<Vec3 >("gravity",2,&_lock);   _retval = getPyNone(); addGravity(flags,vel,gravity);  _args.check(); } pbFinalizePlugin(parent,"addGravity" ); return _retval; } catch(std::exception& e) { pbSetError("addGravity",e.what()); return 0; } } static const Pb::Register _RP_addGravity ("","addGravity",_W_1); 

//! add Buoyancy force based on smoke density
 struct KnAddBuoyancy : public KernelBase { KnAddBuoyancy(FlagGrid& flags, Grid<Real>& density, MACGrid& vel, Vec3 strength) :  KernelBase(&flags,1) ,flags(flags),density(density),vel(vel),strength(strength)   { run(); }  inline void op(int i, int j, int k, FlagGrid& flags, Grid<Real>& density, MACGrid& vel, Vec3 strength )  {    
	if (!flags.isFluid(i,j,k)) return;
	if (flags.isFluid(i-1,j,k))
		vel(i,j,k).x += (0.5 * strength.x) * (density(i,j,k)+density(i-1,j,k));
	if (flags.isFluid(i,j-1,k))
		vel(i,j,k).y += (0.5 * strength.y) * (density(i,j,k)+density(i,j-1,k));
	if (vel.is3D() && flags.isFluid(i,j,k-1))
		vel(i,j,k).z += (0.5 * strength.z) * (density(i,j,k)+density(i,j,k-1));    
}   inline FlagGrid& getArg0() { return flags; } typedef FlagGrid type0;inline Grid<Real>& getArg1() { return density; } typedef Grid<Real> type1;inline MACGrid& getArg2() { return vel; } typedef MACGrid type2;inline Vec3& getArg3() { return strength; } typedef Vec3 type3; void run() {  const int _maxX = maxX; const int _maxY = maxY; for (int k=minZ; k< maxZ; k++) for (int j=1; j< _maxY; j++) for (int i=1; i< _maxX; i++) op(i,j,k, flags,density,vel,strength);  } FlagGrid& flags; Grid<Real>& density; MACGrid& vel; Vec3 strength;   };

//! add Buoyancy force based on smoke density
void addBuoyancy(FlagGrid& flags, Grid<Real>& density, MACGrid& vel, Vec3 gravity) {
	Vec3 f = - gravity * flags.getParent()->getDt() / flags.getParent()->getDx();
	KnAddBuoyancy(flags,density, vel, f);
} static PyObject* _W_2 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "addBuoyancy" ); PyObject *_retval = 0; { ArgLocker _lock; FlagGrid& flags = *_args.getPtr<FlagGrid >("flags",0,&_lock); Grid<Real>& density = *_args.getPtr<Grid<Real> >("density",1,&_lock); MACGrid& vel = *_args.getPtr<MACGrid >("vel",2,&_lock); Vec3 gravity = _args.get<Vec3 >("gravity",3,&_lock);   _retval = getPyNone(); addBuoyancy(flags,density,vel,gravity);  _args.check(); } pbFinalizePlugin(parent,"addBuoyancy" ); return _retval; } catch(std::exception& e) { pbSetError("addBuoyancy",e.what()); return 0; } } static const Pb::Register _RP_addBuoyancy ("","addBuoyancy",_W_2); 

		
//! set no-stick wall boundary condition between ob/fl and ob/ob cells
 struct KnSetWallBcs : public KernelBase { KnSetWallBcs(FlagGrid& flags, MACGrid& vel, Vector3D<bool> lo, Vector3D<bool> up, bool admm) :  KernelBase(&flags,0) ,flags(flags),vel(vel),lo(lo),up(up),admm(admm)   { run(); }  inline void op(int i, int j, int k, FlagGrid& flags, MACGrid& vel, Vector3D<bool> lo, Vector3D<bool> up, bool admm )  {

	bool curFluid = flags.isFluid(i,j,k);
    bool curObstacle = flags.isObstacle(i,j,k);
	if (!curFluid && !curObstacle) return;

	// MLE 2014-07-04
	// if not admm, leave it as in orig
	// if openBound, don't correct anything (solid is as empty)
	// if admm, correct if vel is pointing outwards
	
	// if "inner" obstacle vel
	if(i>0 && curObstacle && !flags.isFluid(i-1,j,k)) vel(i,j,k).x = 0;
	if(j>0 && curObstacle && !flags.isFluid(i,j-1,k)) vel(i,j,k).y = 0;

	// check lo.x
	if(!lo.x && i>0 && curFluid && flags.isObstacle(i-1,j,k) && ((admm&&vel(i,j,k).x<0)||!admm)) vel(i,j,k).x = 0;
	// check up.x
	if(!up.x && i>0 && curObstacle && flags.isFluid(i-1,j,k) && ((admm&&vel(i,j,k).x>0)||!admm)) vel(i,j,k).x = 0;
	// check lo.y
	if(!lo.y && j>0 && curFluid && flags.isObstacle(i,j-1,k) && ((admm&&vel(i,j,k).y<0)||!admm)) vel(i,j,k).y = 0;
	// check up.y
	if(!up.y && j>0 && curObstacle && flags.isFluid(i,j-1,k) && ((admm&&vel(i,j,k).y>0)||!admm)) vel(i,j,k).y = 0;
	// check lo.z
	if(!lo.z && k>0 && curFluid && flags.isObstacle(i,j,k-1) && ((admm&&vel(i,j,k).z<0)||!admm)) vel(i,j,k).z = 0;
	// check up.z
	if(!up.z && k>0 && curObstacle && flags.isFluid(i,j,k-1) && ((admm&&vel(i,j,k).z>0)||!admm)) vel(i,j,k).z = 0;
	

	/* MLE consider later	
	if (curFluid) {
		if ((i>0 && flags.isStick(i-1,j,k)) || (i<flags.getSizeX()-1 && flags.isStick(i+1,j,k)))
			vel(i,j,k).y = vel(i,j,k).z = 0;
		if ((j>0 && flags.isStick(i,j-1,k)) || (j<flags.getSizeY()-1 && flags.isStick(i,j+1,k)))
			vel(i,j,k).x = vel(i,j,k).z = 0;
		if (vel.is3D() && ((k>0 && flags.isStick(i,j,k-1)) || (k<flags.getSizeZ()-1 && flags.isStick(i,j,k+1))))
			vel(i,j,k).x = vel(i,j,k).y = 0;
	}
	*/
}   inline FlagGrid& getArg0() { return flags; } typedef FlagGrid type0;inline MACGrid& getArg1() { return vel; } typedef MACGrid type1;inline Vector3D<bool> & getArg2() { return lo; } typedef Vector3D<bool>  type2;inline Vector3D<bool> & getArg3() { return up; } typedef Vector3D<bool>  type3;inline bool& getArg4() { return admm; } typedef bool type4; void run() {  const int _maxX = maxX; const int _maxY = maxY; for (int k=minZ; k< maxZ; k++) for (int j=0; j< _maxY; j++) for (int i=0; i< _maxX; i++) op(i,j,k, flags,vel,lo,up,admm);  } FlagGrid& flags; MACGrid& vel; Vector3D<bool>  lo; Vector3D<bool>  up; bool admm;   };

// MLE 2014-07-04
//! set no-stick boundary condition on walls
void setWallBcs(FlagGrid& flags, MACGrid& vel, string openBound="", bool admm=false) {
	Vector3D<bool> lo, up;
    convertDescToVec(openBound, lo, up);
    KnSetWallBcs(flags, vel, lo, up, admm);
} static PyObject* _W_3 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "setWallBcs" ); PyObject *_retval = 0; { ArgLocker _lock; FlagGrid& flags = *_args.getPtr<FlagGrid >("flags",0,&_lock); MACGrid& vel = *_args.getPtr<MACGrid >("vel",1,&_lock); string openBound = _args.getOpt<string >("openBound",2,"",&_lock); bool admm = _args.getOpt<bool >("admm",3,false,&_lock);   _retval = getPyNone(); setWallBcs(flags,vel,openBound,admm);  _args.check(); } pbFinalizePlugin(parent,"setWallBcs" ); return _retval; } catch(std::exception& e) { pbSetError("setWallBcs",e.what()); return 0; } } static const Pb::Register _RP_setWallBcs ("","setWallBcs",_W_3);  
//! Kernel: gradient norm operator
 struct KnConfForce : public KernelBase { KnConfForce(Grid<Vec3>& force, const Grid<Real>& grid, const Grid<Vec3>& curl, Real str) :  KernelBase(&force,1) ,force(force),grid(grid),curl(curl),str(str)   { run(); }  inline void op(int i, int j, int k, Grid<Vec3>& force, const Grid<Real>& grid, const Grid<Vec3>& curl, Real str )  {
	Vec3 grad = 0.5 * Vec3(        grid(i+1,j,k)-grid(i-1,j,k), 
								   grid(i,j+1,k)-grid(i,j-1,k), 0.);
	if(grid.is3D()) grad[2]= 0.5*( grid(i,j,k+1)-grid(i,j,k-1) );
	normalize(grad);
	force(i,j,k) = str * cross(grad, curl(i,j,k));
}   inline Grid<Vec3>& getArg0() { return force; } typedef Grid<Vec3> type0;inline const Grid<Real>& getArg1() { return grid; } typedef Grid<Real> type1;inline const Grid<Vec3>& getArg2() { return curl; } typedef Grid<Vec3> type2;inline Real& getArg3() { return str; } typedef Real type3; void run() {  const int _maxX = maxX; const int _maxY = maxY; for (int k=minZ; k< maxZ; k++) for (int j=1; j< _maxY; j++) for (int i=1; i< _maxX; i++) op(i,j,k, force,grid,curl,str);  } Grid<Vec3>& force; const Grid<Real>& grid; const Grid<Vec3>& curl; Real str;   };

void vorticityConfinement(MACGrid& vel, FlagGrid& flags, Real strength) {
	Grid<Vec3> velCenter(flags.getParent()), curl(flags.getParent()), force(flags.getParent());
	Grid<Real> norm(flags.getParent());
	
	GetCentered(velCenter, vel);
	CurlOp(velCenter, curl);
	GridNorm(norm, curl);
	KnConfForce(force, norm, curl, strength);
	KnAddForceField(flags, vel, force);
} static PyObject* _W_4 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "vorticityConfinement" ); PyObject *_retval = 0; { ArgLocker _lock; MACGrid& vel = *_args.getPtr<MACGrid >("vel",0,&_lock); FlagGrid& flags = *_args.getPtr<FlagGrid >("flags",1,&_lock); Real strength = _args.get<Real >("strength",2,&_lock);   _retval = getPyNone(); vorticityConfinement(vel,flags,strength);  _args.check(); } pbFinalizePlugin(parent,"vorticityConfinement" ); return _retval; } catch(std::exception& e) { pbSetError("vorticityConfinement",e.what()); return 0; } } static const Pb::Register _RP_vorticityConfinement ("","vorticityConfinement",_W_4); 

//! enforce a constant inflow/outflow at the grid boundaries
 struct KnSetInflow : public KernelBase { KnSetInflow(MACGrid& vel, int dim, int p0, const Vec3& val) :  KernelBase(&vel,0) ,vel(vel),dim(dim),p0(p0),val(val)   { run(); }  inline void op(int i, int j, int k, MACGrid& vel, int dim, int p0, const Vec3& val )  {
	Vec3i p(i,j,k);
	if (p[dim] == p0 || p[dim] == p0+1)
		vel(i,j,k) = val;
}   inline MACGrid& getArg0() { return vel; } typedef MACGrid type0;inline int& getArg1() { return dim; } typedef int type1;inline int& getArg2() { return p0; } typedef int type2;inline const Vec3& getArg3() { return val; } typedef Vec3 type3; void run() {  const int _maxX = maxX; const int _maxY = maxY; for (int k=minZ; k< maxZ; k++) for (int j=0; j< _maxY; j++) for (int i=0; i< _maxX; i++) op(i,j,k, vel,dim,p0,val);  } MACGrid& vel; int dim; int p0; const Vec3& val;   };

//! enforce a constant inflow/outflow at the grid boundaries
void setInflowBcs(MACGrid& vel, string dir, Vec3 value) {
	for(size_t i=0; i<dir.size(); i++) {
		if (dir[i] >= 'x' && dir[i] <= 'z') { 
			int dim = dir[i]-'x';
			KnSetInflow(vel,dim,0,value);
		} else if (dir[i] >= 'X' && dir[i] <= 'Z') {
			int dim = dir[i]-'X';
			KnSetInflow(vel,dim,vel.getSize()[dim]-1,value);
		} else 
			errMsg("invalid character in direction string. Only [xyzXYZ] allowed.");
	}
} static PyObject* _W_5 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "setInflowBcs" ); PyObject *_retval = 0; { ArgLocker _lock; MACGrid& vel = *_args.getPtr<MACGrid >("vel",0,&_lock); string dir = _args.get<string >("dir",1,&_lock); Vec3 value = _args.get<Vec3 >("value",2,&_lock);   _retval = getPyNone(); setInflowBcs(vel,dir,value);  _args.check(); } pbFinalizePlugin(parent,"setInflowBcs" ); return _retval; } catch(std::exception& e) { pbSetError("setInflowBcs",e.what()); return 0; } } static const Pb::Register _RP_setInflowBcs ("","setInflowBcs",_W_5); 

} // namespace


