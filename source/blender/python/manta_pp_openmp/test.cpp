




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/home/user/Developer/mantaflowgit/source/test.cpp"
/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Use this file to test new functionality
 *
 ******************************************************************************/

#include "levelset.h"
#include "commonkernels.h"
#include "particle.h"
#include <cmath>

using namespace std;

namespace Manta {



 struct reductionTest : public KernelBase { reductionTest(const Grid<Real>& v) :  KernelBase(&v,0) ,v(v) ,sum(0)  { run(); }  inline void op(int idx, const Grid<Real>& v ,double& sum)  {
	sum += v[idx];
}   inline operator double () { return sum; } inline double  & getRet() { return sum; }  inline const Grid<Real>& getArg0() { return v; } typedef Grid<Real> type0; void run() {  const int _sz = size; for (int i=0; i < _sz; i++) op(i, v,sum);  } const Grid<Real>& v;  double sum;  };



 struct minReduction : public KernelBase { minReduction(const Grid<Real>& v) :  KernelBase(&v,0) ,v(v) ,sum(0)  { run(); }  inline void op(int idx, const Grid<Real>& v ,double& sum)  {
	if (sum < v[idx])
		sum = v[idx];
}   inline operator double () { return sum; } inline double  & getRet() { return sum; }  inline const Grid<Real>& getArg0() { return v; } typedef Grid<Real> type0; void run() {  const int _sz = size; for (int i=0; i < _sz; i++) op(i, v,sum);  } const Grid<Real>& v;  double sum;  };


void getCurl(MACGrid& vel, Grid<Real>& vort, int comp) {
	Grid<Vec3> velCenter(vel.getParent()), curl(vel.getParent());
	
	GetCentered(velCenter, vel);
	CurlOp(velCenter, curl);
	GetComponent(curl, vort, comp);
} static PyObject* _W_0 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "getCurl" ); PyObject *_retval = 0; { ArgLocker _lock; MACGrid& vel = *_args.getPtr<MACGrid >("vel",0,&_lock); Grid<Real>& vort = *_args.getPtr<Grid<Real> >("vort",1,&_lock); int comp = _args.get<int >("comp",2,&_lock);   _retval = getPyNone(); getCurl(vel,vort,comp);  _args.check(); } pbFinalizePlugin(parent,"getCurl" ); return _retval; } catch(std::exception& e) { pbSetError("getCurl",e.what()); return 0; } } static const Pb::Register _RP_getCurl ("","getCurl",_W_0); 

void setinflow(FlagGrid& flags, MACGrid& vel, LevelsetGrid& phi, Real h) {
	FOR_IJK(vel) {
		if (i<=2) {
			if (j < h*flags.getSizeY()) {
				vel(i,j,k).x = 1;            
				if (!flags.isObstacle(i,j,k)) { 
					flags(i,j,k) = 1;        
					phi(i,j,k) = -1;
				}                
			} else {
				vel(i,j,k).x = 0;                            
				if (!flags.isObstacle(i,j,k)) { 
					flags(i,j,k) = 4;
					phi(i,j,k) = 1;
				}
			}
		}
		else if (i>=flags.getSizeX()-2) {
			vel(i,j,k).x = 1;            
		}
	}
} static PyObject* _W_1 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "setinflow" ); PyObject *_retval = 0; { ArgLocker _lock; FlagGrid& flags = *_args.getPtr<FlagGrid >("flags",0,&_lock); MACGrid& vel = *_args.getPtr<MACGrid >("vel",1,&_lock); LevelsetGrid& phi = *_args.getPtr<LevelsetGrid >("phi",2,&_lock); Real h = _args.get<Real >("h",3,&_lock);   _retval = getPyNone(); setinflow(flags,vel,phi,h);  _args.check(); } pbFinalizePlugin(parent,"setinflow" ); return _retval; } catch(std::exception& e) { pbSetError("setinflow",e.what()); return 0; } } static const Pb::Register _RP_setinflow ("","setinflow",_W_1); 
	
void testDiscardNth(BasicParticleSystem& parts, int skip=1) { 
	//knSetPdataConst<Real>(pd,value); 
	for(int i=0; i<parts.size(); ++i) {
		if(i%(skip+1) == skip) { // keep 
		} else {
			parts.setPos(i, Vec3(-100000) );
		}
	}
} static PyObject* _W_2 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "testDiscardNth" ); PyObject *_retval = 0; { ArgLocker _lock; BasicParticleSystem& parts = *_args.getPtr<BasicParticleSystem >("parts",0,&_lock); int skip = _args.getOpt<int >("skip",1,1,&_lock);   _retval = getPyNone(); testDiscardNth(parts,skip);  _args.check(); } pbFinalizePlugin(parent,"testDiscardNth" ); return _retval; } catch(std::exception& e) { pbSetError("testDiscardNth",e.what()); return 0; } } static const Pb::Register _RP_testDiscardNth ("","testDiscardNth",_W_2); 



} //namespace



