




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
 * Use this file to test new functionality
 *
 ******************************************************************************/

#include "levelset.h"
#include "commonkernels.h"
#include "particle.h"
#include <cmath>

using namespace std;

namespace Manta {



template <class S> void addToGrid(Grid<S>& a, S v) {
	FOR_IDX(a) a[idx] += v;
}template <class S>  static PyObject* _W_T_0 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "addToGrid" ); PyObject *_retval = 0; { ArgLocker _lock; Grid<S>& a = *_args.getPtr<Grid<S> >("a",0,&_lock); S v = _args.get<S >("v",1,&_lock);   _retval = getPyNone(); addToGrid(a,v);  _args.check(); } pbFinalizePlugin(parent,"addToGrid" ); return _retval; } catch(std::exception& e) { pbSetError("addToGrid",e.what()); return 0; } }template <class S>  static bool _K_0 (PbArgs& A) { return A.typeCheck<Grid<S> >(0,"a") && A.typeCheck<S >(1,"v"); }static PyObject* _W_0 (PyObject* s, PyObject* l, PyObject* kw) { PbArgs args(l, kw); int hits=0; PyObject* (*call)(PyObject*,PyObject*,PyObject*); if (_K_0<int>(args)) {hits++; call = _W_T_0<int>; }if (_K_0<Real>(args)) {hits++; call = _W_T_0<Real>; }if (_K_0<Vec3>(args)) {hits++; call = _W_T_0<Vec3>; } if (hits == 1) return call(s,l,kw); if (hits == 0) pbSetError("addToGrid", "Can't deduce template parameters"); else pbSetError("addToGrid", "Argument matches multiple templates"); return  0  ; } static const Pb::Register _RP_addToGrid ("","addToGrid",_W_0); 


//! Kernel: get component (not shifted)
/*KERNEL(idx) returns(Grid<Real> ret(parent))
Grid<Real> GetComponent2(const Grid<Vec3>& grid, int dim) {
	ret[idx] = grid[idx][dim];
};

PYTHON void testp(Grid<Vec3>& b) {
	Grid<Real> d(b.getParent());
	b(20,20,20) = Vec3(21,22,23); 
	{
		cout <<"middle" << endl;        
		Grid<Real> a = GetComponent2(b,0);
		cout << a(20,20,20) << endl;        
		cout <<"middle" << endl;        
	}
	cout << "end" << endl;errMsg("f");
}
*/



 struct ddtest : public KernelBase { ddtest(const Grid<Real>& v) :  KernelBase(&v,0) ,v(v) ,sum(0)  { run(); }  inline void op(int idx, const Grid<Real>& v ,double& sum)  {
	sum += v[idx];
}   inline operator double () { return sum; } inline double  & getRet() { return sum; }  inline const Grid<Real>& getArg0() { return v; } typedef Grid<Real> type0; void run() {  const int _sz = size; for (int i=0; i < _sz; i++) op(i, v,sum);  } const Grid<Real>& v;  double sum;  };



 struct detest : public KernelBase { detest(const Grid<Real>& v) :  KernelBase(&v,0) ,v(v) ,sum(0)  { run(); }  inline void op(int idx, const Grid<Real>& v ,double& sum)  {
	if (sum < v[idx])
		sum = v[idx];
}   inline operator double () { return sum; } inline double  & getRet() { return sum; }  inline const Grid<Real>& getArg0() { return v; } typedef Grid<Real> type0; void run() {  const int _sz = size; for (int i=0; i < _sz; i++) op(i, v,sum);  } const Grid<Real>& v;  double sum;  };

void checkGrids(Grid<int>& flags1, Grid<int>& flags2, Grid<Real>& phi1, Grid<Real>& phi2, Grid<Vec3>& vel1, Grid<Vec3>& vel2) {
	FOR_IJK(flags1) {
		assertMsg(flags1(i,j,k) == flags2(i,j,k), "flags mismatch");
		assertMsg(norm(vel1(i,j,k)-vel2(i,j,k)) < 1e-1, "vel mismatch");
		assertMsg( fabs(phi1(i,j,k)-phi2(i,j,k)) < 1e-4, "phi mismatch");
	}
} static PyObject* _W_1 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "checkGrids" ); PyObject *_retval = 0; { ArgLocker _lock; Grid<int>& flags1 = *_args.getPtr<Grid<int> >("flags1",0,&_lock); Grid<int>& flags2 = *_args.getPtr<Grid<int> >("flags2",1,&_lock); Grid<Real>& phi1 = *_args.getPtr<Grid<Real> >("phi1",2,&_lock); Grid<Real>& phi2 = *_args.getPtr<Grid<Real> >("phi2",3,&_lock); Grid<Vec3>& vel1 = *_args.getPtr<Grid<Vec3> >("vel1",4,&_lock); Grid<Vec3>& vel2 = *_args.getPtr<Grid<Vec3> >("vel2",5,&_lock);   _retval = getPyNone(); checkGrids(flags1,flags2,phi1,phi2,vel1,vel2);  _args.check(); } pbFinalizePlugin(parent,"checkGrids" ); return _retval; } catch(std::exception& e) { pbSetError("checkGrids",e.what()); return 0; } } static const Pb::Register _RP_checkGrids ("","checkGrids",_W_1); 


struct myvec {
	myvec(int n) : x(n) { cout << "constructor" << endl; };
	myvec(const myvec& a) : x(a.x) { cout << "copy constructor" << endl; }
	myvec& operator=(const myvec& a) { x=a.x; cout << "copy operator" << endl; return *this;}
	int& operator[](int idx) { return x[idx]; }
	
	vector<int> x;
};


 struct testy : public KernelBase { testy(vector<int>& a) :  KernelBase(a.size()) ,a(a) ,vec((size))  { run(); }  inline void op(int idx, vector<int>& a ,myvec& vec)  {
	vec[idx] = a[idx];
}   inline operator myvec () { return vec; } inline myvec  & getRet() { return vec; }  inline vector<int>& getArg0() { return a; } typedef vector<int> type0; void run() {  const int _sz = size; for (int i=0; i < _sz; i++) op(i, a,vec);  } vector<int>& a;  myvec vec;  };

void kernelTest() {
	cout << "kernel test" << endl;
	vector<int> a(10);
	for (int i=0;i<10;i++) a[i]=i;
	
	//testy xx(a);
	myvec b = testy(a);
	for (int i=0;i<10;i++) cout << b[i] << endl;
	cout << "kernel end" << endl;
} static PyObject* _W_2 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "kernelTest" ); PyObject *_retval = 0; { ArgLocker _lock;   _retval = getPyNone(); kernelTest();  _args.check(); } pbFinalizePlugin(parent,"kernelTest" ); return _retval; } catch(std::exception& e) { pbSetError("kernelTest",e.what()); return 0; } } static const Pb::Register _RP_kernelTest ("","kernelTest",_W_2); 

void getCurl(MACGrid& vel, Grid<Real>& vort, int comp) {
	Grid<Vec3> velCenter(vel.getParent()), curl(vel.getParent());
	
	GetCentered(velCenter, vel);
	CurlOp(velCenter, curl);
	GetComponent(curl, vort, comp);
} static PyObject* _W_3 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "getCurl" ); PyObject *_retval = 0; { ArgLocker _lock; MACGrid& vel = *_args.getPtr<MACGrid >("vel",0,&_lock); Grid<Real>& vort = *_args.getPtr<Grid<Real> >("vort",1,&_lock); int comp = _args.get<int >("comp",2,&_lock);   _retval = getPyNone(); getCurl(vel,vort,comp);  _args.check(); } pbFinalizePlugin(parent,"getCurl" ); return _retval; } catch(std::exception& e) { pbSetError("getCurl",e.what()); return 0; } } static const Pb::Register _RP_getCurl ("","getCurl",_W_3); 

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
			/*if (j < 30-12) {
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
			}*/
		}
	}
} static PyObject* _W_4 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "setinflow" ); PyObject *_retval = 0; { ArgLocker _lock; FlagGrid& flags = *_args.getPtr<FlagGrid >("flags",0,&_lock); MACGrid& vel = *_args.getPtr<MACGrid >("vel",1,&_lock); LevelsetGrid& phi = *_args.getPtr<LevelsetGrid >("phi",2,&_lock); Real h = _args.get<Real >("h",3,&_lock);   _retval = getPyNone(); setinflow(flags,vel,phi,h);  _args.check(); } pbFinalizePlugin(parent,"setinflow" ); return _retval; } catch(std::exception& e) { pbSetError("setinflow",e.what()); return 0; } } static const Pb::Register _RP_setinflow ("","setinflow",_W_4); 
	
void testDiscardNth(BasicParticleSystem& parts, int skip=1) { 
	//knSetPdataConst<Real>(pd,value); 
	for(int i=0; i<parts.size(); ++i) {
		if(i%(skip+1) == skip) { // keep 
		} else {
			parts.setPos(i, Vec3(-100000) );
		}
	}
} static PyObject* _W_5 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "testDiscardNth" ); PyObject *_retval = 0; { ArgLocker _lock; BasicParticleSystem& parts = *_args.getPtr<BasicParticleSystem >("parts",0,&_lock); int skip = _args.getOpt<int >("skip",1,1,&_lock);   _retval = getPyNone(); testDiscardNth(parts,skip);  _args.check(); } pbFinalizePlugin(parent,"testDiscardNth" ); return _retval; } catch(std::exception& e) { pbSetError("testDiscardNth",e.what()); return 0; } } static const Pb::Register _RP_testDiscardNth ("","testDiscardNth",_W_5); 

} //namespace



