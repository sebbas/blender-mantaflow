




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
 * Tools to setup fields and inflows
 *
 ******************************************************************************/

#include "vectorbase.h"
#include "shapes.h"
#include "commonkernels.h"
#include "particle.h"
#include "noisefield.h"
#include "mesh.h"

using namespace std;

namespace Manta {
	
	//! Apply noise to grid
	
	
	struct KnApplyNoise : public KernelBase { KnApplyNoise(FlagGrid& flags, Grid<Real>& density, WaveletNoiseField& noise, Grid<Real>& sdf, Real scale, Real sigma) :  KernelBase(&flags,0) ,flags(flags),density(density),noise(noise),sdf(sdf),scale(scale),sigma(sigma)   { run(); }  inline void op(int i, int j, int k, FlagGrid& flags, Grid<Real>& density, WaveletNoiseField& noise, Grid<Real>& sdf, Real scale, Real sigma )  {
		if (!flags.isFluid(i,j,k) || sdf(i,j,k) > sigma) return;
		Real factor = clamp(1.0-0.5/sigma * (sdf(i,j,k)+sigma), 0.0, 1.0);
		
		Real target = noise.evaluate(Vec3(i,j,k)) * scale * factor;
		if (density(i,j,k) < target)
			density(i,j,k) = target;
	}   inline FlagGrid& getArg0() { return flags; } typedef FlagGrid type0;inline Grid<Real>& getArg1() { return density; } typedef Grid<Real> type1;inline WaveletNoiseField& getArg2() { return noise; } typedef WaveletNoiseField type2;inline Grid<Real>& getArg3() { return sdf; } typedef Grid<Real> type3;inline Real& getArg4() { return scale; } typedef Real type4;inline Real& getArg5() { return sigma; } typedef Real type5; void run() {  const int _maxX = maxX; const int _maxY = maxY; for (int k=minZ; k< maxZ; k++) for (int j=0; j< _maxY; j++) for (int i=0; i< _maxX; i++) op(i,j,k, flags,density,noise,sdf,scale,sigma);  } FlagGrid& flags; Grid<Real>& density; WaveletNoiseField& noise; Grid<Real>& sdf; Real scale; Real sigma;   };
	
	
	
	struct KnApplyDensity : public KernelBase { KnApplyDensity(FlagGrid& flags, Grid<Real>& density, Grid<Real>& sdf, Real value, Real sigma) :  KernelBase(&flags,0) ,flags(flags),density(density),sdf(sdf),value(value),sigma(sigma)   { run(); }  inline void op(int i, int j, int k, FlagGrid& flags, Grid<Real>& density, Grid<Real>& sdf, Real value, Real sigma )  {
		if (!flags.isFluid(i,j,k) || sdf(i,j,k) > sigma) return;
		density(i,j,k) = value;
	}   inline FlagGrid& getArg0() { return flags; } typedef FlagGrid type0;inline Grid<Real>& getArg1() { return density; } typedef Grid<Real> type1;inline Grid<Real>& getArg2() { return sdf; } typedef Grid<Real> type2;inline Real& getArg3() { return value; } typedef Real type3;inline Real& getArg4() { return sigma; } typedef Real type4; void run() {  const int _maxX = maxX; const int _maxY = maxY; for (int k=minZ; k< maxZ; k++) for (int j=0; j< _maxY; j++) for (int i=0; i< _maxX; i++) op(i,j,k, flags,density,sdf,value,sigma);  } FlagGrid& flags; Grid<Real>& density; Grid<Real>& sdf; Real value; Real sigma;   };
	//! Init noise-modulated density inside shape
	
	void densityInflow(FlagGrid& flags, Grid<Real>& density, WaveletNoiseField& noise, Shape* shape, Real scale=1.0, Real sigma=0) {
		Grid<Real> sdf = shape->computeLevelset();
		KnApplyNoise(flags, density, noise, sdf, scale, sigma);
	} static PyObject* _W_0 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "densityInflow" ); PyObject *_retval = 0; { ArgLocker _lock; FlagGrid& flags = *_args.getPtr<FlagGrid >("flags",0,&_lock); Grid<Real>& density = *_args.getPtr<Grid<Real> >("density",1,&_lock); WaveletNoiseField& noise = *_args.getPtr<WaveletNoiseField >("noise",2,&_lock); Shape* shape = _args.getPtr<Shape >("shape",3,&_lock); Real scale = _args.getOpt<Real >("scale",4,1.0,&_lock); Real sigma = _args.getOpt<Real >("sigma",5,0,&_lock);   _retval = getPyNone(); densityInflow(flags,density,noise,shape,scale,sigma);  _args.check(); } pbFinalizePlugin(parent,"densityInflow" ); return _retval; } catch(std::exception& e) { pbSetError("densityInflow",e.what()); return 0; } } static const Pb::Register _RP_densityInflow ("","densityInflow",_W_0); 
	
	
	//! Init noise-modulated density inside mesh
	
	void densityInflowMeshNoise(FlagGrid& flags, Grid<Real>& density, WaveletNoiseField& noise, Mesh* mesh, Real scale=1.0, Real sigma=0) {
		FluidSolver dummy(density.getSize());
		LevelsetGrid sdf(&dummy, false);
		mesh->meshSDF(*mesh, sdf, 1.);
		KnApplyNoise(flags, density, noise, sdf, scale, sigma);
	} static PyObject* _W_1 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "densityInflowMeshNoise" ); PyObject *_retval = 0; { ArgLocker _lock; FlagGrid& flags = *_args.getPtr<FlagGrid >("flags",0,&_lock); Grid<Real>& density = *_args.getPtr<Grid<Real> >("density",1,&_lock); WaveletNoiseField& noise = *_args.getPtr<WaveletNoiseField >("noise",2,&_lock); Mesh* mesh = _args.getPtr<Mesh >("mesh",3,&_lock); Real scale = _args.getOpt<Real >("scale",4,1.0,&_lock); Real sigma = _args.getOpt<Real >("sigma",5,0,&_lock);   _retval = getPyNone(); densityInflowMeshNoise(flags,density,noise,mesh,scale,sigma);  _args.check(); } pbFinalizePlugin(parent,"densityInflowMeshNoise" ); return _retval; } catch(std::exception& e) { pbSetError("densityInflowMeshNoise",e.what()); return 0; } } static const Pb::Register _RP_densityInflowMeshNoise ("","densityInflowMeshNoise",_W_1); 
	//! Init still density inside mesh
	
	void densityInflowMesh(FlagGrid& flags, Grid<Real>& density, Mesh* mesh, Real value=1., Real cutoff = 7, Real sigma=0) {
		FluidSolver dummy(density.getSize());
		LevelsetGrid sdf(&dummy, false);
		mesh->meshSDF(*mesh, sdf, 2.,cutoff);
		KnApplyDensity(flags, density, sdf, value, sigma);
	} static PyObject* _W_2 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "densityInflowMesh" ); PyObject *_retval = 0; { ArgLocker _lock; FlagGrid& flags = *_args.getPtr<FlagGrid >("flags",0,&_lock); Grid<Real>& density = *_args.getPtr<Grid<Real> >("density",1,&_lock); Mesh* mesh = _args.getPtr<Mesh >("mesh",2,&_lock); Real value = _args.getOpt<Real >("value",3,1.,&_lock); Real cutoff = _args.getOpt<Real >("cutoff",4,7,&_lock); Real sigma = _args.getOpt<Real >("sigma",5,0,&_lock);   _retval = getPyNone(); densityInflowMesh(flags,density,mesh,value,cutoff,sigma);  _args.check(); } pbFinalizePlugin(parent,"densityInflowMesh" ); return _retval; } catch(std::exception& e) { pbSetError("densityInflowMesh",e.what()); return 0; } } static const Pb::Register _RP_densityInflowMesh ("","densityInflowMesh",_W_2); 
	//! sample noise field and set pdata with its values (for convenience, scale the noise values)
	
	template <class T>  struct knSetPdataNoise : public KernelBase { knSetPdataNoise(BasicParticleSystem& parts, ParticleDataImpl<T>& pdata, WaveletNoiseField& noise, Real scale) :  KernelBase(parts.size()) ,parts(parts),pdata(pdata),noise(noise),scale(scale)   { run(); }  inline void op(int idx, BasicParticleSystem& parts, ParticleDataImpl<T>& pdata, WaveletNoiseField& noise, Real scale )  {
		pdata[idx] = noise.evaluate( parts.getPos(idx) ) * scale;
	}   inline BasicParticleSystem& getArg0() { return parts; } typedef BasicParticleSystem type0;inline ParticleDataImpl<T>& getArg1() { return pdata; } typedef ParticleDataImpl<T> type1;inline WaveletNoiseField& getArg2() { return noise; } typedef WaveletNoiseField type2;inline Real& getArg3() { return scale; } typedef Real type3; void run() {  const int _sz = size; for (int i=0; i < _sz; i++) op(i, parts,pdata,noise,scale);  } BasicParticleSystem& parts; ParticleDataImpl<T>& pdata; WaveletNoiseField& noise; Real scale;   };
	
	template <class T>  struct knSetPdataNoiseVec : public KernelBase { knSetPdataNoiseVec(BasicParticleSystem& parts, ParticleDataImpl<T>& pdata, WaveletNoiseField& noise, Real scale) :  KernelBase(parts.size()) ,parts(parts),pdata(pdata),noise(noise),scale(scale)   { run(); }  inline void op(int idx, BasicParticleSystem& parts, ParticleDataImpl<T>& pdata, WaveletNoiseField& noise, Real scale )  {
		pdata[idx] = noise.evaluateVec( parts.getPos(idx) ) * scale;
	}   inline BasicParticleSystem& getArg0() { return parts; } typedef BasicParticleSystem type0;inline ParticleDataImpl<T>& getArg1() { return pdata; } typedef ParticleDataImpl<T> type1;inline WaveletNoiseField& getArg2() { return noise; } typedef WaveletNoiseField type2;inline Real& getArg3() { return scale; } typedef Real type3; void run() {  const int _sz = size; for (int i=0; i < _sz; i++) op(i, parts,pdata,noise,scale);  } BasicParticleSystem& parts; ParticleDataImpl<T>& pdata; WaveletNoiseField& noise; Real scale;   };
	void setNoisePdata(BasicParticleSystem& parts, ParticleDataImpl<Real>& pd, WaveletNoiseField& noise, Real scale=1.) { knSetPdataNoise<Real>(parts, pd,noise,scale); } static PyObject* _W_3 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "setNoisePdata" ); PyObject *_retval = 0; { ArgLocker _lock; BasicParticleSystem& parts = *_args.getPtr<BasicParticleSystem >("parts",0,&_lock); ParticleDataImpl<Real>& pd = *_args.getPtr<ParticleDataImpl<Real> >("pd",1,&_lock); WaveletNoiseField& noise = *_args.getPtr<WaveletNoiseField >("noise",2,&_lock); Real scale = _args.getOpt<Real >("scale",3,1.,&_lock);   _retval = getPyNone(); setNoisePdata(parts,pd,noise,scale);  _args.check(); } pbFinalizePlugin(parent,"setNoisePdata" ); return _retval; } catch(std::exception& e) { pbSetError("setNoisePdata",e.what()); return 0; } } static const Pb::Register _RP_setNoisePdata ("","setNoisePdata",_W_3); 
	void setNoisePdataVec3(BasicParticleSystem& parts, ParticleDataImpl<Vec3>& pd, WaveletNoiseField& noise, Real scale=1.) { knSetPdataNoiseVec<Vec3>(parts, pd,noise,scale); } static PyObject* _W_4 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "setNoisePdataVec3" ); PyObject *_retval = 0; { ArgLocker _lock; BasicParticleSystem& parts = *_args.getPtr<BasicParticleSystem >("parts",0,&_lock); ParticleDataImpl<Vec3>& pd = *_args.getPtr<ParticleDataImpl<Vec3> >("pd",1,&_lock); WaveletNoiseField& noise = *_args.getPtr<WaveletNoiseField >("noise",2,&_lock); Real scale = _args.getOpt<Real >("scale",3,1.,&_lock);   _retval = getPyNone(); setNoisePdataVec3(parts,pd,noise,scale);  _args.check(); } pbFinalizePlugin(parent,"setNoisePdataVec3" ); return _retval; } catch(std::exception& e) { pbSetError("setNoisePdataVec3",e.what()); return 0; } } static const Pb::Register _RP_setNoisePdataVec3 ("","setNoisePdataVec3",_W_4); 
	void setNoisePdataInt(BasicParticleSystem& parts, ParticleDataImpl<int >& pd, WaveletNoiseField& noise, Real scale=1.) { knSetPdataNoise<int> (parts, pd,noise,scale); } static PyObject* _W_5 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "setNoisePdataInt" ); PyObject *_retval = 0; { ArgLocker _lock; BasicParticleSystem& parts = *_args.getPtr<BasicParticleSystem >("parts",0,&_lock); ParticleDataImpl<int >& pd = *_args.getPtr<ParticleDataImpl<int > >("pd",1,&_lock); WaveletNoiseField& noise = *_args.getPtr<WaveletNoiseField >("noise",2,&_lock); Real scale = _args.getOpt<Real >("scale",3,1.,&_lock);   _retval = getPyNone(); setNoisePdataInt(parts,pd,noise,scale);  _args.check(); } pbFinalizePlugin(parent,"setNoisePdataInt" ); return _retval; } catch(std::exception& e) { pbSetError("setNoisePdataInt",e.what()); return 0; } } static const Pb::Register _RP_setNoisePdataInt ("","setNoisePdataInt",_W_5); 
	
	//! SDF gradient from obstacle flags
	Grid<Vec3> obstacleGradient(FlagGrid& flags) {
		LevelsetGrid levelset(flags.getParent(),false);
		Grid<Vec3> gradient(flags.getParent());
		
		// rebuild obstacle levelset
		FOR_IDX(levelset) {
			levelset[idx] = flags.isObstacle(idx) ? -0.5 : 0.5;
		}
		levelset.reinitMarching(flags, 6.0, 0, true, false, FlagGrid::TypeReserved);
		
		// build levelset gradient
		GradientOp(gradient, levelset);
		
		FOR_IDX(levelset) {
			Vec3 grad = gradient[idx];
			Real s = normalize(grad);
			if (s <= 0.1 || levelset[idx] >= 0) 
				grad=Vec3(0.);        
			gradient[idx] = grad * levelset[idx];
		}
		
		return gradient;
	} static PyObject* _W_6 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "obstacleGradient" ); PyObject *_retval = 0; { ArgLocker _lock; FlagGrid& flags = *_args.getPtr<FlagGrid >("flags",0,&_lock);   _retval = toPy(obstacleGradient(flags));  _args.check(); } pbFinalizePlugin(parent,"obstacleGradient" ); return _retval; } catch(std::exception& e) { pbSetError("obstacleGradient",e.what()); return 0; } } static const Pb::Register _RP_obstacleGradient ("","obstacleGradient",_W_6); 
	
	LevelsetGrid obstacleLevelset(FlagGrid& flags) {
		LevelsetGrid levelset(flags.getParent(),false);
		Grid<Vec3> gradient(flags.getParent());
		
		// rebuild obstacle levelset
		FOR_IDX(levelset) {
			levelset[idx] = flags.isObstacle(idx) ? -0.5 : 0.5;
		}
		levelset.reinitMarching(flags, 6.0, 0, true, false, FlagGrid::TypeReserved);
		
		return levelset;
	} static PyObject* _W_7 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "obstacleLevelset" ); PyObject *_retval = 0; { ArgLocker _lock; FlagGrid& flags = *_args.getPtr<FlagGrid >("flags",0,&_lock);   _retval = toPy(obstacleLevelset(flags));  _args.check(); } pbFinalizePlugin(parent,"obstacleLevelset" ); return _retval; } catch(std::exception& e) { pbSetError("obstacleLevelset",e.what()); return 0; } } static const Pb::Register _RP_obstacleLevelset ("","obstacleLevelset",_W_7);     
	
	
} // namespace


