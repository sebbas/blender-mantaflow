




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/Users/sbarschkis/Developer/Mantaflow/blenderIntegration/mantaflowgit/source/plugin/initplugins.cpp"
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
#include "simpleimage.h"
#include "mesh.h"

using namespace std;

namespace Manta {
	
//! Apply noise to grid


 struct KnApplyNoiseInfl : public KernelBase { KnApplyNoiseInfl(FlagGrid& flags, Grid<Real>& density, WaveletNoiseField& noise, Grid<Real>& sdf, Real scale, Real sigma) :  KernelBase(&flags,0) ,flags(flags),density(density),noise(noise),sdf(sdf),scale(scale),sigma(sigma)   { runMessage(); run(); }  inline void op(int i, int j, int k, FlagGrid& flags, Grid<Real>& density, WaveletNoiseField& noise, Grid<Real>& sdf, Real scale, Real sigma )  {
	if (!flags.isFluid(i,j,k) || sdf(i,j,k) > sigma) return;
	Real factor = clamp(1.0-0.5/sigma * (sdf(i,j,k)+sigma), 0.0, 1.0);
	
	Real target = noise.evaluate(Vec3(i,j,k)) * scale * factor;
	if (density(i,j,k) < target)
		density(i,j,k) = target;
}   inline FlagGrid& getArg0() { return flags; } typedef FlagGrid type0;inline Grid<Real>& getArg1() { return density; } typedef Grid<Real> type1;inline WaveletNoiseField& getArg2() { return noise; } typedef WaveletNoiseField type2;inline Grid<Real>& getArg3() { return sdf; } typedef Grid<Real> type3;inline Real& getArg4() { return scale; } typedef Real type4;inline Real& getArg5() { return sigma; } typedef Real type5; void runMessage() { debMsg("Executing kernel KnApplyNoiseInfl ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int k=minZ; k < maxZ; k++) for (int j=0; j < _maxY; j++) for (int i=0; i < _maxX; i++) op(i,j,k,flags,density,noise,sdf,scale,sigma);  } } else { const int k=0; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int j=0; j < _maxY; j++) for (int i=0; i < _maxX; i++) op(i,j,k,flags,density,noise,sdf,scale,sigma);  } }  } FlagGrid& flags; Grid<Real>& density; WaveletNoiseField& noise; Grid<Real>& sdf; Real scale; Real sigma;   };
#line 29 "plugin/initplugins.cpp"



//! Init noise-modulated density inside shape

void densityInflow(FlagGrid& flags, Grid<Real>& density, WaveletNoiseField& noise, Shape* shape, Real scale=1.0, Real sigma=0) {
	Grid<Real> sdf = shape->computeLevelset();
	KnApplyNoiseInfl(flags, density, noise, sdf, scale, sigma);
} static PyObject* _W_0 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "densityInflow" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; FlagGrid& flags = *_args.getPtr<FlagGrid >("flags",0,&_lock); Grid<Real>& density = *_args.getPtr<Grid<Real> >("density",1,&_lock); WaveletNoiseField& noise = *_args.getPtr<WaveletNoiseField >("noise",2,&_lock); Shape* shape = _args.getPtr<Shape >("shape",3,&_lock); Real scale = _args.getOpt<Real >("scale",4,1.0,&_lock); Real sigma = _args.getOpt<Real >("sigma",5,0,&_lock);   _retval = getPyNone(); densityInflow(flags,density,noise,shape,scale,sigma);  _args.check(); } pbFinalizePlugin(parent,"densityInflow", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("densityInflow",e.what()); return 0; } } static const Pb::Register _RP_densityInflow ("","densityInflow",_W_0); 
//! Apply noise to real grid based on an SDF
 struct KnAddNoise : public KernelBase { KnAddNoise(FlagGrid& flags, Grid<Real>& density, WaveletNoiseField& noise, Grid<Real>* sdf, Real scale) :  KernelBase(&flags,0) ,flags(flags),density(density),noise(noise),sdf(sdf),scale(scale)   { runMessage(); run(); }  inline void op(int i, int j, int k, FlagGrid& flags, Grid<Real>& density, WaveletNoiseField& noise, Grid<Real>* sdf, Real scale )  {
	if (!flags.isFluid(i,j,k) || (sdf && (*sdf)(i,j,k) > 0.) ) return;
	density(i,j,k) += noise.evaluate(Vec3(i,j,k)) * scale;
}   inline FlagGrid& getArg0() { return flags; } typedef FlagGrid type0;inline Grid<Real>& getArg1() { return density; } typedef Grid<Real> type1;inline WaveletNoiseField& getArg2() { return noise; } typedef WaveletNoiseField type2;inline Grid<Real>* getArg3() { return sdf; } typedef Grid<Real> type3;inline Real& getArg4() { return scale; } typedef Real type4; void runMessage() { debMsg("Executing kernel KnAddNoise ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int k=minZ; k < maxZ; k++) for (int j=0; j < _maxY; j++) for (int i=0; i < _maxX; i++) op(i,j,k,flags,density,noise,sdf,scale);  } } else { const int k=0; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int j=0; j < _maxY; j++) for (int i=0; i < _maxX; i++) op(i,j,k,flags,density,noise,sdf,scale);  } }  } FlagGrid& flags; Grid<Real>& density; WaveletNoiseField& noise; Grid<Real>* sdf; Real scale;   };
#line 45 "plugin/initplugins.cpp"


void addNoise(FlagGrid& flags, Grid<Real>& density, WaveletNoiseField& noise, Grid<Real>* sdf=NULL, Real scale=1.0 ) {
	KnAddNoise(flags, density, noise, sdf, scale );
} static PyObject* _W_1 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "addNoise" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; FlagGrid& flags = *_args.getPtr<FlagGrid >("flags",0,&_lock); Grid<Real>& density = *_args.getPtr<Grid<Real> >("density",1,&_lock); WaveletNoiseField& noise = *_args.getPtr<WaveletNoiseField >("noise",2,&_lock); Grid<Real>* sdf = _args.getPtrOpt<Grid<Real> >("sdf",3,NULL,&_lock); Real scale = _args.getOpt<Real >("scale",4,1.0 ,&_lock);   _retval = getPyNone(); addNoise(flags,density,noise,sdf,scale);  _args.check(); } pbFinalizePlugin(parent,"addNoise", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("addNoise",e.what()); return 0; } } static const Pb::Register _RP_addNoise ("","addNoise",_W_1); 

//! sample noise field and set pdata with its values (for convenience, scale the noise values)

template <class T>  struct knSetPdataNoise : public KernelBase { knSetPdataNoise(BasicParticleSystem& parts, ParticleDataImpl<T>& pdata, WaveletNoiseField& noise, Real scale) :  KernelBase(parts.size()) ,parts(parts),pdata(pdata),noise(noise),scale(scale)   { runMessage(); run(); }   inline void op(IndexInt idx, BasicParticleSystem& parts, ParticleDataImpl<T>& pdata, WaveletNoiseField& noise, Real scale )  {
	pdata[idx] = noise.evaluate( parts.getPos(idx) ) * scale;
}    inline BasicParticleSystem& getArg0() { return parts; } typedef BasicParticleSystem type0;inline ParticleDataImpl<T>& getArg1() { return pdata; } typedef ParticleDataImpl<T> type1;inline WaveletNoiseField& getArg2() { return noise; } typedef WaveletNoiseField type2;inline Real& getArg3() { return scale; } typedef Real type3; void runMessage() { debMsg("Executing kernel knSetPdataNoise ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (IndexInt i = 0; i < _sz; i++) op(i,parts,pdata,noise,scale);  }   } BasicParticleSystem& parts; ParticleDataImpl<T>& pdata; WaveletNoiseField& noise; Real scale;   };
#line 55 "plugin/initplugins.cpp"



template <class T>  struct knSetPdataNoiseVec : public KernelBase { knSetPdataNoiseVec(BasicParticleSystem& parts, ParticleDataImpl<T>& pdata, WaveletNoiseField& noise, Real scale) :  KernelBase(parts.size()) ,parts(parts),pdata(pdata),noise(noise),scale(scale)   { runMessage(); run(); }   inline void op(IndexInt idx, BasicParticleSystem& parts, ParticleDataImpl<T>& pdata, WaveletNoiseField& noise, Real scale )  {
	pdata[idx] = noise.evaluateVec( parts.getPos(idx) ) * scale;
}    inline BasicParticleSystem& getArg0() { return parts; } typedef BasicParticleSystem type0;inline ParticleDataImpl<T>& getArg1() { return pdata; } typedef ParticleDataImpl<T> type1;inline WaveletNoiseField& getArg2() { return noise; } typedef WaveletNoiseField type2;inline Real& getArg3() { return scale; } typedef Real type3; void runMessage() { debMsg("Executing kernel knSetPdataNoiseVec ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (IndexInt i = 0; i < _sz; i++) op(i,parts,pdata,noise,scale);  }   } BasicParticleSystem& parts; ParticleDataImpl<T>& pdata; WaveletNoiseField& noise; Real scale;   };
#line 59 "plugin/initplugins.cpp"


void setNoisePdata(BasicParticleSystem& parts, ParticleDataImpl<Real>& pd, WaveletNoiseField& noise, Real scale=1.) { knSetPdataNoise<Real>(parts, pd,noise,scale); } static PyObject* _W_2 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "setNoisePdata" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; BasicParticleSystem& parts = *_args.getPtr<BasicParticleSystem >("parts",0,&_lock); ParticleDataImpl<Real>& pd = *_args.getPtr<ParticleDataImpl<Real> >("pd",1,&_lock); WaveletNoiseField& noise = *_args.getPtr<WaveletNoiseField >("noise",2,&_lock); Real scale = _args.getOpt<Real >("scale",3,1.,&_lock);   _retval = getPyNone(); setNoisePdata(parts,pd,noise,scale);  _args.check(); } pbFinalizePlugin(parent,"setNoisePdata", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("setNoisePdata",e.what()); return 0; } } static const Pb::Register _RP_setNoisePdata ("","setNoisePdata",_W_2); 
void setNoisePdataVec3(BasicParticleSystem& parts, ParticleDataImpl<Vec3>& pd, WaveletNoiseField& noise, Real scale=1.) { knSetPdataNoiseVec<Vec3>(parts, pd,noise,scale); } static PyObject* _W_3 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "setNoisePdataVec3" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; BasicParticleSystem& parts = *_args.getPtr<BasicParticleSystem >("parts",0,&_lock); ParticleDataImpl<Vec3>& pd = *_args.getPtr<ParticleDataImpl<Vec3> >("pd",1,&_lock); WaveletNoiseField& noise = *_args.getPtr<WaveletNoiseField >("noise",2,&_lock); Real scale = _args.getOpt<Real >("scale",3,1.,&_lock);   _retval = getPyNone(); setNoisePdataVec3(parts,pd,noise,scale);  _args.check(); } pbFinalizePlugin(parent,"setNoisePdataVec3", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("setNoisePdataVec3",e.what()); return 0; } } static const Pb::Register _RP_setNoisePdataVec3 ("","setNoisePdataVec3",_W_3); 
void setNoisePdataInt(BasicParticleSystem& parts, ParticleDataImpl<int >& pd, WaveletNoiseField& noise, Real scale=1.) { knSetPdataNoise<int> (parts, pd,noise,scale); } static PyObject* _W_4 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "setNoisePdataInt" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; BasicParticleSystem& parts = *_args.getPtr<BasicParticleSystem >("parts",0,&_lock); ParticleDataImpl<int >& pd = *_args.getPtr<ParticleDataImpl<int > >("pd",1,&_lock); WaveletNoiseField& noise = *_args.getPtr<WaveletNoiseField >("noise",2,&_lock); Real scale = _args.getOpt<Real >("scale",3,1.,&_lock);   _retval = getPyNone(); setNoisePdataInt(parts,pd,noise,scale);  _args.check(); } pbFinalizePlugin(parent,"setNoisePdataInt", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("setNoisePdataInt",e.what()); return 0; } } static const Pb::Register _RP_setNoisePdataInt ("","setNoisePdataInt",_W_4); 

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
} static PyObject* _W_5 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "obstacleGradient" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; FlagGrid& flags = *_args.getPtr<FlagGrid >("flags",0,&_lock);   _retval = toPy(obstacleGradient(flags));  _args.check(); } pbFinalizePlugin(parent,"obstacleGradient", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("obstacleGradient",e.what()); return 0; } } static const Pb::Register _RP_obstacleGradient ("","obstacleGradient",_W_5); 

LevelsetGrid obstacleLevelset(FlagGrid& flags) {
   LevelsetGrid levelset(flags.getParent(),false);
	Grid<Vec3> gradient(flags.getParent());

	// rebuild obstacle levelset
	FOR_IDX(levelset) {
		levelset[idx] = flags.isObstacle(idx) ? -0.5 : 0.5;
	}
	levelset.reinitMarching(flags, 6.0, 0, true, false, FlagGrid::TypeReserved);

	return levelset;
} static PyObject* _W_6 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "obstacleLevelset" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; FlagGrid& flags = *_args.getPtr<FlagGrid >("flags",0,&_lock);   _retval = toPy(obstacleLevelset(flags));  _args.check(); } pbFinalizePlugin(parent,"obstacleLevelset", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("obstacleLevelset",e.what()); return 0; } } static const Pb::Register _RP_obstacleLevelset ("","obstacleLevelset",_W_6);     


//*****************************************************************************
// blender init functions 



 struct KnApplyEmission : public KernelBase { KnApplyEmission(FlagGrid& flags, Grid<Real>& density, Grid<Real>& emission, bool isAbsolute) :  KernelBase(&flags,0) ,flags(flags),density(density),emission(emission),isAbsolute(isAbsolute)   { runMessage(); run(); }  inline void op(int i, int j, int k, FlagGrid& flags, Grid<Real>& density, Grid<Real>& emission, bool isAbsolute )  {
	if (!flags.isFluid(i,j,k) || emission(i,j,k) == 0.) return;
	if (isAbsolute)
		density(i,j,k) = emission(i,j,k);
	else
		density(i,j,k) += emission(i,j,k);
}   inline FlagGrid& getArg0() { return flags; } typedef FlagGrid type0;inline Grid<Real>& getArg1() { return density; } typedef Grid<Real> type1;inline Grid<Real>& getArg2() { return emission; } typedef Grid<Real> type2;inline bool& getArg3() { return isAbsolute; } typedef bool type3; void runMessage() { debMsg("Executing kernel KnApplyEmission ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int k=minZ; k < maxZ; k++) for (int j=0; j < _maxY; j++) for (int i=0; i < _maxX; i++) op(i,j,k,flags,density,emission,isAbsolute);  } } else { const int k=0; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int j=0; j < _maxY; j++) for (int i=0; i < _maxX; i++) op(i,j,k,flags,density,emission,isAbsolute);  } }  } FlagGrid& flags; Grid<Real>& density; Grid<Real>& emission; bool isAbsolute;   };
#line 110 "plugin/initplugins.cpp"



//! Add emission values
//isAbsolute: whether to add emission values to existing, or replace
void applyEmission(FlagGrid& flags, Grid<Real>& density, Grid<Real>& emission, bool isAbsolute) {
	KnApplyEmission(flags, density, emission, isAbsolute);
} static PyObject* _W_7 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "applyEmission" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; FlagGrid& flags = *_args.getPtr<FlagGrid >("flags",0,&_lock); Grid<Real>& density = *_args.getPtr<Grid<Real> >("density",1,&_lock); Grid<Real>& emission = *_args.getPtr<Grid<Real> >("emission",2,&_lock); bool isAbsolute = _args.get<bool >("isAbsolute",3,&_lock);   _retval = getPyNone(); applyEmission(flags,density,emission,isAbsolute);  _args.check(); } pbFinalizePlugin(parent,"applyEmission", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("applyEmission",e.what()); return 0; } } static const Pb::Register _RP_applyEmission ("","applyEmission",_W_7); 

// blender init functions for meshes



 struct KnApplyDensity : public KernelBase { KnApplyDensity(FlagGrid& flags, Grid<Real>& density, Grid<Real>& sdf, Real value, Real sigma) :  KernelBase(&flags,0) ,flags(flags),density(density),sdf(sdf),value(value),sigma(sigma)   { runMessage(); run(); }  inline void op(int i, int j, int k, FlagGrid& flags, Grid<Real>& density, Grid<Real>& sdf, Real value, Real sigma )  {
	if (!flags.isFluid(i,j,k) || sdf(i,j,k) > sigma) return;
	density(i,j,k) = value;
}   inline FlagGrid& getArg0() { return flags; } typedef FlagGrid type0;inline Grid<Real>& getArg1() { return density; } typedef Grid<Real> type1;inline Grid<Real>& getArg2() { return sdf; } typedef Grid<Real> type2;inline Real& getArg3() { return value; } typedef Real type3;inline Real& getArg4() { return sigma; } typedef Real type4; void runMessage() { debMsg("Executing kernel KnApplyDensity ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int k=minZ; k < maxZ; k++) for (int j=0; j < _maxY; j++) for (int i=0; i < _maxX; i++) op(i,j,k,flags,density,sdf,value,sigma);  } } else { const int k=0; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int j=0; j < _maxY; j++) for (int i=0; i < _maxX; i++) op(i,j,k,flags,density,sdf,value,sigma);  } }  } FlagGrid& flags; Grid<Real>& density; Grid<Real>& sdf; Real value; Real sigma;   };
#line 128 "plugin/initplugins.cpp"


//! Init noise-modulated density inside mesh

void densityInflowMeshNoise(FlagGrid& flags, Grid<Real>& density, WaveletNoiseField& noise, Mesh* mesh, Real scale=1.0, Real sigma=0) {
	LevelsetGrid sdf(density.getParent(), false);
	mesh->computeLevelset(sdf, 1.);
	KnApplyNoiseInfl(flags, density, noise, sdf, scale, sigma);
} static PyObject* _W_8 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "densityInflowMeshNoise" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; FlagGrid& flags = *_args.getPtr<FlagGrid >("flags",0,&_lock); Grid<Real>& density = *_args.getPtr<Grid<Real> >("density",1,&_lock); WaveletNoiseField& noise = *_args.getPtr<WaveletNoiseField >("noise",2,&_lock); Mesh* mesh = _args.getPtr<Mesh >("mesh",3,&_lock); Real scale = _args.getOpt<Real >("scale",4,1.0,&_lock); Real sigma = _args.getOpt<Real >("sigma",5,0,&_lock);   _retval = getPyNone(); densityInflowMeshNoise(flags,density,noise,mesh,scale,sigma);  _args.check(); } pbFinalizePlugin(parent,"densityInflowMeshNoise", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("densityInflowMeshNoise",e.what()); return 0; } } static const Pb::Register _RP_densityInflowMeshNoise ("","densityInflowMeshNoise",_W_8); 

//! Init constant density inside mesh

void densityInflowMesh(FlagGrid& flags, Grid<Real>& density, Mesh* mesh, Real value=1., Real cutoff = 7, Real sigma=0) {
	LevelsetGrid sdf(density.getParent(), false);
	mesh->computeLevelset(sdf, 2., cutoff);
	KnApplyDensity(flags, density, sdf, value, sigma);
} static PyObject* _W_9 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "densityInflowMesh" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; FlagGrid& flags = *_args.getPtr<FlagGrid >("flags",0,&_lock); Grid<Real>& density = *_args.getPtr<Grid<Real> >("density",1,&_lock); Mesh* mesh = _args.getPtr<Mesh >("mesh",2,&_lock); Real value = _args.getOpt<Real >("value",3,1.,&_lock); Real cutoff = _args.getOpt<Real >("cutoff",4,7,&_lock); Real sigma = _args.getOpt<Real >("sigma",5,0,&_lock);   _retval = getPyNone(); densityInflowMesh(flags,density,mesh,value,cutoff,sigma);  _args.check(); } pbFinalizePlugin(parent,"densityInflowMesh", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("densityInflowMesh",e.what()); return 0; } } static const Pb::Register _RP_densityInflowMesh ("","densityInflowMesh",_W_9); 


//*****************************************************************************

//! check for symmetry , optionally enfore by copying

void checkSymmetry( Grid<Real>& a, Grid<Real>* err=NULL, bool symmetrize=false, int axis=0, int bound=0) {
	const int c  = axis; 
	const int s = a.getSize()[c];
	FOR_IJK(a) { 
		Vec3i idx(i,j,k), mdx(i,j,k);
		mdx[c] = s-1-idx[c];
		if( bound>0 && ((!a.isInBounds(idx,bound)) || (!a.isInBounds(mdx,bound))) ) continue;

		if(err) (*err)(idx) = fabs( (double)(a(idx) - a(mdx) ) ); 
		if(symmetrize && (idx[c]<s/2)) {
			a(idx) = a(mdx);
		}
	}
} static PyObject* _W_10 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "checkSymmetry" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; Grid<Real>& a = *_args.getPtr<Grid<Real> >("a",0,&_lock); Grid<Real>* err = _args.getPtrOpt<Grid<Real> >("err",1,NULL,&_lock); bool symmetrize = _args.getOpt<bool >("symmetrize",2,false,&_lock); int axis = _args.getOpt<int >("axis",3,0,&_lock); int bound = _args.getOpt<int >("bound",4,0,&_lock);   _retval = getPyNone(); checkSymmetry(a,err,symmetrize,axis,bound);  _args.check(); } pbFinalizePlugin(parent,"checkSymmetry", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("checkSymmetry",e.what()); return 0; } } static const Pb::Register _RP_checkSymmetry ("","checkSymmetry",_W_10); 
//! check for symmetry , mac grid version


void checkSymmetryVec3( Grid<Vec3>& a, Grid<Real>* err=NULL, bool symmetrize=false , int axis=0, int bound=0, int disable=0) {
	if(err) err->setConst(0.);

	// each dimension is measured separately for flexibility (could be combined)
	const int c  = axis;
	const int o1 = (c+1)%3;
	const int o2 = (c+2)%3;

	// x
	if(! (disable&1) ) {
		const int s = a.getSize()[c]+1; 
		FOR_IJK(a) { 
			Vec3i idx(i,j,k), mdx(i,j,k);
			mdx[c] = s-1-idx[c]; 
			if(mdx[c] >= a.getSize()[c]) continue; 
			if( bound>0 && ((!a.isInBounds(idx,bound)) || (!a.isInBounds(mdx,bound))) ) continue;

			// special case: center "line" of values , should be zero!
			if(mdx[c] == idx[c] ) {
				if(err) (*err)(idx) += fabs( (double)( a(idx)[c] ) ); 
				if(symmetrize) a(idx)[c] = 0.;
				continue; 
			}

			// note - the a(mdx) component needs to be inverted here!
			if(err) (*err)(idx) += fabs( (double)( a(idx)[c]- (a(mdx)[c]*-1.) ) ); 
			if(symmetrize && (idx[c]<s/2)) {
				a(idx)[c] = a(mdx)[c] * -1.;
			}
		}
	}

	// y
	if(! (disable&2) ) {
		const int s = a.getSize()[c];
		FOR_IJK(a) { 
			Vec3i idx(i,j,k), mdx(i,j,k);
			mdx[c] = s-1-idx[c]; 
			if( bound>0 && ((!a.isInBounds(idx,bound)) || (!a.isInBounds(mdx,bound))) ) continue;

			if(err) (*err)(idx) += fabs( (double)( a(idx)[o1]-a(mdx)[o1] ) ); 
			if(symmetrize && (idx[c]<s/2)) {
				a(idx)[o1] = a(mdx)[o1];
			}
		}
	} 

	// z
	if(! (disable&4) ) {
		const int s = a.getSize()[c];
		FOR_IJK(a) { 
			Vec3i idx(i,j,k), mdx(i,j,k);
			mdx[c] = s-1-idx[c]; 
			if( bound>0 && ((!a.isInBounds(idx,bound)) || (!a.isInBounds(mdx,bound))) ) continue;

			if(err) (*err)(idx) += fabs( (double)( a(idx)[o2]-a(mdx)[o2] ) ); 
			if(symmetrize && (idx[c]<s/2)) {
				a(idx)[o2] = a(mdx)[o2];
			}
		}
	} 

} static PyObject* _W_11 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "checkSymmetryVec3" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; Grid<Vec3>& a = *_args.getPtr<Grid<Vec3> >("a",0,&_lock); Grid<Real>* err = _args.getPtrOpt<Grid<Real> >("err",1,NULL,&_lock); bool symmetrize = _args.getOpt<bool >("symmetrize",2,false ,&_lock); int axis = _args.getOpt<int >("axis",3,0,&_lock); int bound = _args.getOpt<int >("bound",4,0,&_lock); int disable = _args.getOpt<int >("disable",5,0,&_lock);   _retval = getPyNone(); checkSymmetryVec3(a,err,symmetrize,axis,bound,disable);  _args.check(); } pbFinalizePlugin(parent,"checkSymmetryVec3", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("checkSymmetryVec3",e.what()); return 0; } } static const Pb::Register _RP_checkSymmetryVec3 ("","checkSymmetryVec3",_W_11); 


// from simpleimage.cpp
void projectImg( SimpleImage& img, Grid<Real>& val, int shadeMode=0, Real scale=1.);

//! output shaded (all 3 axes at once for 3D)
//! shading modes: 0 smoke, 1 surfaces

void projectPpmFull( Grid<Real>& val, string name, int shadeMode=0, Real scale=1.) {
	SimpleImage img;
	projectImg( img, val, shadeMode, scale );
	img.writePpm( name );
} static PyObject* _W_12 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "projectPpmFull" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; Grid<Real>& val = *_args.getPtr<Grid<Real> >("val",0,&_lock); string name = _args.get<string >("name",1,&_lock); int shadeMode = _args.getOpt<int >("shadeMode",2,0,&_lock); Real scale = _args.getOpt<Real >("scale",3,1.,&_lock);   _retval = getPyNone(); projectPpmFull(val,name,shadeMode,scale);  _args.check(); } pbFinalizePlugin(parent,"projectPpmFull", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("projectPpmFull",e.what()); return 0; } } static const Pb::Register _RP_projectPpmFull ("","projectPpmFull",_W_12); 

// helper functions for pdata operator tests

//! init some test particles at the origin

void addTestParts( BasicParticleSystem& parts, int num) {
	for(int i=0; i<num; ++i)
		parts.addBuffered( Vec3(0,0,0) );

	parts.doCompress();
	parts.insertBufferedParticles();
} static PyObject* _W_13 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "addTestParts" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; BasicParticleSystem& parts = *_args.getPtr<BasicParticleSystem >("parts",0,&_lock); int num = _args.get<int >("num",1,&_lock);   _retval = getPyNone(); addTestParts(parts,num);  _args.check(); } pbFinalizePlugin(parent,"addTestParts", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("addTestParts",e.what()); return 0; } } static const Pb::Register _RP_addTestParts ("","addTestParts",_W_13); 

// calculate the difference between two pdata fields (note - slow!, not parallelized)

Real pdataMaxDiff( ParticleDataBase* a, ParticleDataBase* b ) {    
	double maxVal = 0.;
	//debMsg(" PD "<< a->getType()<<"  as"<<a->getSizeSlow()<<"  bs"<<b->getSizeSlow() , 1);
	assertMsg(a->getType()     == b->getType()    , "pdataMaxDiff problem - different pdata types!");
	assertMsg(a->getSizeSlow() == b->getSizeSlow(), "pdataMaxDiff problem - different pdata sizes!");
	
	if (a->getType() & ParticleDataBase::TypeReal) 
	{
		ParticleDataImpl<Real>& av = *dynamic_cast<ParticleDataImpl<Real>*>(a);
		ParticleDataImpl<Real>& bv = *dynamic_cast<ParticleDataImpl<Real>*>(b);
		FOR_PARTS(av) {
			maxVal = std::max(maxVal, (double)fabs( av[idx]-bv[idx] ));
		}
	} else if (a->getType() & ParticleDataBase::TypeInt) 
	{
		ParticleDataImpl<int>& av = *dynamic_cast<ParticleDataImpl<int>*>(a);
		ParticleDataImpl<int>& bv = *dynamic_cast<ParticleDataImpl<int>*>(b);
		FOR_PARTS(av) {
			maxVal = std::max(maxVal, (double)fabs( (double)av[idx]-bv[idx] ));
		}
	} else if (a->getType() & ParticleDataBase::TypeVec3) {
		ParticleDataImpl<Vec3>& av = *dynamic_cast<ParticleDataImpl<Vec3>*>(a);
		ParticleDataImpl<Vec3>& bv = *dynamic_cast<ParticleDataImpl<Vec3>*>(b);
		FOR_PARTS(av) {
			double d = 0.;
			for(int c=0; c<3; ++c) { 
				d += fabs( (double)av[idx][c] - (double)bv[idx][c] );
			}
			maxVal = std::max(maxVal, d );
		}
	} else {
		errMsg("pdataMaxDiff: Grid Type is not supported (only Real, Vec3, int)");    
	}

	return maxVal;
} static PyObject* _W_14 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "pdataMaxDiff" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; ParticleDataBase* a = _args.getPtr<ParticleDataBase >("a",0,&_lock); ParticleDataBase* b = _args.getPtr<ParticleDataBase >("b",1,&_lock);   _retval = toPy(pdataMaxDiff(a,b));  _args.check(); } pbFinalizePlugin(parent,"pdataMaxDiff", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("pdataMaxDiff",e.what()); return 0; } } static const Pb::Register _RP_pdataMaxDiff ("","pdataMaxDiff",_W_14); 

//*****************************************************************************
// helper functions for volume fractions (which are needed for second order obstacle boundaries)


 struct kninitVortexVelocity : public KernelBase { kninitVortexVelocity(Grid<Real> &phiObs, MACGrid& vel, const Vec3 &center, const Real &radius) :  KernelBase(&phiObs,0) ,phiObs(phiObs),vel(vel),center(center),radius(radius)   { runMessage(); run(); }  inline void op(int i, int j, int k, Grid<Real> &phiObs, MACGrid& vel, const Vec3 &center, const Real &radius )  {
	
	if(phiObs(i,j,k) >= -1.) {

		Real dx = i - center.x; if(dx>=0) dx -= .5; else dx += .5;
		Real dy = j - center.y;
		Real r = std::sqrt(dx*dx+dy*dy);
		Real alpha = atan2(dy,dx);

		vel(i,j,k).x = -std::sin(alpha)*(r/radius);

		dx = i - center.x;
		dy = j - center.y; if(dy>=0) dy -= .5; else dy += .5;
		r = std::sqrt(dx*dx+dy*dy);
		alpha = atan2(dy,dx);

		vel(i,j,k).y = std::cos(alpha)*(r/radius);

	}

}   inline Grid<Real> & getArg0() { return phiObs; } typedef Grid<Real>  type0;inline MACGrid& getArg1() { return vel; } typedef MACGrid type1;inline const Vec3& getArg2() { return center; } typedef Vec3 type2;inline const Real& getArg3() { return radius; } typedef Real type3; void runMessage() { debMsg("Executing kernel kninitVortexVelocity ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int k=minZ; k < maxZ; k++) for (int j=0; j < _maxY; j++) for (int i=0; i < _maxX; i++) op(i,j,k,phiObs,vel,center,radius);  } } else { const int k=0; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int j=0; j < _maxY; j++) for (int i=0; i < _maxX; i++) op(i,j,k,phiObs,vel,center,radius);  } }  } Grid<Real> & phiObs; MACGrid& vel; const Vec3& center; const Real& radius;   };
#line 302 "plugin/initplugins.cpp"



void initVortexVelocity(Grid<Real> &phiObs, MACGrid& vel, const Vec3 &center, const Real &radius) {
	kninitVortexVelocity(phiObs,  vel, center, radius);
} static PyObject* _W_15 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "initVortexVelocity" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; Grid<Real> & phiObs = *_args.getPtr<Grid<Real>  >("phiObs",0,&_lock); MACGrid& vel = *_args.getPtr<MACGrid >("vel",1,&_lock); const Vec3& center = _args.get<Vec3 >("center",2,&_lock); const Real& radius = _args.get<Real >("radius",3,&_lock);   _retval = getPyNone(); initVortexVelocity(phiObs,vel,center,radius);  _args.check(); } pbFinalizePlugin(parent,"initVortexVelocity", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("initVortexVelocity",e.what()); return 0; } } static const Pb::Register _RP_initVortexVelocity ("","initVortexVelocity",_W_15); 

inline static Real calcFraction(Real phi1, Real phi2)
{
	if(phi1>0. && phi2>0.) return 1.;
	if(phi1<0. && phi2<0.) return 0.;

	// make sure phi1 < phi2
	if (phi2<phi1) { Real t = phi1; phi1= phi2; phi2 = t; }
	Real denom = phi1-phi2;
	if (denom > -1e-04) return 0.5; 

	Real frac = 1. - phi1/denom;
	if(frac<0.01) frac = 0.; // skip , dont mark as fluid
	return std::max(Real(0), std::min(Real(1), frac ));
}


 struct KnUpdateFractions : public KernelBase { KnUpdateFractions(FlagGrid& flags, Grid<Real>& phiObs, MACGrid& fractions, const int &boundaryWidth) :  KernelBase(&flags,1) ,flags(flags),phiObs(phiObs),fractions(fractions),boundaryWidth(boundaryWidth)   { runMessage(); run(); }  inline void op(int i, int j, int k, FlagGrid& flags, Grid<Real>& phiObs, MACGrid& fractions, const int &boundaryWidth )  {

	// walls at domain bounds and inner objects
	fractions(i,j,k).x = calcFraction( phiObs(i,j,k) , phiObs(i-1,j,k));
	fractions(i,j,k).y = calcFraction( phiObs(i,j,k) , phiObs(i,j-1,k));
    if(phiObs.is3D()) {
	fractions(i,j,k).z = calcFraction( phiObs(i,j,k) , phiObs(i,j,k-1));
	}

	// remaining BCs at the domain boundaries 
	const int w = boundaryWidth;
	// only set if not in obstacle
 	if(phiObs(i,j,k)<0.) return;

	// x-direction boundaries
	if(i <= w+1) {                     //min x
		if( (flags.isInflow(i-1,j,k)) ||
			(flags.isOutflow(i-1,j,k)) ||
			(flags.isOpen(i-1,j,k)) ) {
				fractions(i,j,k).x = fractions(i,j,k).y = 1.; if(flags.is3D()) fractions(i,j,k).z = 1.;
		}
	}
	if(i >= flags.getSizeX()-w-2) {    //max x
		if(	(flags.isInflow(i+1,j,k)) ||
			(flags.isOutflow(i+1,j,k)) ||
			(flags.isOpen(i+1,j,k)) ) {
			fractions(i+1,j,k).x = fractions(i+1,j,k).y = 1.; if(flags.is3D()) fractions(i+1,j,k).z = 1.;
		}
	}
	// y-direction boundaries
 	if(j <= w+1) {                     //min y
		if(	(flags.isInflow(i,j-1,k)) ||
			(flags.isOutflow(i,j-1,k)) ||
			(flags.isOpen(i,j-1,k)) ) {
			fractions(i,j,k).x = fractions(i,j,k).y = 1.; if(flags.is3D()) fractions(i,j,k).z = 1.;
		}
 	}
 	if(j >= flags.getSizeY()-w-2) {      //max y
		if(	(flags.isInflow(i,j+1,k)) ||
			(flags.isOutflow(i,j+1,k)) ||
			(flags.isOpen(i,j+1,k)) ) {
			fractions(i,j+1,k).x = fractions(i,j+1,k).y = 1.; if(flags.is3D()) fractions(i,j+1,k).z = 1.;
		}
 	}
	// z-direction boundaries
	if(flags.is3D()) {
	if(k <= w+1) {                 //min z
		if(	(flags.isInflow(i,j,k-1)) ||
			(flags.isOutflow(i,j,k-1)) ||
			(flags.isOpen(i,j,k-1)) ) {
			fractions(i,j,k).x = fractions(i,j,k).y = 1.; if(flags.is3D()) fractions(i,j,k).z = 1.;
		}
	}
	if(j >= flags.getSizeZ()-w-2) { //max z
		if(	(flags.isInflow(i,j,k+1)) ||
			(flags.isOutflow(i,j,k+1)) ||
			(flags.isOpen(i,j,k+1)) ) {
			fractions(i,j,k+1).x = fractions(i,j,k+1).y = 1.; if(flags.is3D()) fractions(i,j,k+1).z = 1.;
		}
	}
	}

}   inline FlagGrid& getArg0() { return flags; } typedef FlagGrid type0;inline Grid<Real>& getArg1() { return phiObs; } typedef Grid<Real> type1;inline MACGrid& getArg2() { return fractions; } typedef MACGrid type2;inline const int& getArg3() { return boundaryWidth; } typedef int type3; void runMessage() { debMsg("Executing kernel KnUpdateFractions ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int k=minZ; k < maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,flags,phiObs,fractions,boundaryWidth);  } } else { const int k=0; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,flags,phiObs,fractions,boundaryWidth);  } }  } FlagGrid& flags; Grid<Real>& phiObs; MACGrid& fractions; const int& boundaryWidth;   };
#line 344 "plugin/initplugins.cpp"



void updateFractions(FlagGrid& flags, Grid<Real>& phiObs, MACGrid& fractions, const int &boundaryWidth=0) {
	fractions.setConst( Vec3(0.) );
	KnUpdateFractions(flags, phiObs, fractions, boundaryWidth);
} static PyObject* _W_16 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "updateFractions" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; FlagGrid& flags = *_args.getPtr<FlagGrid >("flags",0,&_lock); Grid<Real>& phiObs = *_args.getPtr<Grid<Real> >("phiObs",1,&_lock); MACGrid& fractions = *_args.getPtr<MACGrid >("fractions",2,&_lock); const int& boundaryWidth = _args.getOpt<int >("boundaryWidth",3,0,&_lock);   _retval = getPyNone(); updateFractions(flags,phiObs,fractions,boundaryWidth);  _args.check(); } pbFinalizePlugin(parent,"updateFractions", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("updateFractions",e.what()); return 0; } } static const Pb::Register _RP_updateFractions ("","updateFractions",_W_16); 


 struct KnUpdateFlags : public KernelBase { KnUpdateFlags(FlagGrid& flags, MACGrid& fractions, Grid<Real>& phiObs) :  KernelBase(&flags,1) ,flags(flags),fractions(fractions),phiObs(phiObs)   { runMessage(); run(); }  inline void op(int i, int j, int k, FlagGrid& flags, MACGrid& fractions, Grid<Real>& phiObs )  {

	Real test = 0.;
	test += fractions.get(i  ,j,k).x;
	test += fractions.get(i+1,j,k).x;
	test += fractions.get(i,j  ,k).y;
	test += fractions.get(i,j+1,k).y;
	if (flags.is3D()) {
	test += fractions.get(i,j,k  ).z;
	test += fractions.get(i,j,k+1).z; }

	if(test==0. && phiObs(i,j,k) < 0.) flags(i,j,k) = FlagGrid::TypeObstacle; 
	else flags(i,j,k) = FlagGrid::TypeEmpty; 
}   inline FlagGrid& getArg0() { return flags; } typedef FlagGrid type0;inline MACGrid& getArg1() { return fractions; } typedef MACGrid type1;inline Grid<Real>& getArg2() { return phiObs; } typedef Grid<Real> type2; void runMessage() { debMsg("Executing kernel KnUpdateFlags ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int k=minZ; k < maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,flags,fractions,phiObs);  } } else { const int k=0; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,flags,fractions,phiObs);  } }  } FlagGrid& flags; MACGrid& fractions; Grid<Real>& phiObs;   };
#line 414 "plugin/initplugins.cpp"



void setObstacleFlags(FlagGrid& flags, MACGrid& fractions, Grid<Real>& phiObs) {
	KnUpdateFlags(flags,fractions, phiObs);
} static PyObject* _W_17 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "setObstacleFlags" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; FlagGrid& flags = *_args.getPtr<FlagGrid >("flags",0,&_lock); MACGrid& fractions = *_args.getPtr<MACGrid >("fractions",1,&_lock); Grid<Real>& phiObs = *_args.getPtr<Grid<Real> >("phiObs",2,&_lock);   _retval = getPyNone(); setObstacleFlags(flags,fractions,phiObs);  _args.check(); } pbFinalizePlugin(parent,"setObstacleFlags", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("setObstacleFlags",e.what()); return 0; } } static const Pb::Register _RP_setObstacleFlags ("","setObstacleFlags",_W_17); 

} // namespace



