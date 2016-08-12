




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/Users/sbarschkis/Developer/Mantaflow/blenderIntegration/mantaflowgit/source/plugin/pressure.cpp"
/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Plugins for pressure correction: solve_pressure, and ghost fluid helpers
 *
 ******************************************************************************/
#include "vectorbase.h"
#include "kernel.h"
#include "conjugategrad.h"

using namespace std;
namespace Manta {

//! Kernel: Construct the right-hand side of the poisson equation



 struct MakeRhs : public KernelBase { MakeRhs(FlagGrid& flags, Grid<Real>& rhs, MACGrid& vel, Grid<Real>* perCellCorr, MACGrid* fractions) :  KernelBase(&flags,1) ,flags(flags),rhs(rhs),vel(vel),perCellCorr(perCellCorr),fractions(fractions) ,cnt(0),sum(0)  { runMessage(); run(); }  inline void op(int i, int j, int k, FlagGrid& flags, Grid<Real>& rhs, MACGrid& vel, Grid<Real>* perCellCorr, MACGrid* fractions ,int& cnt,double& sum)  {
	if (!flags.isFluid(i,j,k)) {
		rhs(i,j,k) = 0;
		return;
	}

	// compute divergence 
	// no flag checks: assumes vel at obstacle interfaces is set to zero
	Real set(0);
	if(!fractions) {
		set =               vel(i,j,k).x - vel(i+1,j,k).x + 
				 			vel(i,j,k).y - vel(i,j+1,k).y; 
		if(vel.is3D()) set+=vel(i,j,k).z - vel(i,j,k+1).z;
	}else{
		set =               (*fractions)(i,j,k).x * vel(i,j,k).x - (*fractions)(i+1,j,k).x * vel(i+1,j,k).x + 
							(*fractions)(i,j,k).y * vel(i,j,k).y - (*fractions)(i,j+1,k).y * vel(i,j+1,k).y; 
		if(vel.is3D()) set+=(*fractions)(i,j,k).z * vel(i,j,k).z - (*fractions)(i,j,k+1).z * vel(i,j,k+1).z;
	}
	
	// per cell divergence correction (optional)
	if(perCellCorr) 
		set += perCellCorr->get(i,j,k);
	
	// obtain sum, cell count
	sum += set;
	cnt++;
	
	rhs(i,j,k) = set;
}   inline FlagGrid& getArg0() { return flags; } typedef FlagGrid type0;inline Grid<Real>& getArg1() { return rhs; } typedef Grid<Real> type1;inline MACGrid& getArg2() { return vel; } typedef MACGrid type2;inline Grid<Real>* getArg3() { return perCellCorr; } typedef Grid<Real> type3;inline MACGrid* getArg4() { return fractions; } typedef MACGrid type4; void runMessage() { debMsg("Executing kernel MakeRhs ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 {  int cnt = 0;double sum = 0; 
#pragma omp for nowait 
  for (int k=minZ; k < maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,flags,rhs,vel,perCellCorr,fractions,cnt,sum); 
#pragma omp critical
{this->cnt += cnt; this->sum += sum; } } } else { const int k=0; 
#pragma omp parallel 
 {  int cnt = 0;double sum = 0; 
#pragma omp for nowait 
  for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,flags,rhs,vel,perCellCorr,fractions,cnt,sum); 
#pragma omp critical
{this->cnt += cnt; this->sum += sum; } } }  } FlagGrid& flags; Grid<Real>& rhs; MACGrid& vel; Grid<Real>* perCellCorr; MACGrid* fractions;  int cnt; double sum;  };
#line 24 "plugin/pressure.cpp"



//! Kernel: Apply velocity update from poisson equation


 struct CorrectVelocity : public KernelBase { CorrectVelocity(FlagGrid& flags, MACGrid& vel, Grid<Real>& pressure) :  KernelBase(&flags,1) ,flags(flags),vel(vel),pressure(pressure)   { runMessage(); run(); }  inline void op(int i, int j, int k, FlagGrid& flags, MACGrid& vel, Grid<Real>& pressure )  {
	IndexInt idx = flags.index(i,j,k);
	if (flags.isFluid(idx))
	{
		if (flags.isFluid(i-1,j,k)) vel[idx].x -= (pressure[idx] - pressure(i-1,j,k));
		if (flags.isFluid(i,j-1,k)) vel[idx].y -= (pressure[idx] - pressure(i,j-1,k));
		if (flags.is3D() && flags.isFluid(i,j,k-1)) vel[idx].z -= (pressure[idx] - pressure(i,j,k-1));
 
		if (flags.isEmpty(i-1,j,k)) vel[idx].x -= pressure[idx];
		if (flags.isEmpty(i,j-1,k)) vel[idx].y -= pressure[idx];
		if (flags.is3D() && flags.isEmpty(i,j,k-1)) vel[idx].z -= pressure[idx];
	}
	else if (flags.isEmpty(idx)&&!flags.isOutflow(idx)) // don't change velocities in outflow cells
	{
		if (flags.isFluid(i-1,j,k)) vel[idx].x += pressure(i-1,j,k);
		else                        vel[idx].x  = 0.f;
		if (flags.isFluid(i,j-1,k)) vel[idx].y += pressure(i,j-1,k);
		else                        vel[idx].y  = 0.f;
		if (flags.is3D() ) {
		if (flags.isFluid(i,j,k-1)) vel[idx].z += pressure(i,j,k-1);
		else                        vel[idx].z  = 0.f;
		}
	}
}   inline FlagGrid& getArg0() { return flags; } typedef FlagGrid type0;inline MACGrid& getArg1() { return vel; } typedef MACGrid type1;inline Grid<Real>& getArg2() { return pressure; } typedef Grid<Real> type2; void runMessage() { debMsg("Executing kernel CorrectVelocity ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int k=minZ; k < maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,flags,vel,pressure);  } } else { const int k=0; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,flags,vel,pressure);  } }  } FlagGrid& flags; MACGrid& vel; Grid<Real>& pressure;   };
#line 57 "plugin/pressure.cpp"



// *****************************************************************************
// Ghost fluid helpers

// calculate fraction filled with liquid (note, assumes inside value is < outside!)
inline static Real thetaHelper(Real inside, Real outside)
{
	Real denom = inside-outside;
	if (denom > -1e-04) return 0.5; // should always be neg, and large enough...
	return std::max(Real(0), std::min(Real(1), inside/denom));
}

// calculate ghost fluid factor, cell at idx should be a fluid cell
inline static Real ghostFluidHelper(IndexInt idx, int offset, const Grid<Real> &phi, Real gfClamp)
{
	Real alpha = thetaHelper(phi[idx], phi[idx+offset]);
	if (alpha < gfClamp) return alpha = gfClamp;
	return (1-(1/alpha)); 
}

//! Kernel: Adapt A0 for ghost fluid


 struct ApplyGhostFluidDiagonal : public KernelBase { ApplyGhostFluidDiagonal(Grid<Real> &A0, const FlagGrid &flags, const Grid<Real> &phi, Real gfClamp) :  KernelBase(&A0,1) ,A0(A0),flags(flags),phi(phi),gfClamp(gfClamp)   { runMessage(); run(); }  inline void op(int i, int j, int k, Grid<Real> &A0, const FlagGrid &flags, const Grid<Real> &phi, Real gfClamp )  {
	const int X = flags.getStrideX(), Y = flags.getStrideY(), Z = flags.getStrideZ();
	IndexInt idx = flags.index(i,j,k);
	if (!flags.isFluid(idx)) return;

	if (flags.isEmpty(i-1,j,k)) A0[idx] -= ghostFluidHelper(idx, -X, phi, gfClamp);
	if (flags.isEmpty(i+1,j,k)) A0[idx] -= ghostFluidHelper(idx, +X, phi, gfClamp);
	if (flags.isEmpty(i,j-1,k)) A0[idx] -= ghostFluidHelper(idx, -Y, phi, gfClamp);
	if (flags.isEmpty(i,j+1,k)) A0[idx] -= ghostFluidHelper(idx, +Y, phi, gfClamp);
	if (flags.is3D()) {
		if (flags.isEmpty(i,j,k-1)) A0[idx] -= ghostFluidHelper(idx, -Z, phi, gfClamp);
		if (flags.isEmpty(i,j,k+1)) A0[idx] -= ghostFluidHelper(idx, +Z, phi, gfClamp);
	}
}   inline Grid<Real> & getArg0() { return A0; } typedef Grid<Real>  type0;inline const FlagGrid& getArg1() { return flags; } typedef FlagGrid type1;inline const Grid<Real> & getArg2() { return phi; } typedef Grid<Real>  type2;inline Real& getArg3() { return gfClamp; } typedef Real type3; void runMessage() { debMsg("Executing kernel ApplyGhostFluidDiagonal ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int k=minZ; k < maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,A0,flags,phi,gfClamp);  } } else { const int k=0; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,A0,flags,phi,gfClamp);  } }  } Grid<Real> & A0; const FlagGrid& flags; const Grid<Real> & phi; Real gfClamp;   };
#line 104 "plugin/pressure.cpp"



//! Kernel: Apply velocity update: ghost fluid contribution


 struct CorrectVelocityGhostFluid : public KernelBase { CorrectVelocityGhostFluid(MACGrid &vel, const FlagGrid &flags, const Grid<Real> &pressure, const Grid<Real> &phi, Real gfClamp) :  KernelBase(&vel,1) ,vel(vel),flags(flags),pressure(pressure),phi(phi),gfClamp(gfClamp)   { runMessage(); run(); }  inline void op(int i, int j, int k, MACGrid &vel, const FlagGrid &flags, const Grid<Real> &pressure, const Grid<Real> &phi, Real gfClamp )  {
	const IndexInt X = flags.getStrideX(), Y = flags.getStrideY(), Z = flags.getStrideZ();
	const IndexInt idx = flags.index(i,j,k);
	if (flags.isFluid(idx))
	{
		if (flags.isEmpty(i-1,j,k)) vel[idx][0] += pressure[idx] * ghostFluidHelper(idx, -X, phi, gfClamp);
		if (flags.isEmpty(i,j-1,k)) vel[idx][1] += pressure[idx] * ghostFluidHelper(idx, -Y, phi, gfClamp);
		if (flags.is3D() && flags.isEmpty(i,j,k-1)) vel[idx][2] += pressure[idx] * ghostFluidHelper(idx, -Z, phi, gfClamp);
	}
	else if (flags.isEmpty(idx)&&!flags.isOutflow(idx)) // do not change velocities in outflow cells
	{
		if (flags.isFluid(i-1,j,k)) vel[idx][0] -= pressure(i-1,j,k) * ghostFluidHelper(idx-X, +X, phi, gfClamp);
		else                        vel[idx].x  = 0.f;
		if (flags.isFluid(i,j-1,k)) vel[idx][1] -= pressure(i,j-1,k) * ghostFluidHelper(idx-Y, +Y, phi, gfClamp);
		else                        vel[idx].y  = 0.f;
		if (flags.is3D() ) {
		if (flags.isFluid(i,j,k-1)) vel[idx][2] -= pressure(i,j,k-1) * ghostFluidHelper(idx-Z, +Z, phi, gfClamp);
		else                        vel[idx].z  = 0.f;
		}
	}
}   inline MACGrid& getArg0() { return vel; } typedef MACGrid type0;inline const FlagGrid& getArg1() { return flags; } typedef FlagGrid type1;inline const Grid<Real> & getArg2() { return pressure; } typedef Grid<Real>  type2;inline const Grid<Real> & getArg3() { return phi; } typedef Grid<Real>  type3;inline Real& getArg4() { return gfClamp; } typedef Real type4; void runMessage() { debMsg("Executing kernel CorrectVelocityGhostFluid ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int k=minZ; k < maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,vel,flags,pressure,phi,gfClamp);  } } else { const int k=0; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,vel,flags,pressure,phi,gfClamp);  } }  } MACGrid& vel; const FlagGrid& flags; const Grid<Real> & pressure; const Grid<Real> & phi; Real gfClamp;   };
#line 122 "plugin/pressure.cpp"




// improve behavior of clamping for large time steps:

inline static Real ghostFluidWasClamped(IndexInt idx, int offset, const Grid<Real> &phi, Real gfClamp)
{
	Real alpha = thetaHelper(phi[idx], phi[idx+offset]);
	if (alpha < gfClamp) return true;
	return false;
}




 struct ReplaceClampedGhostFluidVels : public KernelBase { ReplaceClampedGhostFluidVels(MACGrid &vel, FlagGrid &flags, const Grid<Real> &pressure, const Grid<Real> &phi, Real gfClamp ) :  KernelBase(&vel,1) ,vel(vel),flags(flags),pressure(pressure),phi(phi),gfClamp(gfClamp)   { runMessage(); run(); }  inline void op(int i, int j, int k, MACGrid &vel, FlagGrid &flags, const Grid<Real> &pressure, const Grid<Real> &phi, Real gfClamp  )  {
	const IndexInt idx = flags.index(i,j,k);
	const IndexInt X   = flags.getStrideX(), Y = flags.getStrideY(), Z = flags.getStrideZ();
	if (!flags.isEmpty(idx)) return;

	if( (flags.isFluid(i-1,j,k)) && ( ghostFluidWasClamped(idx-X, +X, phi, gfClamp) ) )
		vel[idx][0] = vel[idx-X][0];
	if( (flags.isFluid(i,j-1,k)) && ( ghostFluidWasClamped(idx-Y, +Y, phi, gfClamp) ) )
		vel[idx][1] = vel[idx-Y][1];
	if( flags.is3D() &&
	  ( (flags.isFluid(i,j,k-1)) && ( ghostFluidWasClamped(idx-Z, +Z, phi, gfClamp) ) ))
		vel[idx][2] = vel[idx-Z][2];

	if( (flags.isFluid(i+1,j,k)) && ( ghostFluidWasClamped(idx+X, -X, phi, gfClamp)) )
		vel[idx][0] = vel[idx+X][0];
	if( (flags.isFluid(i,j+1,k)) && ( ghostFluidWasClamped(idx+Y, -Y, phi, gfClamp)) )
		vel[idx][1] = vel[idx+Y][1];
	if( flags.is3D() && 
	   (flags.isFluid(i,j,k+1))  && ( ghostFluidWasClamped(idx+Z, -Z, phi, gfClamp)) )
		vel[idx][2] = vel[idx+Z][2];
}   inline MACGrid& getArg0() { return vel; } typedef MACGrid type0;inline FlagGrid& getArg1() { return flags; } typedef FlagGrid type1;inline const Grid<Real> & getArg2() { return pressure; } typedef Grid<Real>  type2;inline const Grid<Real> & getArg3() { return phi; } typedef Grid<Real>  type3;inline Real& getArg4() { return gfClamp; } typedef Real type4; void runMessage() { debMsg("Executing kernel ReplaceClampedGhostFluidVels ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int k=minZ; k < maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,vel,flags,pressure,phi,gfClamp);  } } else { const int k=0; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,vel,flags,pressure,phi,gfClamp);  } }  } MACGrid& vel; FlagGrid& flags; const Grid<Real> & pressure; const Grid<Real> & phi; Real gfClamp;   };
#line 157 "plugin/pressure.cpp"



//! Kernel: Compute min value of Real grid

 struct CountEmptyCells : public KernelBase { CountEmptyCells(FlagGrid& flags) :  KernelBase(&flags,0) ,flags(flags) ,numEmpty(0)  { runMessage(); run(); }   inline void op(IndexInt idx, FlagGrid& flags ,int& numEmpty)  {
	if (flags.isEmpty(idx) ) numEmpty++;
}    inline operator int () { return numEmpty; } inline int  & getRet() { return numEmpty; }  inline FlagGrid& getArg0() { return flags; } typedef FlagGrid type0; void runMessage() { debMsg("Executing kernel CountEmptyCells ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  int numEmpty = 0; 
#pragma omp for nowait 
  for (IndexInt i = 0; i < _sz; i++) op(i,flags,numEmpty); 
#pragma omp critical
{this->numEmpty += numEmpty; } }   } FlagGrid& flags;  int numEmpty;  };
#line 181 "plugin/pressure.cpp"



// *****************************************************************************
// Main pressure solve

IndexInt solvePressureBase(MACGrid& vel, Grid<Real>& pressure, FlagGrid& flags, Grid<Real>& rhs, Real cgAccuracy = 1e-3,
	Grid<Real>* phi = 0,
	Grid<Real>* perCellCorr = 0,
	MACGrid* fractions = 0,
	Real gfClamp = 1e-04,
	Real cgMaxIterFac = 1.5,
	bool precondition = true,
	bool enforceCompatibility = false,
	bool useL2Norm = false,
	bool zeroPressureFixing = false,
	Grid<Real>* retRhs = NULL)
{
	// reserve temp grids
	FluidSolver* parent = flags.getParent();
	Grid<Real> residual(parent);
	Grid<Real> search(parent);
	Grid<Real> A0(parent);
	Grid<Real> Ai(parent);
	Grid<Real> Aj(parent);
	Grid<Real> Ak(parent);
	Grid<Real> tmp(parent);
	Grid<Real> pca0(parent);
	Grid<Real> pca1(parent);
	Grid<Real> pca2(parent);
	Grid<Real> pca3(parent);
		
	// setup matrix and boundaries 
	MakeLaplaceMatrix (flags, A0, Ai, Aj, Ak, fractions);

	if (phi) {
		ApplyGhostFluidDiagonal(A0, flags, *phi, gfClamp);
	}
	
	// compute divergence and init right hand side
	MakeRhs kernMakeRhs (flags, rhs, vel, perCellCorr, fractions);
	
	if (enforceCompatibility)
		rhs += (Real)(-kernMakeRhs.sum / (Real)kernMakeRhs.cnt);
	
	// check whether we need to fix some pressure value...
	// (manually enable, or automatically for high accuracy, can cause asymmetries otherwise)
	IndexInt fixPidx = -1;
	if(zeroPressureFixing || cgAccuracy<1e-07) 
	{
		if(FLOATINGPOINT_PRECISION==1) debMsg("Warning - high CG accuracy with single-precision floating point accuracy might not converge...", 2);

		int numEmpty = CountEmptyCells(flags);
		if(numEmpty==0) {
			FOR_IJK_BND(flags,1) {
				if(flags.isFluid(i,j,k)) {
					fixPidx = flags.index(i,j,k);
					break;
				}
			}
			//debMsg("No empty cells! Fixing pressure of cell "<<fixPidx<<" to zero",1);
		}
		if(fixPidx>=0) {
			// adjustment for approx. symmetric zeroPressureFixing cell (top center)
			fixPidx = flags.index(flags.getSizeX() / 2, flags.getSizeY() - 2, flags.is3D() ? flags.getSizeZ() / 2 : 0);

			flags[fixPidx] |= FlagGrid::TypeZeroPressure;
			rhs[fixPidx] = 0.; 
			debMsg("Pinning pressure of cell "<<fixPidx<<" to zero", 2);
		}
	}

	// CG setup
	// note: the last factor increases the max iterations for 2d, which right now can't use a preconditioner 
	const int maxIter = (int)(cgMaxIterFac * flags.getSize().max()) * (flags.is3D() ? 1 : 4);
	GridCgInterface *gcg;
	if (vel.is3D())
		gcg = new GridCg<ApplyMatrix>  (pressure, rhs, residual, search, flags, tmp, &A0, &Ai, &Aj, &Ak );
	else
		gcg = new GridCg<ApplyMatrix2D>(pressure, rhs, residual, search, flags, tmp, &A0, &Ai, &Aj, &Ak );
	
	gcg->setAccuracy( cgAccuracy ); 
	gcg->setUseL2Norm( useL2Norm );

	// optional preconditioning
	gcg->setPreconditioner( precondition ? GridCgInterface::PC_mICP : GridCgInterface::PC_None, &pca0, &pca1, &pca2, &pca3);

	for (int iter=0; iter<maxIter; iter++) {
		if (!gcg->iterate()) iter=maxIter;
	} 
	//debMsg("FluidSolver::solvePressureBase iterations:"<<gcg->getIterations()<<", res:"<<gcg->getSigma(), 1);
	delete gcg;
	return fixPidx;
}

//! Perform pressure projection of the velocity grid











void solvePressure(MACGrid& vel, Grid<Real>& pressure, FlagGrid& flags, Real cgAccuracy = 1e-3, Grid<Real>* phi = 0, Grid<Real>* perCellCorr = 0, MACGrid* fractions = 0, Real gfClamp = 1e-04, Real cgMaxIterFac = 1.5, bool precondition = true, bool enforceCompatibility = false, bool useL2Norm = false, bool zeroPressureFixing = false, Grid<Real>* retRhs = NULL ) {
	Grid<Real> rhs(vel.getParent());
	IndexInt fixPidx = solvePressureBase(vel, pressure, flags, rhs, cgAccuracy, phi, perCellCorr, fractions, gfClamp, cgMaxIterFac, precondition, enforceCompatibility, useL2Norm, zeroPressureFixing, retRhs);
	
	CorrectVelocity(flags, vel, pressure ); 
	if (phi) {
		CorrectVelocityGhostFluid (vel, flags, pressure, *phi, gfClamp);
		// improve behavior of clamping for large time steps:
		ReplaceClampedGhostFluidVels (vel, flags, pressure, *phi, gfClamp);
	}

	if(fixPidx>=0)
		flags[fixPidx] &= ~FlagGrid::TypeZeroPressure;

	// optionally , return RHS
	if(retRhs) {
		retRhs->copyFrom( rhs );
	}
} static PyObject* _W_0 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "solvePressure" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; MACGrid& vel = *_args.getPtr<MACGrid >("vel",0,&_lock); Grid<Real>& pressure = *_args.getPtr<Grid<Real> >("pressure",1,&_lock); FlagGrid& flags = *_args.getPtr<FlagGrid >("flags",2,&_lock); Real cgAccuracy = _args.getOpt<Real >("cgAccuracy",3,1e-3,&_lock); Grid<Real>* phi = _args.getPtrOpt<Grid<Real> >("phi",4,0,&_lock); Grid<Real>* perCellCorr = _args.getPtrOpt<Grid<Real> >("perCellCorr",5,0,&_lock); MACGrid* fractions = _args.getPtrOpt<MACGrid >("fractions",6,0,&_lock); Real gfClamp = _args.getOpt<Real >("gfClamp",7,1e-04,&_lock); Real cgMaxIterFac = _args.getOpt<Real >("cgMaxIterFac",8,1.5,&_lock); bool precondition = _args.getOpt<bool >("precondition",9,true,&_lock); bool enforceCompatibility = _args.getOpt<bool >("enforceCompatibility",10,false,&_lock); bool useL2Norm = _args.getOpt<bool >("useL2Norm",11,false,&_lock); bool zeroPressureFixing = _args.getOpt<bool >("zeroPressureFixing",12,false,&_lock); Grid<Real>* retRhs = _args.getPtrOpt<Grid<Real> >("retRhs",13,NULL ,&_lock);   _retval = getPyNone(); solvePressure(vel,pressure,flags,cgAccuracy,phi,perCellCorr,fractions,gfClamp,cgMaxIterFac,precondition,enforceCompatibility,useL2Norm,zeroPressureFixing,retRhs);  _args.check(); } pbFinalizePlugin(parent,"solvePressure", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("solvePressure",e.what()); return 0; } } static const Pb::Register _RP_solvePressure ("","solvePressure",_W_0); 

} // end namespace



