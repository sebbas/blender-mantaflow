




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
 * Plugins for pressure correction: solve_pressure, and ghost fluid helpers
 *
 ******************************************************************************/
#include "vectorbase.h"
#include "kernel.h"
#include "conjugategrad.h"

using namespace std;
namespace Manta {

//! Kernel: Construct the right-hand side of the poisson equation



 struct MakeRhs : public KernelBase { MakeRhs(FlagGrid& flags, Grid<Real>& rhs, MACGrid& vel, Grid<Real>* perCellCorr, MACGrid* fractions) :  KernelBase(&flags,1) ,flags(flags),rhs(rhs),vel(vel),perCellCorr(perCellCorr),fractions(fractions) ,cnt(0),sum(0)  { run(); }  inline void op(int i, int j, int k, FlagGrid& flags, Grid<Real>& rhs, MACGrid& vel, Grid<Real>* perCellCorr, MACGrid* fractions ,int& cnt,double& sum)  {
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
		set =               fractions->get(i,j,k).x*vel(i,j,k).x - fractions->get(i+1,j,k).x*vel(i+1,j,k).x + 
							fractions->get(i,j,k).y*vel(i,j,k).y - fractions->get(i,j+1,k).y*vel(i,j+1,k).y; 
		if(vel.is3D()) set+=fractions->get(i,j,k).z*vel(i,j,k).z - fractions->get(i,j,k+1).z*vel(i,j,k+1).z;
	}
	
	// per cell divergence correction
	if(perCellCorr) 
		set += perCellCorr->get(i,j,k);
	
	// obtain sum, cell count
	sum += set;
	cnt++;
	
	rhs(i,j,k) = set;
}   inline FlagGrid& getArg0() { return flags; } typedef FlagGrid type0;inline Grid<Real>& getArg1() { return rhs; } typedef Grid<Real> type1;inline MACGrid& getArg2() { return vel; } typedef MACGrid type2;inline Grid<Real>* getArg3() { return perCellCorr; } typedef Grid<Real> type3;inline MACGrid* getArg4() { return fractions; } typedef MACGrid type4; void run() {  const int _maxX = maxX; const int _maxY = maxY; for (int k=minZ; k< maxZ; k++) for (int j=1; j< _maxY; j++) for (int i=1; i< _maxX; i++) op(i,j,k, flags,rhs,vel,perCellCorr,fractions,cnt,sum);  } FlagGrid& flags; Grid<Real>& rhs; MACGrid& vel; Grid<Real>* perCellCorr; MACGrid* fractions;  int cnt; double sum;  };

//! Kernel: Apply velocity update from poisson equation


 struct CorrectVelocity : public KernelBase { CorrectVelocity(FlagGrid& flags, MACGrid& vel, Grid<Real>& pressure) :  KernelBase(&flags,1) ,flags(flags),vel(vel),pressure(pressure)   { run(); }  inline void op(int i, int j, int k, FlagGrid& flags, MACGrid& vel, Grid<Real>& pressure )  {
	int idx = flags.index(i,j,k);
	if (flags.isFluid(idx))
	{
		if (flags.isFluid(i-1,j,k)) vel[idx].x -= (pressure[idx] - pressure(i-1,j,k));
		if (flags.isFluid(i,j-1,k)) vel[idx].y -= (pressure[idx] - pressure(i,j-1,k));
		if (flags.is3D() && flags.isFluid(i,j,k-1)) vel[idx].z -= (pressure[idx] - pressure(i,j,k-1));
 
		if (flags.isEmpty(i-1,j,k)) vel[idx].x -= pressure[idx];
		if (flags.isEmpty(i,j-1,k)) vel[idx].y -= pressure[idx];
		if (flags.is3D() && flags.isEmpty(i,j,k-1)) vel[idx].z -= pressure[idx];
	}
	else if (flags.isEmpty(idx))
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
}   inline FlagGrid& getArg0() { return flags; } typedef FlagGrid type0;inline MACGrid& getArg1() { return vel; } typedef MACGrid type1;inline Grid<Real>& getArg2() { return pressure; } typedef Grid<Real> type2; void run() {  const int _maxX = maxX; const int _maxY = maxY; for (int k=minZ; k< maxZ; k++) for (int j=1; j< _maxY; j++) for (int i=1; i< _maxX; i++) op(i,j,k, flags,vel,pressure);  } FlagGrid& flags; MACGrid& vel; Grid<Real>& pressure;   };




 struct SetOpenBound : public KernelBase { SetOpenBound(Grid<Real> &A0,Grid<Real> &Ai,Grid<Real> &Aj,Grid<Real> &Ak,FlagGrid& flags,MACGrid& vel, Vector3D<bool> lo, Vector3D<bool> up) :  KernelBase(&A0,0) ,A0(A0),Ai(Ai),Aj(Aj),Ak(Ak),flags(flags),vel(vel),lo(lo),up(up)   { run(); }  inline void op(int i, int j, int k, Grid<Real> &A0,Grid<Real> &Ai,Grid<Real> &Aj,Grid<Real> &Ak,FlagGrid& flags,MACGrid& vel, Vector3D<bool> lo, Vector3D<bool> up )  {

	if (!flags.isFluid(i,j,k))
		return;
	
	int b = flags.getBoundaryWidth();

	// set matrix stencil in and at boundary to empty
	if((lo.x && i <= b+1)||(up.x && i >= maxX-b-2)||(lo.y && j <= b+1)||(up.y && j >= maxY-b-2))
		A0(i,j,k) = (flags.is3D()) ? 6. : 4.;
	
	if ((lo.x && i <= b)||(up.x && i >= maxX-b-2))					 Ai(i,j,k) = .0;
	if ((lo.y && j <= b)||(up.y && j >= maxY-b-2))					 Aj(i,j,k) = .0;
	if (flags.is3D() && ((lo.z && k <= b)||(up.z && k >= maxZ-b-2))) Ak(i,j,k) = .0;

	// set velocity boundary conditions
	if (lo.x && i == b)				vel(b,j,k) = vel(b+1,j,k);
	if (lo.y && j == b)				vel(i,b,k) = vel(i,b+1,k);
	if (up.x && i == maxX-b-1)		vel(maxX-b-1,j,k) = vel(maxX-b-2,j,k);
	if (up.y && j == maxY-b-1)		vel(i,maxY-b-1,k) = vel(i,maxY-b-2,k);
	if(flags.is3D()) {
		if (lo.z && k == b)			vel(i,j,b) = vel(i,j,b+1);
		if (up.z && k == maxZ-b-1)	vel(i,j,maxZ-b-1) = vel(i,j,maxZ-b-2); 
	}
}   inline Grid<Real> & getArg0() { return A0; } typedef Grid<Real>  type0;inline Grid<Real> & getArg1() { return Ai; } typedef Grid<Real>  type1;inline Grid<Real> & getArg2() { return Aj; } typedef Grid<Real>  type2;inline Grid<Real> & getArg3() { return Ak; } typedef Grid<Real>  type3;inline FlagGrid& getArg4() { return flags; } typedef FlagGrid type4;inline MACGrid& getArg5() { return vel; } typedef MACGrid type5;inline Vector3D<bool> & getArg6() { return lo; } typedef Vector3D<bool>  type6;inline Vector3D<bool> & getArg7() { return up; } typedef Vector3D<bool>  type7; void run() {  const int _maxX = maxX; const int _maxY = maxY; for (int k=minZ; k< maxZ; k++) for (int j=0; j< _maxY; j++) for (int i=0; i< _maxX; i++) op(i,j,k, A0,Ai,Aj,Ak,flags,vel,lo,up);  } Grid<Real> & A0; Grid<Real> & Ai; Grid<Real> & Aj; Grid<Real> & Ak; FlagGrid& flags; MACGrid& vel; Vector3D<bool>  lo; Vector3D<bool>  up;   };


//! Kernel: Set matrix rhs for outflow

 struct SetOutflow : public KernelBase { SetOutflow(Grid<Real>& rhs, Vector3D<bool> lowerBound, Vector3D<bool> upperBound, int height) :  KernelBase(&rhs,0) ,rhs(rhs),lowerBound(lowerBound),upperBound(upperBound),height(height)   { run(); }  inline void op(int i, int j, int k, Grid<Real>& rhs, Vector3D<bool> lowerBound, Vector3D<bool> upperBound, int height )  {
	if ((lowerBound.x && i < height) || (upperBound.x && i >= maxX-1-height) ||
		(lowerBound.y && j < height) || (upperBound.y && j >= maxY-1-height) ||
		(lowerBound.z && k < height) || (upperBound.z && k >= maxZ-1-height))
		rhs(i,j,k) = 0;
}   inline Grid<Real>& getArg0() { return rhs; } typedef Grid<Real> type0;inline Vector3D<bool> & getArg1() { return lowerBound; } typedef Vector3D<bool>  type1;inline Vector3D<bool> & getArg2() { return upperBound; } typedef Vector3D<bool>  type2;inline int& getArg3() { return height; } typedef int type3; void run() {  const int _maxX = maxX; const int _maxY = maxY; for (int k=minZ; k< maxZ; k++) for (int j=0; j< _maxY; j++) for (int i=0; i< _maxX; i++) op(i,j,k, rhs,lowerBound,upperBound,height);  } Grid<Real>& rhs; Vector3D<bool>  lowerBound; Vector3D<bool>  upperBound; int height;   };


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
inline static Real ghostFluidHelper(int idx, int offset, const Grid<Real> &phi, Real gfClamp)
{
	Real alpha = thetaHelper(phi[idx], phi[idx+offset]);
	if (alpha < gfClamp) return alpha = gfClamp;
	return (1-(1/alpha)); 
}

//! Kernel: Adapt A0 for ghost fluid


 struct ApplyGhostFluidDiagonal : public KernelBase { ApplyGhostFluidDiagonal(Grid<Real> &A0, const FlagGrid &flags, const Grid<Real> &phi, Real gfClamp) :  KernelBase(&A0,1) ,A0(A0),flags(flags),phi(phi),gfClamp(gfClamp)   { run(); }  inline void op(int i, int j, int k, Grid<Real> &A0, const FlagGrid &flags, const Grid<Real> &phi, Real gfClamp )  {
	const int X = flags.getStrideX(), Y = flags.getStrideY(), Z = flags.getStrideZ();
	int idx = flags.index(i,j,k);
	if (!flags.isFluid(idx)) return;

	if (flags.isEmpty(i-1,j,k)) A0[idx] -= ghostFluidHelper(idx, -X, phi, gfClamp);
	if (flags.isEmpty(i+1,j,k)) A0[idx] -= ghostFluidHelper(idx, +X, phi, gfClamp);
	if (flags.isEmpty(i,j-1,k)) A0[idx] -= ghostFluidHelper(idx, -Y, phi, gfClamp);
	if (flags.isEmpty(i,j+1,k)) A0[idx] -= ghostFluidHelper(idx, +Y, phi, gfClamp);
	if (flags.is3D()) {
		if (flags.isEmpty(i,j,k-1)) A0[idx] -= ghostFluidHelper(idx, -Z, phi, gfClamp);
		if (flags.isEmpty(i,j,k+1)) A0[idx] -= ghostFluidHelper(idx, +Z, phi, gfClamp);
	}
}   inline Grid<Real> & getArg0() { return A0; } typedef Grid<Real>  type0;inline const FlagGrid& getArg1() { return flags; } typedef FlagGrid type1;inline const Grid<Real> & getArg2() { return phi; } typedef Grid<Real>  type2;inline Real& getArg3() { return gfClamp; } typedef Real type3; void run() {  const int _maxX = maxX; const int _maxY = maxY; for (int k=minZ; k< maxZ; k++) for (int j=1; j< _maxY; j++) for (int i=1; i< _maxX; i++) op(i,j,k, A0,flags,phi,gfClamp);  } Grid<Real> & A0; const FlagGrid& flags; const Grid<Real> & phi; Real gfClamp;   };

//! Kernel: Apply velocity update: ghost fluid contribution


 struct CorrectVelocityGhostFluid : public KernelBase { CorrectVelocityGhostFluid(MACGrid &vel, const FlagGrid &flags, const Grid<Real> &pressure, const Grid<Real> &phi, Real gfClamp) :  KernelBase(&vel,1) ,vel(vel),flags(flags),pressure(pressure),phi(phi),gfClamp(gfClamp)   { run(); }  inline void op(int i, int j, int k, MACGrid &vel, const FlagGrid &flags, const Grid<Real> &pressure, const Grid<Real> &phi, Real gfClamp )  {
	const int X = flags.getStrideX(), Y = flags.getStrideY(), Z = flags.getStrideZ();
	const int idx = flags.index(i,j,k);
	if (flags.isFluid(idx))
	{
		if (flags.isEmpty(i-1,j,k)) vel[idx][0] += pressure[idx] * ghostFluidHelper(idx, -X, phi, gfClamp);
		if (flags.isEmpty(i,j-1,k)) vel[idx][1] += pressure[idx] * ghostFluidHelper(idx, -Y, phi, gfClamp);
		if (flags.is3D() && flags.isEmpty(i,j,k-1)) vel[idx][2] += pressure[idx] * ghostFluidHelper(idx, -Z, phi, gfClamp);
	}
	else if (flags.isEmpty(idx))
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
}   inline MACGrid& getArg0() { return vel; } typedef MACGrid type0;inline const FlagGrid& getArg1() { return flags; } typedef FlagGrid type1;inline const Grid<Real> & getArg2() { return pressure; } typedef Grid<Real>  type2;inline const Grid<Real> & getArg3() { return phi; } typedef Grid<Real>  type3;inline Real& getArg4() { return gfClamp; } typedef Real type4; void run() {  const int _maxX = maxX; const int _maxY = maxY; for (int k=minZ; k< maxZ; k++) for (int j=1; j< _maxY; j++) for (int i=1; i< _maxX; i++) op(i,j,k, vel,flags,pressure,phi,gfClamp);  } MACGrid& vel; const FlagGrid& flags; const Grid<Real> & pressure; const Grid<Real> & phi; Real gfClamp;   };


// improve behavior of clamping for large time steps:

inline static Real ghostFluidWasClamped(int idx, int offset, const Grid<Real> &phi, Real gfClamp)
{
	Real alpha = thetaHelper(phi[idx], phi[idx+offset]);
	if (alpha < gfClamp) return true;
	return false;
}




 struct ReplaceClampedGhostFluidVels : public KernelBase { ReplaceClampedGhostFluidVels(MACGrid &vel, FlagGrid &flags, const Grid<Real> &pressure, const Grid<Real> &phi, Real gfClamp ) :  KernelBase(&vel,1) ,vel(vel),flags(flags),pressure(pressure),phi(phi),gfClamp(gfClamp)   { run(); }  inline void op(int i, int j, int k, MACGrid &vel, FlagGrid &flags, const Grid<Real> &pressure, const Grid<Real> &phi, Real gfClamp  )  {
	const int X = flags.getStrideX(), Y = flags.getStrideY(), Z = flags.getStrideZ();
	const int idx = flags.index(i,j,k);
	if (flags.isFluid(idx))
	{
		if( (flags.isEmpty(i-1,j,k)) && (ghostFluidWasClamped(idx, -X, phi, gfClamp)) )
			vel[idx-X][0] = vel[idx][0];
		if( (flags.isEmpty(i,j-1,k)) && (ghostFluidWasClamped(idx, -Y, phi, gfClamp)) )
			vel[idx-Y][1] = vel[idx][1];
		if( flags.is3D() && 
		   (flags.isEmpty(i,j,k-1)) && (ghostFluidWasClamped(idx, -Z, phi, gfClamp)) )
			vel[idx-Z][2] = vel[idx][2];
	}
	else if (flags.isEmpty(idx))
	{
		if( (i>-1) && (flags.isFluid(i-1,j,k)) && ( ghostFluidWasClamped(idx-X, +X, phi, gfClamp) ) )
			vel[idx][0] = vel[idx-X][0];
		if( (j>-1) && (flags.isFluid(i,j-1,k)) && ( ghostFluidWasClamped(idx-Y, +Y, phi, gfClamp) ) )
			vel[idx][1] = vel[idx-Y][1];
		if( flags.is3D() &&
		  ( (k>-1) && (flags.isFluid(i,j,k-1)) && ( ghostFluidWasClamped(idx-Z, +Z, phi, gfClamp) ) ))
			vel[idx][2] = vel[idx-Z][2];
	}
}   inline MACGrid& getArg0() { return vel; } typedef MACGrid type0;inline FlagGrid& getArg1() { return flags; } typedef FlagGrid type1;inline const Grid<Real> & getArg2() { return pressure; } typedef Grid<Real>  type2;inline const Grid<Real> & getArg3() { return phi; } typedef Grid<Real>  type3;inline Real& getArg4() { return gfClamp; } typedef Real type4; void run() {  const int _maxX = maxX; const int _maxY = maxY; for (int k=minZ; k< maxZ; k++) for (int j=1; j< _maxY; j++) for (int i=1; i< _maxX; i++) op(i,j,k, vel,flags,pressure,phi,gfClamp);  } MACGrid& vel; FlagGrid& flags; const Grid<Real> & pressure; const Grid<Real> & phi; Real gfClamp;   };

//! Kernel: Compute min value of Real grid

 struct CountEmptyCells : public KernelBase { CountEmptyCells(FlagGrid& flags) :  KernelBase(&flags,0) ,flags(flags) ,numEmpty(0)  { run(); }  inline void op(int idx, FlagGrid& flags ,int& numEmpty)  {
	if (flags.isEmpty(idx) ) numEmpty++;
}   inline operator int () { return numEmpty; } inline int  & getRet() { return numEmpty; }  inline FlagGrid& getArg0() { return flags; } typedef FlagGrid type0; void run() {  const int _sz = size; for (int i=0; i < _sz; i++) op(i, flags,numEmpty);  } FlagGrid& flags;  int numEmpty;  };

// *****************************************************************************
// Main pressure solve

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

 struct KnupdateFractions : public KernelBase { KnupdateFractions(FlagGrid& flags, Grid<Real>& phi, MACGrid& fractions) :  KernelBase(&flags,1) ,flags(flags),phi(phi),fractions(fractions)   { run(); }  inline void op(int i, int j, int k, FlagGrid& flags, Grid<Real>& phi, MACGrid& fractions )  {

	fractions(i,j,k).x = fractions(i,j,k).y = fractions(i,j,k).z = static_cast<float>(flags(i,j,k) & 1);

	Real tmp = 0.;
	tmp = fabs((phi(i,j,k) + phi(i-1,j,k))/2.);
    if(tmp < 1.0 && !((phi(i,j,k)<0) == (phi(i-1,j,k)<0)) ) {
    	fractions(i,j,k).x = tmp;
    }else if( (flags(i-1,j,k) == 2 && flags(i,j,k) == 1) || (flags(i,j,k) == 2 && flags(i-1,j,k) == 1)) {
		fractions(i,j,k).x = 0.;
    }
    tmp = fabs((phi(i,j,k) + phi(i,j-1,k))/2.);
    if(tmp < 1.0 && !((phi(i,j,k)<0) == (phi(i,j-1,k)<0)) ) {
    	fractions(i,j,k).y = tmp;
    }else if( (flags(i,j-1,k) == 2 && flags(i,j,k) == 1) || (flags(i,j,k) == 2 && flags(i,j-1,k) == 1)) {
		fractions(i,j,k).y = 0.;
    }
    if(flags.is3D()) {
	    tmp = fabs((phi(i,j,k) + phi(i,j,k-1))/2.);
	    if(tmp < 1.0 && !((phi(i,j,k)<0) == (phi(i,j,k-1)<0)) ) {
	    	fractions(i,j,k).z = tmp;
	    }else if( (flags(i,j,k-1) == 2 && flags(i,j,k) == 1) || (flags(i,j,k) == 2 && flags(i,j,k-1) == 1)) {
			fractions(i,j,k).z = 0.;
	    }
	}

}   inline FlagGrid& getArg0() { return flags; } typedef FlagGrid type0;inline Grid<Real>& getArg1() { return phi; } typedef Grid<Real> type1;inline MACGrid& getArg2() { return fractions; } typedef MACGrid type2; void run() {  const int _maxX = maxX; const int _maxY = maxY; for (int k=minZ; k< maxZ; k++) for (int j=1; j< _maxY; j++) for (int i=1; i< _maxX; i++) op(i,j,k, flags,phi,fractions);  } FlagGrid& flags; Grid<Real>& phi; MACGrid& fractions;   };

void updateFractions(FlagGrid& flags, Grid<Real>& phi, MACGrid& fractions) {
	KnupdateFractions(flags, phi, fractions);
} static PyObject* _W_0 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "updateFractions" ); PyObject *_retval = 0; { ArgLocker _lock; FlagGrid& flags = *_args.getPtr<FlagGrid >("flags",0,&_lock); Grid<Real>& phi = *_args.getPtr<Grid<Real> >("phi",1,&_lock); MACGrid& fractions = *_args.getPtr<MACGrid >("fractions",2,&_lock);   _retval = getPyNone(); updateFractions(flags,phi,fractions);  _args.check(); } pbFinalizePlugin(parent,"updateFractions" ); return _retval; } catch(std::exception& e) { pbSetError("updateFractions",e.what()); return 0; } } static const Pb::Register _RP_updateFractions ("","updateFractions",_W_0); 

//! Perform pressure projection of the velocity grid













void solvePressure(MACGrid& vel, Grid<Real>& pressure, FlagGrid& flags, string openBound="", Grid<Real>* phi = 0, Grid<Real>* perCellCorr = 0, MACGrid* fractions = 0, Real gfClamp = 1e-04, Real cgMaxIterFac = 1.5, Real cgAccuracy = 1e-3, string outflow = "", int outflowHeight = 1, bool precondition = true, bool enforceCompatibility = false, bool useL2Norm = false, Grid<Real>* retRhs = NULL ) {
	// parse strings
	Vector3D<bool> loOpenBound, upOpenBound, loOutflow, upOutflow;
	convertDescToVec(openBound, loOpenBound, upOpenBound);
	convertDescToVec(outflow, loOutflow, upOutflow);
	if (vel.is2D() && (loOpenBound.z || upOpenBound.z))
		errMsg("open boundaries for z specified for 2D grid");
	
	// reserve temp grids
	FluidSolver* parent = flags.getParent();
	Grid<Real> rhs(parent);
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

	SetOpenBound (A0, Ai, Aj, Ak, flags, vel, loOpenBound, upOpenBound);
	
	if (phi) {
		ApplyGhostFluidDiagonal(A0, flags, *phi, gfClamp);
	}
	
	// compute divergence and init right hand side
	MakeRhs kernMakeRhs (flags, rhs, vel, perCellCorr, fractions);
	
	if (!outflow.empty())
		SetOutflow (rhs, loOutflow, upOutflow, outflowHeight);
	
	if (enforceCompatibility)
		rhs += (Real)(-kernMakeRhs.sum / (Real)kernMakeRhs.cnt);
	
	// check whether we need to fix some pressure value...
	int fixPidx = -1;
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
		flags[fixPidx] |= FlagGrid::TypeZeroPressure;
		rhs[fixPidx] = 0.; 
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
	debMsg("FluidSolver::solvePressure iterations:"<<gcg->getIterations()<<", res:"<<gcg->getSigma(), 1);
	delete gcg;
	
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
} static PyObject* _W_1 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "solvePressure" ); PyObject *_retval = 0; { ArgLocker _lock; MACGrid& vel = *_args.getPtr<MACGrid >("vel",0,&_lock); Grid<Real>& pressure = *_args.getPtr<Grid<Real> >("pressure",1,&_lock); FlagGrid& flags = *_args.getPtr<FlagGrid >("flags",2,&_lock); string openBound = _args.getOpt<string >("openBound",3,"",&_lock); Grid<Real>* phi = _args.getPtrOpt<Grid<Real> >("phi",4,0,&_lock); Grid<Real>* perCellCorr = _args.getPtrOpt<Grid<Real> >("perCellCorr",5,0,&_lock); MACGrid* fractions = _args.getPtrOpt<MACGrid >("fractions",6,0,&_lock); Real gfClamp = _args.getOpt<Real >("gfClamp",7,1e-04,&_lock); Real cgMaxIterFac = _args.getOpt<Real >("cgMaxIterFac",8,1.5,&_lock); Real cgAccuracy = _args.getOpt<Real >("cgAccuracy",9,1e-3,&_lock); string outflow = _args.getOpt<string >("outflow",10,"",&_lock); int outflowHeight = _args.getOpt<int >("outflowHeight",11,1,&_lock); bool precondition = _args.getOpt<bool >("precondition",12,true,&_lock); bool enforceCompatibility = _args.getOpt<bool >("enforceCompatibility",13,false,&_lock); bool useL2Norm = _args.getOpt<bool >("useL2Norm",14,false,&_lock); Grid<Real>* retRhs = _args.getPtrOpt<Grid<Real> >("retRhs",15,NULL ,&_lock);   _retval = getPyNone(); solvePressure(vel,pressure,flags,openBound,phi,perCellCorr,fractions,gfClamp,cgMaxIterFac,cgAccuracy,outflow,outflowHeight,precondition,enforceCompatibility,useL2Norm,retRhs);  _args.check(); } pbFinalizePlugin(parent,"solvePressure" ); return _retval; } catch(std::exception& e) { pbSetError("solvePressure",e.what()); return 0; } } static const Pb::Register _RP_solvePressure ("","solvePressure",_W_1); 

} // end namespace



