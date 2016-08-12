




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/Users/sbarschkis/Developer/Mantaflow/blenderIntegration/mantaflowgit/source/conjugategrad.h"
/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Conjugate gradient solver
 *
 ******************************************************************************/

#ifndef _CONJUGATEGRADIENT_H
#define _CONJUGATEGRADIENT_H

#include "vectorbase.h"
#include "grid.h"
#include "kernel.h"


namespace Manta { 

static const bool CG_DEBUG = false;

//! Basic CG interface 
class GridCgInterface {
	public:
		enum PreconditionType { PC_None=0, PC_ICP, PC_mICP };
		
		GridCgInterface() : mUseL2Norm(true) {};
		virtual ~GridCgInterface() {};

		// solving functions
		virtual bool iterate() = 0;
		virtual void solve(int maxIter) = 0;

		// precond
		virtual void setPreconditioner(PreconditionType method, Grid<Real> *A0, Grid<Real> *Ai, Grid<Real> *Aj, Grid<Real> *Ak) = 0;

		// access
		virtual Real getSigma() const = 0;
		virtual Real getIterations() const = 0;
		virtual Real getResNorm() const = 0;
		virtual void setAccuracy(Real set) = 0;
		virtual Real getAccuracy() const = 0;

		void setUseL2Norm(bool set) { mUseL2Norm = set; }

	protected:

		// use l2 norm of residualfor threshold? (otherwise uses max norm)
		bool mUseL2Norm; 
};


//! Run single iteration of the cg solver
/*! the template argument determines the type of matrix multiplication,
	typically a ApplyMatrix kernel, another one is needed e.g. for the
	mesh-based wave equation solver */
template<class APPLYMAT>
class GridCg : public GridCgInterface {
	public:
		//! constructor
		GridCg(Grid<Real>& dst, Grid<Real>& rhs, Grid<Real>& residual, Grid<Real>& search, FlagGrid& flags, Grid<Real>& tmp, 
				Grid<Real>* A0, Grid<Real>* pAi, Grid<Real>* pAj, Grid<Real>* pAk);
		~GridCg() {}
		
		void doInit();
		bool iterate();
		void solve(int maxIter);
		//! init pointers, and copy values from "normal" matrix
		void setPreconditioner(PreconditionType method, Grid<Real> *A0, Grid<Real> *Ai, Grid<Real> *Aj, Grid<Real> *Ak);
		
		// Accessors        
		Real getSigma() const { return mSigma; }
		Real getIterations() const { return mIterations; }

		Real getResNorm() const { return mResNorm; }

		void setAccuracy(Real set) { mAccuracy=set; }
		Real getAccuracy() const { return mAccuracy; }

	protected:
		bool mInited;
		int mIterations;
		// grids
		Grid<Real>& mDst;
		Grid<Real>& mRhs;
		Grid<Real>& mResidual;
		Grid<Real>& mSearch;
		FlagGrid& mFlags;
		Grid<Real>& mTmp;

		Grid<Real> *mpA0, *mpAi, *mpAj, *mpAk;

		PreconditionType mPcMethod;
		//! preconditioning grids
		Grid<Real> *mpPCA0, *mpPCAi, *mpPCAj, *mpPCAk;

		//! sigma / residual
		Real mSigma;
		//! accuracy of solver (max. residuum)
		Real mAccuracy;
		//! norm of the residual
		Real mResNorm;
}; // GridCg


//! Kernel: Apply symmetric stored Matrix



 struct ApplyMatrix : public KernelBase { ApplyMatrix(FlagGrid& flags, Grid<Real>& dst, Grid<Real>& src, Grid<Real>& A0, Grid<Real>& Ai, Grid<Real>& Aj, Grid<Real>& Ak) :  KernelBase(&flags,0) ,flags(flags),dst(dst),src(src),A0(A0),Ai(Ai),Aj(Aj),Ak(Ak)   { runMessage(); run(); }   inline void op(IndexInt idx, FlagGrid& flags, Grid<Real>& dst, Grid<Real>& src, Grid<Real>& A0, Grid<Real>& Ai, Grid<Real>& Aj, Grid<Real>& Ak )  {
	if (flags[idx] & FlagGrid::TypeZeroPressure) { 
		dst[idx]=0.; return; 
	} 
	if (!flags.isFluid(idx)) {
		dst[idx] = src[idx]; return;
	}    

	dst[idx] =  src[idx] * A0[idx]
				+ src[idx-X] * Ai[idx-X]
				+ src[idx+X] * Ai[idx]
				+ src[idx-Y] * Aj[idx-Y]
				+ src[idx+Y] * Aj[idx]
				+ src[idx-Z] * Ak[idx-Z] 
				+ src[idx+Z] * Ak[idx];
}    inline FlagGrid& getArg0() { return flags; } typedef FlagGrid type0;inline Grid<Real>& getArg1() { return dst; } typedef Grid<Real> type1;inline Grid<Real>& getArg2() { return src; } typedef Grid<Real> type2;inline Grid<Real>& getArg3() { return A0; } typedef Grid<Real> type3;inline Grid<Real>& getArg4() { return Ai; } typedef Grid<Real> type4;inline Grid<Real>& getArg5() { return Aj; } typedef Grid<Real> type5;inline Grid<Real>& getArg6() { return Ak; } typedef Grid<Real> type6; void runMessage() { debMsg("Executing kernel ApplyMatrix ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (IndexInt i = 0; i < _sz; i++) op(i,flags,dst,src,A0,Ai,Aj,Ak);  }   } FlagGrid& flags; Grid<Real>& dst; Grid<Real>& src; Grid<Real>& A0; Grid<Real>& Ai; Grid<Real>& Aj; Grid<Real>& Ak;   };
#line 114 "conjugategrad.h"



//! Kernel: Apply symmetric stored Matrix. 2D version



 struct ApplyMatrix2D : public KernelBase { ApplyMatrix2D(FlagGrid& flags, Grid<Real>& dst, Grid<Real>& src, Grid<Real>& A0, Grid<Real>& Ai, Grid<Real>& Aj, Grid<Real>& Ak) :  KernelBase(&flags,0) ,flags(flags),dst(dst),src(src),A0(A0),Ai(Ai),Aj(Aj),Ak(Ak)   { runMessage(); run(); }   inline void op(IndexInt idx, FlagGrid& flags, Grid<Real>& dst, Grid<Real>& src, Grid<Real>& A0, Grid<Real>& Ai, Grid<Real>& Aj, Grid<Real>& Ak )  {
	unusedParameter(Ak); // only there for parameter compatibility with ApplyMatrix
	
	if (flags[idx] & FlagGrid::TypeZeroPressure) { 
		dst[idx]=0.; return; 
	}

	if (!flags.isFluid(idx)) {
		dst[idx] = src[idx]; return;
	}    

	dst[idx] =  src[idx] * A0[idx]
				+ src[idx-X] * Ai[idx-X]
				+ src[idx+X] * Ai[idx]
				+ src[idx-Y] * Aj[idx-Y]
				+ src[idx+Y] * Aj[idx];
}    inline FlagGrid& getArg0() { return flags; } typedef FlagGrid type0;inline Grid<Real>& getArg1() { return dst; } typedef Grid<Real> type1;inline Grid<Real>& getArg2() { return src; } typedef Grid<Real> type2;inline Grid<Real>& getArg3() { return A0; } typedef Grid<Real> type3;inline Grid<Real>& getArg4() { return Ai; } typedef Grid<Real> type4;inline Grid<Real>& getArg5() { return Aj; } typedef Grid<Real> type5;inline Grid<Real>& getArg6() { return Ak; } typedef Grid<Real> type6; void runMessage() { debMsg("Executing kernel ApplyMatrix2D ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (IndexInt i = 0; i < _sz; i++) op(i,flags,dst,src,A0,Ai,Aj,Ak);  }   } FlagGrid& flags; Grid<Real>& dst; Grid<Real>& src; Grid<Real>& A0; Grid<Real>& Ai; Grid<Real>& Aj; Grid<Real>& Ak;   };
#line 135 "conjugategrad.h"



//! Kernel: Construct the matrix for the poisson equation

 struct MakeLaplaceMatrix : public KernelBase { MakeLaplaceMatrix(FlagGrid& flags, Grid<Real>& A0, Grid<Real>& Ai, Grid<Real>& Aj, Grid<Real>& Ak, MACGrid* fractions = 0) :  KernelBase(&flags,1) ,flags(flags),A0(A0),Ai(Ai),Aj(Aj),Ak(Ak),fractions(fractions)   { runMessage(); run(); }  inline void op(int i, int j, int k, FlagGrid& flags, Grid<Real>& A0, Grid<Real>& Ai, Grid<Real>& Aj, Grid<Real>& Ak, MACGrid* fractions = 0 )  {
	if (!flags.isFluid(i,j,k))
		return;
	
	if(!fractions) {
		// diagonal, A0
		if (!flags.isObstacle(i-1,j,k)) A0(i,j,k) += 1.;
		if (!flags.isObstacle(i+1,j,k)) A0(i,j,k) += 1.;
		if (!flags.isObstacle(i,j-1,k)) A0(i,j,k) += 1.;
		if (!flags.isObstacle(i,j+1,k)) A0(i,j,k) += 1.;
		if (flags.is3D() && !flags.isObstacle(i,j,k-1)) A0(i,j,k) += 1.;
		if (flags.is3D() && !flags.isObstacle(i,j,k+1)) A0(i,j,k) += 1.;
		
		// off-diagonal entries
		if (flags.isFluid(i+1,j,k)) Ai(i,j,k) = -1.;
		if (flags.isFluid(i,j+1,k)) Aj(i,j,k) = -1.;
		if (flags.is3D() && flags.isFluid(i,j,k+1)) Ak(i,j,k) = -1.;
	} else {
		// diagonal
		A0(i,j,k) += fractions->get(i,j,k).x;
		A0(i,j,k) += fractions->get(i+1,j,k).x;
		A0(i,j,k) += fractions->get(i,j,k).y;
		A0(i,j,k) += fractions->get(i,j+1,k).y;
		if (flags.is3D()) A0(i,j,k) += fractions->get(i,j,k).z;
		if (flags.is3D()) A0(i,j,k) += fractions->get(i,j,k+1).z;

		// off-diagonal entries
		Ai(i,j,k) = -fractions->get(i+1,j,k).x;
		Aj(i,j,k) = -fractions->get(i,j+1,k).y;
		if (flags.is3D()) Ak(i,j,k) = -fractions->get(i,j,k+1).z;
	}

}   inline FlagGrid& getArg0() { return flags; } typedef FlagGrid type0;inline Grid<Real>& getArg1() { return A0; } typedef Grid<Real> type1;inline Grid<Real>& getArg2() { return Ai; } typedef Grid<Real> type2;inline Grid<Real>& getArg3() { return Aj; } typedef Grid<Real> type3;inline Grid<Real>& getArg4() { return Ak; } typedef Grid<Real> type4;inline MACGrid* getArg5() { return fractions; } typedef MACGrid type5; void runMessage() { debMsg("Executing kernel MakeLaplaceMatrix ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int k=minZ; k < maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,flags,A0,Ai,Aj,Ak,fractions);  } } else { const int k=0; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,flags,A0,Ai,Aj,Ak,fractions);  } }  } FlagGrid& flags; Grid<Real>& A0; Grid<Real>& Ai; Grid<Real>& Aj; Grid<Real>& Ak; MACGrid* fractions;   };
#line 155 "conjugategrad.h"






} // namespace

#endif 


