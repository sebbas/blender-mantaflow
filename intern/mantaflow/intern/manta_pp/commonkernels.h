




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/Users/sbarschkis/Developer/Mantaflow/blenderIntegration/mantaflowgit/source/commonkernels.h"
/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Common grid kernels
 *
 ******************************************************************************/

#ifndef _COMMONKERNELS_H
#define _COMMONKERNELS_H

#include "general.h"
#include "kernel.h"
#include "grid.h"

namespace Manta {
   
//! Kernel: Invert real values, if positive and fluid


 struct InvertCheckFluid : public KernelBase { InvertCheckFluid(FlagGrid& flags, Grid<Real>& grid) :  KernelBase(&flags,0) ,flags(flags),grid(grid)   { runMessage(); run(); }   inline void op(IndexInt idx, FlagGrid& flags, Grid<Real>& grid )  {
	if (flags.isFluid(idx) && grid[idx] > 0)
		grid[idx] = 1.0 / grid[idx];
}    inline FlagGrid& getArg0() { return flags; } typedef FlagGrid type0;inline Grid<Real>& getArg1() { return grid; } typedef Grid<Real> type1; void runMessage() { debMsg("Executing kernel InvertCheckFluid ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (IndexInt i = 0; i < _sz; i++) op(i,flags,grid);  }   } FlagGrid& flags; Grid<Real>& grid;   };
#line 26 "commonkernels.h"



//! Kernel: Squared sum over grid

 struct GridSumSqr : public KernelBase { GridSumSqr(Grid<Real>& grid) :  KernelBase(&grid,0) ,grid(grid) ,sum(0)  { runMessage(); run(); }   inline void op(IndexInt idx, Grid<Real>& grid ,double& sum)  {
	sum += square((double)grid[idx]);    
}    inline operator double () { return sum; } inline double  & getRet() { return sum; }  inline Grid<Real>& getArg0() { return grid; } typedef Grid<Real> type0; void runMessage() { debMsg("Executing kernel GridSumSqr ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  double sum = 0; 
#pragma omp for nowait 
  for (IndexInt i = 0; i < _sz; i++) op(i,grid,sum); 
#pragma omp critical
{this->sum += sum; } }   } Grid<Real>& grid;  double sum;  };
#line 33 "commonkernels.h"



//! Kernel: rotation operator \nabla x v for centered vector fields

 struct CurlOp : public KernelBase { CurlOp(const Grid<Vec3>& grid, Grid<Vec3>& dst) :  KernelBase(&grid,1) ,grid(grid),dst(dst)   { runMessage(); run(); }  inline void op(int i, int j, int k, const Grid<Vec3>& grid, Grid<Vec3>& dst )  {
	Vec3 v = Vec3(0. , 0. , 
			   0.5*((grid(i+1,j,k).y - grid(i-1,j,k).y) - (grid(i,j+1,k).x - grid(i,j-1,k).x)) );
	if(dst.is3D()) {
		v[0] = 0.5*((grid(i,j+1,k).z - grid(i,j-1,k).z) - (grid(i,j,k+1).y - grid(i,j,k-1).y));
		v[1] = 0.5*((grid(i,j,k+1).x - grid(i,j,k-1).x) - (grid(i+1,j,k).z - grid(i-1,j,k).z));
	}
	dst(i,j,k) = v;
}   inline const Grid<Vec3>& getArg0() { return grid; } typedef Grid<Vec3> type0;inline Grid<Vec3>& getArg1() { return dst; } typedef Grid<Vec3> type1; void runMessage() { debMsg("Executing kernel CurlOp ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int k=minZ; k < maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,grid,dst);  } } else { const int k=0; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,grid,dst);  } }  } const Grid<Vec3>& grid; Grid<Vec3>& dst;   };
#line 39 "commonkernels.h"

;

//! Kernel: divergence operator (from MAC grid)

 struct DivergenceOpMAC : public KernelBase { DivergenceOpMAC(Grid<Real>& div, const MACGrid& grid) :  KernelBase(&div,1) ,div(div),grid(grid)   { runMessage(); run(); }  inline void op(int i, int j, int k, Grid<Real>& div, const MACGrid& grid )  {
	Vec3 del = Vec3(grid(i+1,j,k).x, grid(i,j+1,k).y, 0.) - grid(i,j,k); 
	if(grid.is3D()) del[2] += grid(i,j,k+1).z;
	else            del[2]  = 0.;
	div(i,j,k) = del.x + del.y + del.z;
}   inline Grid<Real>& getArg0() { return div; } typedef Grid<Real> type0;inline const MACGrid& getArg1() { return grid; } typedef MACGrid type1; void runMessage() { debMsg("Executing kernel DivergenceOpMAC ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int k=minZ; k < maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,div,grid);  } } else { const int k=0; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,div,grid);  } }  } Grid<Real>& div; const MACGrid& grid;   };
#line 51 "commonkernels.h"



//! Kernel: gradient operator for MAC grid
 struct GradientOpMAC : public KernelBase { GradientOpMAC(MACGrid& gradient, const Grid<Real>& grid) :  KernelBase(&gradient,1) ,gradient(gradient),grid(grid)   { runMessage(); run(); }  inline void op(int i, int j, int k, MACGrid& gradient, const Grid<Real>& grid )  {
	Vec3 grad = (Vec3(grid(i,j,k)) - Vec3(grid(i-1,j,k), grid(i,j-1,k), 0. ));
	if(grid.is3D()) grad[2] -= grid(i,j,k-1);
	else            grad[2]  = 0.;
	gradient(i,j,k) = grad;
}   inline MACGrid& getArg0() { return gradient; } typedef MACGrid type0;inline const Grid<Real>& getArg1() { return grid; } typedef Grid<Real> type1; void runMessage() { debMsg("Executing kernel GradientOpMAC ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int k=minZ; k < maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,gradient,grid);  } } else { const int k=0; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,gradient,grid);  } }  } MACGrid& gradient; const Grid<Real>& grid;   };
#line 59 "commonkernels.h"



//! Kernel: centered gradient operator 
 struct GradientOp : public KernelBase { GradientOp(Grid<Vec3>& gradient, const Grid<Real>& grid) :  KernelBase(&gradient,1) ,gradient(gradient),grid(grid)   { runMessage(); run(); }  inline void op(int i, int j, int k, Grid<Vec3>& gradient, const Grid<Real>& grid )  {
	Vec3 grad = 0.5 * Vec3(        grid(i+1,j,k)-grid(i-1,j,k), 
								   grid(i,j+1,k)-grid(i,j-1,k), 0.);
	if(grid.is3D()) grad[2]= 0.5*( grid(i,j,k+1)-grid(i,j,k-1) );
	gradient(i,j,k) = grad;
}   inline Grid<Vec3>& getArg0() { return gradient; } typedef Grid<Vec3> type0;inline const Grid<Real>& getArg1() { return grid; } typedef Grid<Real> type1; void runMessage() { debMsg("Executing kernel GradientOp ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int k=minZ; k < maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,gradient,grid);  } } else { const int k=0; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,gradient,grid);  } }  } Grid<Vec3>& gradient; const Grid<Real>& grid;   };
#line 67 "commonkernels.h"



//! Kernel: Laplace operator
 struct LaplaceOp : public KernelBase { LaplaceOp(Grid<Real>& laplace, const Grid<Real>& grid) :  KernelBase(&laplace,1) ,laplace(laplace),grid(grid)   { runMessage(); run(); }  inline void op(int i, int j, int k, Grid<Real>& laplace, const Grid<Real>& grid )  {
	laplace(i,j,k) = -(6.0*grid(i,j,k)-grid(i+1,j,k)-grid(i-1,j,k)-grid(i,j+1,k)-grid(i,j-1,k)-grid(i,j,k+1)-grid(i,j,k-1));
}   inline Grid<Real>& getArg0() { return laplace; } typedef Grid<Real> type0;inline const Grid<Real>& getArg1() { return grid; } typedef Grid<Real> type1; void runMessage() { debMsg("Executing kernel LaplaceOp ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int k=minZ; k < maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,laplace,grid);  } } else { const int k=0; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,laplace,grid);  } }  } Grid<Real>& laplace; const Grid<Real>& grid;   };
#line 75 "commonkernels.h"



//! Kernel: get component at MAC positions
 struct GetShiftedComponent : public KernelBase { GetShiftedComponent(const Grid<Vec3>& grid, Grid<Real>& comp, int dim) :  KernelBase(&grid,1) ,grid(grid),comp(comp),dim(dim)   { runMessage(); run(); }  inline void op(int i, int j, int k, const Grid<Vec3>& grid, Grid<Real>& comp, int dim )  {
	Vec3i ishift(i,j,k);
	ishift[dim]--;
	comp(i,j,k) = 0.5*(grid(i,j,k)[dim] + grid(ishift)[dim]);
}   inline const Grid<Vec3>& getArg0() { return grid; } typedef Grid<Vec3> type0;inline Grid<Real>& getArg1() { return comp; } typedef Grid<Real> type1;inline int& getArg2() { return dim; } typedef int type2; void runMessage() { debMsg("Executing kernel GetShiftedComponent ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int k=minZ; k < maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,grid,comp,dim);  } } else { const int k=0; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,grid,comp,dim);  } }  } const Grid<Vec3>& grid; Grid<Real>& comp; int dim;   };
#line 80 "commonkernels.h"

;

//! Kernel: get component (not shifted)
 struct GetComponent : public KernelBase { GetComponent(const Grid<Vec3>& grid, Grid<Real>& comp, int dim) :  KernelBase(&grid,0) ,grid(grid),comp(comp),dim(dim)   { runMessage(); run(); }   inline void op(IndexInt idx, const Grid<Vec3>& grid, Grid<Real>& comp, int dim )  {
	comp[idx] = grid[idx][dim];
}    inline const Grid<Vec3>& getArg0() { return grid; } typedef Grid<Vec3> type0;inline Grid<Real>& getArg1() { return comp; } typedef Grid<Real> type1;inline int& getArg2() { return dim; } typedef int type2; void runMessage() { debMsg("Executing kernel GetComponent ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (IndexInt i = 0; i < _sz; i++) op(i,grid,comp,dim);  }   } const Grid<Vec3>& grid; Grid<Real>& comp; int dim;   };
#line 87 "commonkernels.h"

;

//! Kernel: get norm of centered grid
 struct GridNorm : public KernelBase { GridNorm(Grid<Real>& n, const Grid<Vec3>& grid) :  KernelBase(&n,0) ,n(n),grid(grid)   { runMessage(); run(); }   inline void op(IndexInt idx, Grid<Real>& n, const Grid<Vec3>& grid )  {
	n[idx] = norm(grid[idx]);
}    inline Grid<Real>& getArg0() { return n; } typedef Grid<Real> type0;inline const Grid<Vec3>& getArg1() { return grid; } typedef Grid<Vec3> type1; void runMessage() { debMsg("Executing kernel GridNorm ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (IndexInt i = 0; i < _sz; i++) op(i,n,grid);  }   } Grid<Real>& n; const Grid<Vec3>& grid;   };
#line 92 "commonkernels.h"

;

//! Kernel: set component (not shifted)
 struct SetComponent : public KernelBase { SetComponent(Grid<Vec3>& grid, const Grid<Real>& comp, int dim) :  KernelBase(&grid,0) ,grid(grid),comp(comp),dim(dim)   { runMessage(); run(); }   inline void op(IndexInt idx, Grid<Vec3>& grid, const Grid<Real>& comp, int dim )  {
	grid[idx][dim] = comp[idx];
}    inline Grid<Vec3>& getArg0() { return grid; } typedef Grid<Vec3> type0;inline const Grid<Real>& getArg1() { return comp; } typedef Grid<Real> type1;inline int& getArg2() { return dim; } typedef int type2; void runMessage() { debMsg("Executing kernel SetComponent ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (IndexInt i = 0; i < _sz; i++) op(i,grid,comp,dim);  }   } Grid<Vec3>& grid; const Grid<Real>& comp; int dim;   };
#line 97 "commonkernels.h"

;

//! Kernel: compute centered velocity field from MAC
 struct GetCentered : public KernelBase { GetCentered(Grid<Vec3>& center, const MACGrid& vel) :  KernelBase(&center,1) ,center(center),vel(vel)   { runMessage(); run(); }  inline void op(int i, int j, int k, Grid<Vec3>& center, const MACGrid& vel )  {
	Vec3 v = 0.5 * ( vel(i,j,k) + Vec3(vel(i+1,j,k).x, vel(i,j+1,k).y, 0. ) );
	if(vel.is3D()) v[2] += 0.5 * vel(i,j,k+1).z;
	else           v[2]  = 0.;
	center(i,j,k) = v;
}   inline Grid<Vec3>& getArg0() { return center; } typedef Grid<Vec3> type0;inline const MACGrid& getArg1() { return vel; } typedef MACGrid type1; void runMessage() { debMsg("Executing kernel GetCentered ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int k=minZ; k < maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,center,vel);  } } else { const int k=0; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,center,vel);  } }  } Grid<Vec3>& center; const MACGrid& vel;   };
#line 102 "commonkernels.h"

;

//! Kernel: compute MAC from centered velocity field
 struct GetMAC : public KernelBase { GetMAC(MACGrid& vel, const Grid<Vec3>& center) :  KernelBase(&vel,1) ,vel(vel),center(center)   { runMessage(); run(); }  inline void op(int i, int j, int k, MACGrid& vel, const Grid<Vec3>& center )  {
	Vec3 v = 0.5*(center(i,j,k) + Vec3(center(i-1,j,k).x, center(i,j-1,k).y, 0. ));
	if(vel.is3D()) v[2] += 0.5 * center(i,j,k-1).z; 
	else           v[2]  = 0.;
	vel(i,j,k) = v;
}   inline MACGrid& getArg0() { return vel; } typedef MACGrid type0;inline const Grid<Vec3>& getArg1() { return center; } typedef Grid<Vec3> type1; void runMessage() { debMsg("Executing kernel GetMAC ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int k=minZ; k < maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,vel,center);  } } else { const int k=0; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,vel,center);  } }  } MACGrid& vel; const Grid<Vec3>& center;   };
#line 110 "commonkernels.h"

;

//! Fill in the domain boundary cells (i,j,k=0/size-1) from the neighboring cells
 struct FillInBoundary : public KernelBase { FillInBoundary(Grid<Vec3>& grid, int g) :  KernelBase(&grid,0) ,grid(grid),g(g)   { runMessage(); run(); }  inline void op(int i, int j, int k, Grid<Vec3>& grid, int g )  {
	if (i==0) grid(i,j,k) = grid(i+1,j,k);
	if (j==0) grid(i,j,k) = grid(i,j+1,k);
	if (k==0) grid(i,j,k) = grid(i,j,k+1);
	if (i==grid.getSizeX()-1) grid(i,j,k) = grid(i-1,j,k);
	if (j==grid.getSizeY()-1) grid(i,j,k) = grid(i,j-1,k);
	if (k==grid.getSizeZ()-1) grid(i,j,k) = grid(i,j,k-1);
}   inline Grid<Vec3>& getArg0() { return grid; } typedef Grid<Vec3> type0;inline int& getArg1() { return g; } typedef int type1; void runMessage() { debMsg("Executing kernel FillInBoundary ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int k=minZ; k < maxZ; k++) for (int j=0; j < _maxY; j++) for (int i=0; i < _maxX; i++) op(i,j,k,grid,g);  } } else { const int k=0; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (int j=0; j < _maxY; j++) for (int i=0; i < _maxX; i++) op(i,j,k,grid,g);  } }  } Grid<Vec3>& grid; int g;   };
#line 118 "commonkernels.h"



} // namespace
#endif

