

// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).

#line 1 "/Users/sebbas/Developer/Mantaflow/mantaflowDevelop/mantaflowgit/source/test.cpp"
/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
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

// two simple example kernels

struct reductionTest : public KernelBase {
  reductionTest(const Grid<Real> &v) : KernelBase(&v, 0), v(v), sum(0)
  {
    runMessage();
    run();
  }
  inline void op(IndexInt idx, const Grid<Real> &v, double &sum)
  {
    sum += v[idx];
  }
  inline operator double()
  {
    return sum;
  }
  inline double &getRet()
  {
    return sum;
  }
  inline const Grid<Real> &getArg0()
  {
    return v;
  }
  typedef Grid<Real> type0;
  void runMessage()
  {
    debMsg("Executing kernel reductionTest ", 3);
    debMsg("Kernel range"
               << " x " << maxX << " y " << maxY << " z " << minZ << " - " << maxZ << " ",
           4);
  };
  void run()
  {
    const IndexInt _sz = size;
#pragma omp parallel
    {
      double sum = 0;
#pragma omp for nowait
      for (IndexInt i = 0; i < _sz; i++)
        op(i, v, sum);
#pragma omp critical
      {
        this->sum += sum;
      }
    }
  }
  const Grid<Real> &v;
  double sum;
};
#line 27 "test.cpp"

struct minReduction : public KernelBase {
  minReduction(const Grid<Real> &v) : KernelBase(&v, 0), v(v), sum(0)
  {
    runMessage();
    run();
  }
  inline void op(IndexInt idx, const Grid<Real> &v, double &sum)
  {
    if (sum < v[idx])
      sum = v[idx];
  }
  inline operator double()
  {
    return sum;
  }
  inline double &getRet()
  {
    return sum;
  }
  inline const Grid<Real> &getArg0()
  {
    return v;
  }
  typedef Grid<Real> type0;
  void runMessage()
  {
    debMsg("Executing kernel minReduction ", 3);
    debMsg("Kernel range"
               << " x " << maxX << " y " << maxY << " z " << minZ << " - " << maxZ << " ",
           4);
  };
  void run()
  {
    const IndexInt _sz = size;
#pragma omp parallel
    {
      double sum = 0;
#pragma omp for nowait
      for (IndexInt i = 0; i < _sz; i++)
        op(i, v, sum);
#pragma omp critical
      {
        this->sum = min(sum, this->sum);
      }
    }
  }
  const Grid<Real> &v;
  double sum;
};
#line 33 "test.cpp"

// ... add more test code here if necessary ...

}  // namespace Manta
