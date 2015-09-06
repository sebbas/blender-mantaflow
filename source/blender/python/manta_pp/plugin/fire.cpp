




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
 * Fire modeling plugin
 *
 ******************************************************************************/

#include <stdio.h>
#include <float.h>
#include "vectorbase.h"
#include "grid.h"
#include "levelset.h"

using namespace std;

namespace Manta {
	
// default domain values
const Real burningRate = 0.75;
const Real flameSmoke = 1.0;
const Real flameVorticity = 0.5;
const Real ignitionPoint = 1.25;
const Real tempMax = 1.75;
const Vec3 flameSmokeColor = Vec3(0.7f, 0.7f, 0.7f);

// default flow values
const Real flowDensity = 1.0;
const Real flowFuelAmount = 1.0;
const Real flowTemp = 1.0;
const Real flowVolumeDensity = 0.0;
const Real flowSurfaceDistance = 1.5;
const Vec3 flowColor = Vec3(0.7f, 0.7f, 0.7f);
	
// other default values
const Real dtDefault = 0.1;
const Real fps = 24.0;
const int absoluteFlow = 1;
const bool withSmoke = true;
const bool withFire = true;















 struct KnProcessBurn : public KernelBase { KnProcessBurn(Grid<Real>& fuel, Grid<Real>& density, Grid<Real>& react, Grid<Real>& heat, Grid<Real>& red, Grid<Real>& green, Grid<Real>& blue, float burningRate, float flameSmoke, float ignitionPoint, float tempMax, float dt, Vec3 flameSmokeColor) :  KernelBase(&fuel,1) ,fuel(fuel),density(density),react(react),heat(heat),red(red),green(green),blue(blue),burningRate(burningRate),flameSmoke(flameSmoke),ignitionPoint(ignitionPoint),tempMax(tempMax),dt(dt),flameSmokeColor(flameSmokeColor)   { run(); }  inline void op(int i, int j, int k, Grid<Real>& fuel, Grid<Real>& density, Grid<Real>& react, Grid<Real>& heat, Grid<Real>& red, Grid<Real>& green, Grid<Real>& blue, float burningRate, float flameSmoke, float ignitionPoint, float tempMax, float dt, Vec3 flameSmokeColor )  {
	// Save initial values
	float origFuel = fuel(i,j,k);
	float origSmoke = density(i,j,k);
	float smokeEmit = 0.0f;
	float flame = 0.0f;
	
	// Process fuel
	fuel(i,j,k) -= burningRate * dt;
	if (fuel(i,j,k) < 0.0f)
	{
		fuel(i,j,k) = 0.0f;
	}
	
	// Process reaction coordinate
	if (origFuel > __FLT_EPSILON__)
	{
		react(i,j,k) *= fuel(i,j,k) / origFuel;
		flame = pow(react(i,j,k), 0.5f);
	}
	else
	{
		react(i,j,k) = 0.0f;
	}
	
	// Set fluid temperature based on fuel burn rate and "flameSmoke" factor
	smokeEmit = (origFuel < 1.0f) ? (1.0f - origFuel) * 0.5f : 0.0f;
	smokeEmit = (smokeEmit + 0.5f) * (origFuel - fuel(i,j,k)) * 0.1f * flameSmoke;
	density(i,j,k) += smokeEmit;
	clamp(density(i,j,k), 0.0f, 1.0f);
	
	// Set fluid temperature from the flame temperature profile
	if (/*heat(i,j,k) &&*/ flame)
	{
		heat(i,j,k) = (1.0f - flame) * ignitionPoint + flame * tempMax;
	}
	
	// Mix new color
	if (smokeEmit > __FLT_EPSILON__)
	{
		float smokeFactor = density(i,j,k) / (origSmoke + smokeEmit);
		red(i,j,k) = (red(i,j,k) + flameSmokeColor.x * smokeEmit) * smokeFactor;
		green(i,j,k) = (green(i,j,k) + flameSmokeColor.y * smokeEmit) * smokeFactor;
		blue(i,j,k) = (blue(i,j,k) + flameSmokeColor.z * smokeEmit) * smokeFactor;
	}
}   inline Grid<Real>& getArg0() { return fuel; } typedef Grid<Real> type0;inline Grid<Real>& getArg1() { return density; } typedef Grid<Real> type1;inline Grid<Real>& getArg2() { return react; } typedef Grid<Real> type2;inline Grid<Real>& getArg3() { return heat; } typedef Grid<Real> type3;inline Grid<Real>& getArg4() { return red; } typedef Grid<Real> type4;inline Grid<Real>& getArg5() { return green; } typedef Grid<Real> type5;inline Grid<Real>& getArg6() { return blue; } typedef Grid<Real> type6;inline float& getArg7() { return burningRate; } typedef float type7;inline float& getArg8() { return flameSmoke; } typedef float type8;inline float& getArg9() { return ignitionPoint; } typedef float type9;inline float& getArg10() { return tempMax; } typedef float type10;inline float& getArg11() { return dt; } typedef float type11;inline Vec3& getArg12() { return flameSmokeColor; } typedef Vec3 type12; void run() {  const int _maxX = maxX; const int _maxY = maxY; for (int k=minZ; k< maxZ; k++) for (int j=1; j< _maxY; j++) for (int i=1; i< _maxX; i++) op(i,j,k, fuel,density,react,heat,red,green,blue,burningRate,flameSmoke,ignitionPoint,tempMax,dt,flameSmokeColor);  } Grid<Real>& fuel; Grid<Real>& density; Grid<Real>& react; Grid<Real>& heat; Grid<Real>& red; Grid<Real>& green; Grid<Real>& blue; float burningRate; float flameSmoke; float ignitionPoint; float tempMax; float dt; Vec3 flameSmokeColor;   };



 struct KnUpdateFlame : public KernelBase { KnUpdateFlame(Grid<Real>& react, Grid<Real>& flame) :  KernelBase(&react,1) ,react(react),flame(flame)   { run(); }  inline void op(int i, int j, int k, Grid<Real>& react, Grid<Real>& flame )  {
	if (react(i,j,k) > 0.0f)
		flame(i,j,k) = pow(react(i,j,k), 0.5f);
	else
		flame(i,j,k) = 0.0f;
}   inline Grid<Real>& getArg0() { return react; } typedef Grid<Real> type0;inline Grid<Real>& getArg1() { return flame; } typedef Grid<Real> type1; void run() {  const int _maxX = maxX; const int _maxY = maxY; for (int k=minZ; k< maxZ; k++) for (int j=1; j< _maxY; j++) for (int i=1; i< _maxX; i++) op(i,j,k, react,flame);  } Grid<Real>& react; Grid<Real>& flame;   };








void processBurn(Grid<Real>& fuel, Grid<Real>& density, Grid<Real>& react, Grid<Real>& heat, Grid<Real>& red, Grid<Real>& green, Grid<Real>& blue) {
	KnProcessBurn(fuel, density, react, heat, red, green, blue, burningRate,
				  flameSmoke, ignitionPoint, tempMax, dtDefault, flameSmokeColor);
} static PyObject* _W_0 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "processBurn" ); PyObject *_retval = 0; { ArgLocker _lock; Grid<Real>& fuel = *_args.getPtr<Grid<Real> >("fuel",0,&_lock); Grid<Real>& density = *_args.getPtr<Grid<Real> >("density",1,&_lock); Grid<Real>& react = *_args.getPtr<Grid<Real> >("react",2,&_lock); Grid<Real>& heat = *_args.getPtr<Grid<Real> >("heat",3,&_lock); Grid<Real>& red = *_args.getPtr<Grid<Real> >("red",4,&_lock); Grid<Real>& green = *_args.getPtr<Grid<Real> >("green",5,&_lock); Grid<Real>& blue = *_args.getPtr<Grid<Real> >("blue",6,&_lock);   _retval = getPyNone(); processBurn(fuel,density,react,heat,red,green,blue);  _args.check(); } pbFinalizePlugin(parent,"processBurn" ); return _retval; } catch(std::exception& e) { pbSetError("processBurn",e.what()); return 0; } } static const Pb::Register _RP_processBurn ("","processBurn",_W_0); 


void updateFlame(Grid<Real>& react, Grid<Real>& flame) {
	KnUpdateFlame(react, flame);
} static PyObject* _W_1 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "updateFlame" ); PyObject *_retval = 0; { ArgLocker _lock; Grid<Real>& react = *_args.getPtr<Grid<Real> >("react",0,&_lock); Grid<Real>& flame = *_args.getPtr<Grid<Real> >("flame",1,&_lock);   _retval = getPyNone(); updateFlame(react,flame);  _args.check(); } pbFinalizePlugin(parent,"updateFlame" ); return _retval; } catch(std::exception& e) { pbSetError("updateFlame",e.what()); return 0; } } static const Pb::Register _RP_updateFlame ("","updateFlame",_W_1); 

} // namespace

