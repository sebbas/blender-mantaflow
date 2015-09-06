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

KERNEL (bnd=1)
void KnProcessBurn(Grid<Real>& fuel,
				   Grid<Real>& density,
				   Grid<Real>& react,
				   Grid<Real>& heat,
				   Grid<Real>& red,
				   Grid<Real>& green,
				   Grid<Real>& blue,
				   float burningRate,
				   float flameSmoke,
				   float ignitionPoint,
				   float tempMax,
				   float dt,
				   Vec3 flameSmokeColor)
{
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
}

KERNEL (bnd=1)
void KnUpdateFlame(Grid<Real>& react, Grid<Real>& flame)
{
	if (react(i,j,k) > 0.0f)
		flame(i,j,k) = pow(react(i,j,k), 0.5f);
	else
		flame(i,j,k) = 0.0f;
}

PYTHON void processBurn(Grid<Real>& fuel,
						Grid<Real>& density,
						Grid<Real>& react,
						Grid<Real>& heat,
						Grid<Real>& red,
						Grid<Real>& green,
						Grid<Real>& blue)
{
	KnProcessBurn(fuel, density, react, heat, red, green, blue, burningRate,
				  flameSmoke, ignitionPoint, tempMax, dtDefault, flameSmokeColor);
}

PYTHON void updateFlame(Grid<Real>& react, Grid<Real>& flame)
{
	KnUpdateFlame(react, flame);
}

} // namespace