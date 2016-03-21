/*
 * ***** BEGIN GPL LICENSE BLOCK *****
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version. 
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * The Original Code is Copyright (C) 2016 Blender Foundation.
 * All rights reserved.
 *
 * Contributor(s): Sebastian Barschkis (sebbas)
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file mantaflow/intern/manta_smoke_API.cpp
 *  \ingroup mantaflow
 */

#include <cmath>

#include "MANTA.h"
#include "manta_smoke_API.h"
#include "spectrum.h"

extern "C" int *smoke_get_manta_flags(struct MANTA *manta) {
	return manta->getMantaFlags();
}

extern "C" MANTA *smoke_init(int *res, struct SmokeModifierData *smd)
{
	MANTA *manta = new MANTA(res, smd);
	return manta;
}

extern "C" void smoke_free(MANTA *manta)
{
	delete manta;
	manta = NULL;
}

extern "C" size_t smoke_get_index(int x, int max_x, int y, int max_y, int z /*, int max_z */)
{
	return x + y * max_x + z * max_x*max_y;
}

extern "C" size_t smoke_get_index2d(int x, int max_x, int y /*, int max_y, int z, int max_z */)
{
	return x + y * max_x;
}

extern "C" void smoke_manta_export(MANTA* manta, SmokeModifierData *smd)
{
	if (!manta && !smd) return;
	manta->exportScript(smd);
	manta->exportGrids(smd);
}

extern "C" void smoke_step(MANTA *manta, SmokeModifierData *smd)
{
	manta->step(smd);
	manta->updatePointers(smd);
	if (manta->usingHighRes())
		manta->updatePointersHigh(smd);
}

static void data_dissolve(float *density, float *heat, float *r, float *g, float *b, int total_cells, int speed, int log)
{
	if (log) {
		/* max density/speed = dydx */
		float fac = 1.0f - (1.0f / (float)speed);

		for(size_t i = 0; i < total_cells; i++)
		{
			/* density */
			density[i] *= fac;

			/* heat */
			if (heat) {
				heat[i] *= fac;
			}

			/* color */
			if (r) {
				r[i] *= fac;
				g[i] *= fac;
				b[i] *= fac;
			}
		}
	}
	else // linear falloff
	{
		/* max density/speed = dydx */
		float dydx = 1.0f / (float)speed;

		for(size_t i = 0; i < total_cells; i++)
		{
			float d = density[i];
			/* density */
			density[i] -= dydx;
			if (density[i] < 0.0f)
				density[i] = 0.0f;

			/* heat */
			if (heat) {
				if      (fabs(heat[i]) < dydx) heat[i] = 0.0f;
				else if (heat[i] > 0.0f) heat[i] -= dydx;
				else if (heat[i] < 0.0f) heat[i] += dydx;
			}

			/* color */
			if (r && d) {
				r[i] *= (density[i]/d);
				g[i] *= (density[i]/d);
				b[i] *= (density[i]/d);
			}
		}
	}
}

extern "C" void smoke_dissolve(MANTA *manta, int speed, int log)
{
	data_dissolve(manta->getDensity(), manta->getHeat(), manta->getColorR(), manta->getColorG(), manta->getColorB(), manta->getTotalCells(), speed, log);
}

extern "C" void smoke_dissolve_wavelet(MANTA *manta, int speed, int log)
{
	data_dissolve(manta->getDensityHigh(), 0, manta->getColorRHigh(), manta->getColorGHigh(), manta->getColorBHigh(), manta->getTotalCellsHigh(), speed, log);
}

extern "C" void smoke_export(MANTA *manta, float *dt, float *dx, float **dens, float **react, float **flame, float **fuel, float **heat, 
							 float **manta_inflow, float **vx, float **vy, float **vz, float **r, float **g, float **b, unsigned char **obstacles)
{
	*dens = manta->getDensity();
	if (fuel)
		*fuel = manta->getFuel();
	if (react)
		*react = manta->getReact();
	if (flame)
		*flame = manta->getFlame();
	if( heat)
		*heat = manta->getHeat();
	if (manta_inflow)
		*manta_inflow = manta->getDensityInflow();
	*vx = manta->getVelocityX();
	*vy = manta->getVelocityY();
	*vz = manta->getVelocityZ();
	if (r)
		*r = manta->getColorR();
	if (g)
		*g = manta->getColorG();
	if (b)
		*b = manta->getColorB();
	*obstacles = manta->getObstacles();
	*dt = 1; //dummy value, not here needed for manta
	*dx = 1; //dummy value, not here needed for manta
}

extern "C" void smoke_turbulence_export(MANTA *manta, float **dens, float **react, float **flame, float **fuel,
                                        float **r, float **g, float **b , float **tcu, float **tcv, float **tcw)
{
	if (!manta && !manta->usingHighRes())
		return;

	*dens = manta->getDensityHigh();
	if (fuel)
		*fuel = manta->getFuelHigh();
	if (react)
		*react = manta->getReactHigh();
	if (flame)
		*flame = manta->getFlameHigh();
	if (r)
		*r = manta->getColorRHigh();
	if (g)
		*g = manta->getColorGHigh();
	if (b)
		*b = manta->getColorBHigh();
	*tcu = manta->getTextureU();
	*tcv = manta->getTextureV();
	*tcw = manta->getTextureW();
}

extern "C" float *smoke_get_density(MANTA *manta)
{
	return manta->getDensity();
}

extern "C" float *smoke_get_fuel(MANTA *manta)
{
	return manta->getFuel();
}

extern "C" float *smoke_get_react(MANTA *manta)
{
	return manta->getReact();
}

extern "C" float *smoke_get_heat(MANTA *manta)
{
	return manta->getHeat();
}

extern "C" float *smoke_get_velocity_x(MANTA *manta)
{
	return manta->getVelocityX();
}

extern "C" float *smoke_get_velocity_y(MANTA *manta)
{
	return manta->getVelocityY();
}

extern "C" float *smoke_get_velocity_z(MANTA *manta)
{
	return manta->getVelocityZ();
}

extern "C" float *smoke_get_force_x(MANTA *manta)
{
	return manta->getForceX();
}

extern "C" float *smoke_get_force_y(MANTA *manta)
{
	return manta->getForceY();
}

extern "C" float *smoke_get_force_z(MANTA *manta)
{
	return manta->getForceZ();
}

extern "C" float *smoke_get_flame(MANTA *manta)
{
	return manta->getFlame();
}

extern "C" float *smoke_get_color_r(MANTA *manta)
{
	return manta->getColorR();
}

extern "C" float *smoke_get_color_g(MANTA *manta)
{
	return manta->getColorG();
}

extern "C" float *smoke_get_color_b(MANTA *manta)
{
	return manta->getColorB();
}

static void get_rgba(float *r, float *g, float *b, float *a, int total_cells, float *data, int sequential)
{
	int i;
	int m = 4, i_g = 1, i_b = 2, i_a = 3;
	/* sequential data */
	if (sequential) {
		m = 1;
		i_g *= total_cells;
		i_b *= total_cells;
		i_a *= total_cells;
	}

	for (i=0; i<total_cells; i++) {
		float alpha = a[i];
		if (alpha) {
			data[i*m  ] = r[i];
			data[i*m+i_g] = g[i];
			data[i*m+i_b] = b[i];
		}
		else {
			data[i*m  ] = data[i*m+i_g] = data[i*m+i_b] = 0.0f;
		}
		data[i*m+i_a] = alpha;
	}
}

extern "C" void smoke_get_rgba(MANTA *manta, float *data, int sequential)
{
	get_rgba(manta->getColorR(), manta->getColorG(), manta->getColorB(), manta->getDensity(), manta->getTotalCells(), data, sequential);
}

extern "C" void smoke_turbulence_get_rgba(MANTA *manta, float *data, int sequential)
{
	get_rgba(manta->getColorRHigh(), manta->getColorGHigh(), manta->getColorBHigh(), manta->getDensityHigh(), manta->getTotalCellsHigh(), data, sequential);
}

/* get a single color premultiplied voxel grid */
static void get_rgba_from_density(float color[3], float *a, int total_cells, float *data, int sequential)
{
	int i;
	int m = 4, i_g = 1, i_b = 2, i_a = 3;
	/* sequential data */
	if (sequential) {
		m = 1;
		i_g *= total_cells;
		i_b *= total_cells;
		i_a *= total_cells;
	}

	for (i=0; i<total_cells; i++) {
		float alpha = a[i];
		if (alpha) {
			data[i*m  ] = color[0] * alpha;
			data[i*m+i_g] = color[1] * alpha;
			data[i*m+i_b] = color[2] * alpha;
		}
		else {
			data[i*m  ] = data[i*m+i_g] = data[i*m+i_b] = 0.0f;
		}
		data[i*m+i_a] = alpha;
	}
}

extern "C" void smoke_get_rgba_from_density(MANTA *manta, float color[3], float *data, int sequential)
{
	get_rgba_from_density(color, manta->getDensity(), manta->getTotalCells(), data, sequential);
}

extern "C" void smoke_turbulence_get_rgba_from_density(MANTA *manta, float color[3], float *data, int sequential)
{
	get_rgba_from_density(color, manta->getDensityHigh(), manta->getTotalCellsHigh(), data, sequential);
}

extern "C" float *smoke_turbulence_get_density(MANTA *manta)
{
	return (manta && manta->usingHighRes()) ? manta->getDensityHigh() : NULL;
}

extern "C" float *smoke_turbulence_get_fuel(MANTA *manta)
{
	return (manta && manta->usingHighRes()) ? manta->getFuelHigh() : NULL;
}

extern "C" float *smoke_turbulence_get_react(MANTA *manta)
{
	return (manta && manta->usingHighRes()) ? manta->getReactHigh() : NULL;
}

extern "C" float *smoke_turbulence_get_color_r(MANTA *manta)
{
	return (manta && manta->usingHighRes()) ? manta->getColorRHigh() : NULL;
}

extern "C" float *smoke_turbulence_get_color_g(MANTA *manta)
{
	return (manta && manta->usingHighRes()) ? manta->getColorGHigh() : NULL;
}

extern "C" float *smoke_turbulence_get_color_b(MANTA *manta)
{
	return (manta && manta->usingHighRes()) ? manta->getColorBHigh() : NULL;
}

extern "C" float *smoke_turbulence_get_flame(MANTA *manta)
{
	return (manta && manta->usingHighRes()) ? manta->getFlameHigh() : NULL;
}

extern "C" void smoke_turbulence_get_res(MANTA *manta, int *res)
{
	if (manta && manta->usingHighRes()) {
		res[0] = manta->getResXHigh();
		res[1] = manta->getResYHigh();
		res[2] = manta->getResZHigh();
	}
}

extern "C" int smoke_turbulence_get_cells(MANTA *manta)
{
	int total_cells_high = manta->getResXHigh() * manta->getResYHigh() * manta->getResZHigh();
	return (manta && manta->usingHighRes()) ? total_cells_high : 0;
}

extern "C" unsigned char *smoke_get_obstacle(MANTA *manta)
{
	return manta->getObstacles();
}

extern "C" void smoke_get_ob_velocity(MANTA *manta, float **x, float **y, float **z)
{
	*x = manta->getObVelocityX();
	*y = manta->getObVelocityY();
	*z = manta->getObVelocityZ();
}

#if 0
extern "C" unsigned char *smoke_get_obstacle_anim(MANTA *manta)
{
	return manta->getObstaclesAnim();
}
#endif

extern "C" void flame_get_spectrum(unsigned char *spec, int width, float t1, float t2)
{
	spectrum(t1, t2, width, spec);
}

extern "C" int smoke_has_heat(MANTA *manta)
{
	return (manta->getHeat()) ? 1 : 0;
}

extern "C" int smoke_has_fuel(MANTA *manta)
{
	return (manta->getFuel()) ? 1 : 0;
}

extern "C" int smoke_has_colors(MANTA *manta)
{
	return (manta->getColorR() && manta->getColorG() && manta->getColorB()) ? 1 : 0;
}

extern "C" int smoke_turbulence_has_fuel(MANTA *manta)
{
	return (manta->getFuelHigh()) ? 1 : 0;
}

extern "C" int smoke_turbulence_has_colors(MANTA *manta)
{
	return (manta->getColorRHigh() && manta->getColorGHigh() && manta->getColorBHigh()) ? 1 : 0;
}

/* additional field initialization */
extern "C" void smoke_ensure_heat(MANTA *manta, struct SmokeModifierData *smd)
{
	if (manta) {
		manta->initHeat(smd);
		manta->updatePointers(smd);
	}
}

extern "C" void smoke_ensure_fire(MANTA *manta, struct SmokeModifierData *smd)
{
	if (manta) {
		manta->initFire(smd);
		manta->updatePointers(smd);
	}
	if (manta && manta->usingHighRes()) {
		manta->initFireHigh(smd);
		manta->updatePointersHigh(smd);
	}
}

extern "C" void smoke_ensure_colors(MANTA *manta, struct SmokeModifierData *smd)
{
	if (manta) {
		manta->initColors(smd);
		manta->updatePointers(smd);
	}
	if (manta && manta->usingHighRes()) {
		manta->initColorsHigh(smd);
		manta->updatePointersHigh(smd);
	}
}

extern "C" float *smoke_get_inflow_grid(MANTA *manta)
{
	return manta->getDensityInflow();
}

extern "C" float *smoke_get_fuel_inflow(MANTA *manta)
{
	return manta->getFuelInflow();
}
