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

#include "MANTA.h"
#include "spectrum.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "manta_smoke_API.h"

extern "C" int *smoke_get_manta_flags(struct MANTA *manta) {
	return manta->_manta_flags;
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

extern "C" void smoke_manta_export(SmokeModifierData *smd)
{
	if (!smd) return;
	MANTA *manta = smd->domain->fluid;
	manta->export_script(smd);
	manta->export_grids(smd);
}

extern "C" void smoke_step(MANTA *manta, SmokeModifierData *smd)
{
	manta->step(smd);
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
				if      (abs(heat[i]) < dydx) heat[i] = 0.0f;
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
	data_dissolve(manta->_density, manta->_heat, manta->_color_r, manta->_color_g, manta->_color_b, manta->_totalCells, speed, log);
}

extern "C" void smoke_dissolve_wavelet(MANTA *manta, int speed, int log)
{
	data_dissolve(manta->_densityBig, 0, manta->_color_rBig, manta->_color_gBig, manta->_color_bBig, manta->_totalCellsBig, speed, log);
}

extern "C" void smoke_export(MANTA *manta, float *dt, float *dx, float **dens, float **react, float **flame, float **fuel, float **heat, 
							 float **manta_inflow, float **vx, float **vy, float **vz, float **r, float **g, float **b, unsigned char **obstacles)
{
	*dens = manta->_density;
	if(fuel)
		*fuel = manta->_fuel;
	if(react)
		*react = manta->_react;
	if(flame)
		*flame = manta->_flame;
	if(heat)
		*heat = manta->_heat;
	if(manta_inflow)
		*manta_inflow = manta->_manta_inflow;
	*vx = manta->_xVelocity;
	*vy = manta->_yVelocity;
	*vz = manta->_zVelocity;
	if(r)
		*r = manta->_color_r;
	if(g)
		*g = manta->_color_g;
	if(b)
		*b = manta->_color_b;
	*obstacles = manta->_obstacles;
	*dt = 1; //dummy value, not here needed for manta
	*dx = 1; //dummy value, not here needed for manta
}

extern "C" void smoke_turbulence_export(MANTA *manta, float **dens, float **react, float **flame, float **fuel,
                                        float **r, float **g, float **b , float **tcu, float **tcv, float **tcw)
{
	if (!manta && !manta->_using_highres)
		return;

	*dens = manta->_densityBig;
	if(fuel)
		*fuel = manta->_fuelBig;
	if(react)
		*react = manta->_reactBig;
	if(flame)
		*flame = manta->_flameBig;
	if(r)
		*r = manta->_color_rBig;
	if(g)
		*g = manta->_color_gBig;
	if(b)
		*b = manta->_color_bBig;
	*tcu = manta->_tcU;
	*tcv = manta->_tcV;
	*tcw = manta->_tcW;
}

extern "C" float *smoke_get_density(MANTA *manta)
{
	return manta->_density;
}

extern "C" float *smoke_get_fuel(MANTA *manta)
{
	return manta->_fuel;
}

extern "C" float *smoke_get_react(MANTA *manta)
{
	return manta->_react;
}

extern "C" float *smoke_get_heat(MANTA *manta)
{
	return manta->_heat;
}

extern "C" float *smoke_get_velocity_x(MANTA *manta)
{
	return manta->_xVelocity;
}

extern "C" float *smoke_get_velocity_y(MANTA *manta)
{
	return manta->_yVelocity;
}

extern "C" float *smoke_get_velocity_z(MANTA *manta)
{
	return manta->_zVelocity;
}

extern "C" float *smoke_get_force_x(MANTA *manta)
{
	return manta->_xForce;
}

extern "C" float *smoke_get_force_y(MANTA *manta)
{
	return manta->_yForce;
}

extern "C" float *smoke_get_force_z(MANTA *manta)
{
	return manta->_zForce;
}

extern "C" float *smoke_get_flame(MANTA *manta)
{
	return manta->_flame;
}

extern "C" float *smoke_get_color_r(MANTA *manta)
{
	return manta->_color_r;
}

extern "C" float *smoke_get_color_g(MANTA *manta)
{
	return manta->_color_g;
}

extern "C" float *smoke_get_color_b(MANTA *manta)
{
	return manta->_color_b;
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
	get_rgba(manta->_color_r, manta->_color_g, manta->_color_b, manta->_density, manta->_totalCells, data, sequential);
}

extern "C" void smoke_turbulence_get_rgba(MANTA *manta, float *data, int sequential)
{
	get_rgba(manta->_color_rBig, manta->_color_gBig, manta->_color_bBig, manta->_densityBig, manta->_totalCellsBig, data, sequential);
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
	get_rgba_from_density(color, manta->_density, manta->_totalCells, data, sequential);
}

extern "C" void smoke_turbulence_get_rgba_from_density(MANTA *manta, float color[3], float *data, int sequential)
{
	get_rgba_from_density(color, manta->_densityBig, manta->_totalCellsBig, data, sequential);
}

extern "C" float *smoke_turbulence_get_density(MANTA *manta)
{
	return (manta && manta->_using_highres) ? manta->getDensityBig() : NULL;
}

extern "C" float *smoke_turbulence_get_fuel(MANTA *manta)
{
	return (manta && manta->_using_highres) ? manta->getFuelBig() : NULL;
}

extern "C" float *smoke_turbulence_get_react(MANTA *manta)
{
	return (manta && manta->_using_highres) ? manta->_reactBig : NULL;
}

extern "C" float *smoke_turbulence_get_color_r(MANTA *manta)
{
	return (manta && manta->_using_highres) ? manta->_color_rBig : NULL;
}

extern "C" float *smoke_turbulence_get_color_g(MANTA *manta)
{
	return (manta && manta->_using_highres) ? manta->_color_gBig : NULL;
}

extern "C" float *smoke_turbulence_get_color_b(MANTA *manta)
{
	return (manta && manta->_using_highres) ? manta->_color_bBig : NULL;
}

extern "C" float *smoke_turbulence_get_flame(MANTA *manta)
{
	return (manta && manta->_using_highres) ? manta->getFlameBig() : NULL;
}

extern "C" void smoke_turbulence_get_res(MANTA *manta, int *res)
{
	if (manta && manta->_using_highres) {
		Vec3Int r = manta->getResBig();
		res[0] = r[0];
		res[1] = r[1];
		res[2] = r[2];
	}
}

extern "C" int smoke_turbulence_get_cells(MANTA *manta)
{
	if (manta && manta->_using_highres) {
		Vec3Int r = manta->getResBig();
		return r[0] * r[1] * r[2];
	}
	return 0;
}

extern "C" unsigned char *smoke_get_obstacle(MANTA *manta)
{
	return manta->_obstacles;
}

extern "C" void smoke_get_ob_velocity(MANTA *manta, float **x, float **y, float **z)
{
	*x = manta->_xVelocityOb;
	*y = manta->_yVelocityOb;
	*z = manta->_zVelocityOb;
}

#if 0
extern "C" unsigned char *smoke_get_obstacle_anim(MANTA *manta)
{
	return manta->_obstaclesAnim;
}
#endif

extern "C" void flame_get_spectrum(unsigned char *spec, int width, float t1, float t2)
{
	spectrum(t1, t2, width, spec);
}

extern "C" int smoke_has_heat(MANTA *manta)
{
	return (manta->_heat) ? 1 : 0;
}

extern "C" int smoke_has_fuel(MANTA *manta)
{
	return (manta->_fuel) ? 1 : 0;
}

extern "C" int smoke_has_colors(MANTA *manta)
{
	return (manta->_color_r && manta->_color_g && manta->_color_b) ? 1 : 0;
}

extern "C" int smoke_turbulence_has_fuel(MANTA *manta)
{
	return (manta->_fuelBig) ? 1 : 0;
}

extern "C" int smoke_turbulence_has_colors(MANTA *manta)
{
	return (manta->_color_rBig && manta->_color_gBig && manta->_color_bBig) ? 1 : 0;
}

/* additional field initialization */
extern "C" void smoke_ensure_heat(MANTA *manta, struct SmokeModifierData *smd)
{
	if (manta) {
		manta->initHeat(smd);
		manta->update_pointers(smd);
	}
}

extern "C" void smoke_ensure_fire(MANTA *manta, struct SmokeModifierData *smd)
{
	if (manta) {
		manta->initFire(smd);
		manta->update_pointers(smd);
	}
	if (manta && manta->_using_highres) {
		manta->initFireHigh(smd);
		manta->update_pointers_high(smd);
	}
}

extern "C" void smoke_ensure_colors(MANTA *manta, struct SmokeModifierData *smd)
{
	if (manta) {
		manta->initColors(smd);
		manta->update_pointers(smd);
	}
	if (manta && manta->_using_highres) {
		manta->initColorsHigh(smd);
		manta->update_pointers_high(smd);
	}
}

extern "C" float *smoke_get_inflow_grid(MANTA *manta)
{
	return manta->_manta_inflow;
}

extern "C" float *smoke_get_fuel_inflow(MANTA *manta)
{
	return manta->_fuel_inflow;
}
