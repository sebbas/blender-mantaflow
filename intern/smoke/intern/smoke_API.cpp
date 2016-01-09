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
 * The Original Code is Copyright (C) 2009 by Daniel Genrich
 * All rights reserved.
 *
 * Contributor(s): Daniel Genrich
 *                 Blender Foundation
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file smoke/intern/smoke_API.cpp
 *  \ingroup smoke
 */

#include "FLUID_3D.h"
#include "WTURBULENCE.h"
#include "spectrum.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MANTA.h"
#include "vectorbase.h" /* mantaflow */
#include "smoke_API.h"  /* to ensure valid prototypes */

extern "C" int *smoke_get_manta_flags(struct FLUID_3D *fluid){
	return fluid->_manta_flags;
}

#ifndef WITH_MANTA
extern "C" FLUID_3D *smoke_init(int *res, float dx, float dtdef, int use_heat, int use_fire, int use_colors )
{
	FLUID_3D *fluid = new FLUID_3D(res, dx, dtdef, use_heat, use_fire, use_colors);
	return fluid;
}

extern "C" WTURBULENCE *smoke_turbulence_init(int *res, int amplify, int noisetype, const char *noisefile_path, int use_fire, int use_colors)
{
	if (amplify)
		return new WTURBULENCE(res[0],res[1],res[2], amplify, noisetype, noisefile_path, use_fire, use_colors);
	else 
		return NULL;
}
//////////////////////////////////////////////////////////////////////
#else /*USING MANTAFLOW STRUCTURES*/
//////////////////////////////////////////////////////////////////////
extern "C" FLUID_3D *smoke_init(int *res, float dx, float dtdef, int use_heat, int use_fire, int use_colors, struct SmokeModifierData *smd )
{
	FLUID_3D *fluid = new FLUID_3D(res, dx, dtdef, use_heat, use_fire, use_colors, smd);
	return fluid;
}

extern "C" WTURBULENCE *smoke_turbulence_init(int *res, int amplify, int noisetype, const char *noisefile_path, int use_fire, int use_colors, struct SmokeDomainSettings *sds)
{
	if (amplify)
		return new WTURBULENCE(res[0],res[1],res[2], amplify, noisetype, noisefile_path, use_fire, use_colors, sds);
	else 
		return NULL;
}
#endif /* WITH_MANTA */

extern "C" void smoke_free(FLUID_3D *fluid)
{
	delete fluid;
	fluid = NULL;
}

extern "C" void smoke_turbulence_free(WTURBULENCE *wt)
{
	 delete wt;
	 wt = NULL;
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
	Manta_API::manta_export_script(smd);
	Manta_API::manta_export_grids(smd);
}

extern "C" void smoke_step(FLUID_3D *fluid, float gravity[3], float dtSubdiv)
{
	if (fluid->_fuel) {
		fluid->processBurn(fluid->_fuel, fluid->_density, fluid->_react, fluid->_heat,
						   fluid->_color_r, fluid->_color_g, fluid->_color_b, fluid->_totalCells, (*fluid->_dtFactor)*dtSubdiv);
	}
	fluid->step(dtSubdiv, gravity);

	if (fluid->_fuel) {
		fluid->updateFlame(fluid->_react, fluid->_flame, fluid->_totalCells);
	}
}

#ifndef WITH_MANTA
extern "C" void smoke_turbulence_step(WTURBULENCE *wt, FLUID_3D *fluid)
{
	if (wt->_fuelBig) {
		fluid->processBurn(wt->_fuelBig, wt->_densityBig, wt->_reactBig, 0,
						   wt->_color_rBig, wt->_color_gBig, wt->_color_bBig, wt->_totalCellsBig, fluid->_dt);
	}
	wt->stepTurbulenceFull(fluid->_dt/fluid->_dx, fluid->_xVelocity, fluid->_yVelocity, fluid->_zVelocity, fluid->_obstacles);

	if (wt->_fuelBig) {
		fluid->updateFlame(wt->_reactBig, wt->_flameBig, wt->_totalCellsBig);
	}
}
//////////////////////////////////////////////////////////////////////
#else /*USING MANTAFLOW STRUCTURES*/
//////////////////////////////////////////////////////////////////////
extern "C" void smoke_turbulence_step(WTURBULENCE *wt, FLUID_3D *fluid)
{
	if (wt->_fuelBig) {
		wt->processBurn(fluid->_burning_rate, fluid->_flame_smoke, fluid->_ignition_temp, fluid->_max_temp, fluid->_dt, fluid->_flame_smoke_color);
	}
	wt->stepTurbulenceFull(fluid->_dt/fluid->_dx, fluid->_xVelocity, fluid->_yVelocity, fluid->_zVelocity, fluid->_obstacles);

	if (wt->_fuelBig) {
		wt->updateFlame();
	}
}
#endif /*WITH MANTA*/

extern "C" void smoke_initBlenderRNA(FLUID_3D *fluid, float *alpha, float *beta, float *dt_factor, float *vorticity, int *border_colli, float *burning_rate,
									 float *flame_smoke, float *flame_smoke_color, float *flame_vorticity, float *flame_ignition_temp, float *flame_max_temp)
{
	fluid->initBlenderRNA(alpha, beta, dt_factor, vorticity, border_colli, burning_rate, flame_smoke, flame_smoke_color, flame_vorticity, flame_ignition_temp, flame_max_temp);
}

extern "C" void smoke_initWaveletBlenderRNA(WTURBULENCE *wt, float *strength)
{
	wt->initBlenderRNA(strength);
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

extern "C" void smoke_dissolve(FLUID_3D *fluid, int speed, int log)
{
	data_dissolve(fluid->_density, fluid->_heat, fluid->_color_r, fluid->_color_g, fluid->_color_b, fluid->_totalCells, speed, log);
}

extern "C" void smoke_dissolve_wavelet(WTURBULENCE *wt, int speed, int log)
{
	data_dissolve(wt->_densityBig, 0, wt->_color_rBig, wt->_color_gBig, wt->_color_bBig, wt->_totalCellsBig, speed, log);
}

#ifndef WITH_MANTA
extern "C" void smoke_export(FLUID_3D *fluid, float *dt, float *dx, float **dens, float **react, float **flame, float **fuel, float **heat, 
							 float **heatold, float **vx, float **vy, float **vz, float **r, float **g, float **b, unsigned char **obstacles)
{
	*dens = fluid->_density;
	if(fuel)
		*fuel = fluid->_fuel;
	if(react)
		*react = fluid->_react;
	if(flame)
		*flame = fluid->_flame;
	if(heat)
		*heat = fluid->_heat;
	if(heatold)
		*heatold = fluid->_heatOld;
	*vx = fluid->_xVelocity;
	*vy = fluid->_yVelocity;
	*vz = fluid->_zVelocity;
	if(r)
		*r = fluid->_color_r;
	if(g)
		*g = fluid->_color_g;
	if(b)
		*b = fluid->_color_b;
	*obstacles = fluid->_obstacles;
	*dt = fluid->_dt;
	*dx = fluid->_dx;
}
//////////////////////////////////////////////////////////////////////
#else /*USING MANTAFLOW STRUCTURES*/
//////////////////////////////////////////////////////////////////////
extern "C" void smoke_export(FLUID_3D *fluid, float *dt, float *dx, float **dens, float **react, float **flame, float **fuel, float **heat, 
							 float **manta_inflow, float **vx, float **vy, float **vz, float **r, float **g, float **b, unsigned char **obstacles)
{
	*dens = fluid->_density;
	if(fuel)
		*fuel = fluid->_fuel;
	if(react)
		*react = fluid->_react;
	if(flame)
		*flame = fluid->_flame;
	if(heat)
		*heat = fluid->_heat;
	if(manta_inflow)
		*manta_inflow = fluid->_manta_inflow;
	*vx = fluid->_xVelocity;
	*vy = fluid->_yVelocity;
	*vz = fluid->_zVelocity;
	if(r)
		*r = fluid->_color_r;
	if(g)
		*g = fluid->_color_g;
	if(b)
		*b = fluid->_color_b;
	*obstacles = fluid->_obstacles;
	*dt = fluid->_dt;
	*dx = fluid->_dx;
}
#endif /* WITH_MANTA */

extern "C" void smoke_turbulence_export(WTURBULENCE *wt, float **dens, float **react, float **flame, float **fuel,
                                        float **r, float **g, float **b , float **tcu, float **tcv, float **tcw)
{
	if (!wt)
		return;

	*dens = wt->_densityBig;
	if(fuel)
		*fuel = wt->_fuelBig;
	if(react)
		*react = wt->_reactBig;
	if(flame)
		*flame = wt->_flameBig;
	if(r)
		*r = wt->_color_rBig;
	if(g)
		*g = wt->_color_gBig;
	if(b)
		*b = wt->_color_bBig;
	*tcu = wt->_tcU;
	*tcv = wt->_tcV;
	*tcw = wt->_tcW;
}

extern "C" float *smoke_get_density(FLUID_3D *fluid)
{
	return fluid->_density;
}

extern "C" float *smoke_get_fuel(FLUID_3D *fluid)
{
	return fluid->_fuel;
}

extern "C" float *smoke_get_react(FLUID_3D *fluid)
{
	return fluid->_react;
}

extern "C" float *smoke_get_heat(FLUID_3D *fluid)
{
	return fluid->_heat;
}

extern "C" float *smoke_get_velocity_x(FLUID_3D *fluid)
{
	return fluid->_xVelocity;
}

extern "C" float *smoke_get_velocity_y(FLUID_3D *fluid)
{
	return fluid->_yVelocity;
}

extern "C" float *smoke_get_velocity_z(FLUID_3D *fluid)
{
	return fluid->_zVelocity;
}

extern "C" float *smoke_get_force_x(FLUID_3D *fluid)
{
	return fluid->_xForce;
}

extern "C" float *smoke_get_force_y(FLUID_3D *fluid)
{
	return fluid->_yForce;
}

extern "C" float *smoke_get_force_z(FLUID_3D *fluid)
{
	return fluid->_zForce;
}

extern "C" float *smoke_get_flame(FLUID_3D *fluid)
{
	return fluid->_flame;
}

extern "C" float *smoke_get_color_r(FLUID_3D *fluid)
{
	return fluid->_color_r;
}

extern "C" float *smoke_get_color_g(FLUID_3D *fluid)
{
	return fluid->_color_g;
}

extern "C" float *smoke_get_color_b(FLUID_3D *fluid)
{
	return fluid->_color_b;
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

extern "C" void smoke_get_rgba(FLUID_3D *fluid, float *data, int sequential)
{
	get_rgba(fluid->_color_r, fluid->_color_g, fluid->_color_b, fluid->_density, fluid->_totalCells, data, sequential);
}

extern "C" void smoke_turbulence_get_rgba(WTURBULENCE *wt, float *data, int sequential)
{
	get_rgba(wt->_color_rBig, wt->_color_gBig, wt->_color_bBig, wt->_densityBig, wt->_totalCellsBig, data, sequential);
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

extern "C" void smoke_get_rgba_from_density(FLUID_3D *fluid, float color[3], float *data, int sequential)
{
	get_rgba_from_density(color, fluid->_density, fluid->_totalCells, data, sequential);
}

extern "C" void smoke_turbulence_get_rgba_from_density(WTURBULENCE *wt, float color[3], float *data, int sequential)
{
	get_rgba_from_density(color, wt->_densityBig, wt->_totalCellsBig, data, sequential);
}

extern "C" float *smoke_turbulence_get_density(WTURBULENCE *wt)
{
	return wt ? wt->getDensityBig() : NULL;
}

extern "C" float *smoke_turbulence_get_fuel(WTURBULENCE *wt)
{
	return wt ? wt->getFuelBig() : NULL;
}

extern "C" float *smoke_turbulence_get_react(WTURBULENCE *wt)
{
	return wt ? wt->_reactBig : NULL;
}

extern "C" float *smoke_turbulence_get_color_r(WTURBULENCE *wt)
{
	return wt ? wt->_color_rBig : NULL;
}

extern "C" float *smoke_turbulence_get_color_g(WTURBULENCE *wt)
{
	return wt ? wt->_color_gBig : NULL;
}

extern "C" float *smoke_turbulence_get_color_b(WTURBULENCE *wt)
{
	return wt ? wt->_color_bBig : NULL;
}

extern "C" float *smoke_turbulence_get_flame(WTURBULENCE *wt)
{
	return wt ? wt->getFlameBig() : NULL;
}

extern "C" void smoke_turbulence_get_res(WTURBULENCE *wt, int *res)
{
	if (wt) {
		Vec3Int r = wt->getResBig();
		res[0] = r[0];
		res[1] = r[1];
		res[2] = r[2];
	}
}

extern "C" int smoke_turbulence_get_cells(WTURBULENCE *wt)
{
	if (wt) {
		Vec3Int r = wt->getResBig();
		return r[0] * r[1] * r[2];
	}
	return 0;
}

extern "C" unsigned char *smoke_get_obstacle(FLUID_3D *fluid)
{
	return fluid->_obstacles;
}

extern "C" void smoke_get_ob_velocity(FLUID_3D *fluid, float **x, float **y, float **z)
{
	*x = fluid->_xVelocityOb;
	*y = fluid->_yVelocityOb;
	*z = fluid->_zVelocityOb;
}

#if 0
extern "C" unsigned char *smoke_get_obstacle_anim(FLUID_3D *fluid)
{
	return fluid->_obstaclesAnim;
}
#endif

extern "C" void smoke_turbulence_set_noise(WTURBULENCE *wt, int type, const char *noisefile_path)
{
	wt->setNoise(type, noisefile_path);
}

extern "C" void flame_get_spectrum(unsigned char *spec, int width, float t1, float t2)
{
	spectrum(t1, t2, width, spec);
}

extern "C" int smoke_has_heat(FLUID_3D *fluid)
{
	return (fluid->_heat) ? 1 : 0;
}

extern "C" int smoke_has_fuel(FLUID_3D *fluid)
{
	return (fluid->_fuel) ? 1 : 0;
}

extern "C" int smoke_has_colors(FLUID_3D *fluid)
{
	return (fluid->_color_r && fluid->_color_g && fluid->_color_b) ? 1 : 0;
}

extern "C" int smoke_turbulence_has_fuel(WTURBULENCE *wt)
{
	return (wt->_fuelBig) ? 1 : 0;
}

extern "C" int smoke_turbulence_has_colors(WTURBULENCE *wt)
{
	return (wt->_color_rBig && wt->_color_gBig && wt->_color_bBig) ? 1 : 0;
}

/* additional field initialization */
extern "C" void smoke_ensure_heat(FLUID_3D *fluid)
{
	if (fluid) {
		fluid->initHeat();
	}
}

extern "C" void smoke_ensure_fire(FLUID_3D *fluid, WTURBULENCE *wt)
{
	if (fluid) {
		fluid->initFire();
	}
	if (wt) {
		wt->initFire();
	}
}

extern "C" void smoke_ensure_colors(FLUID_3D *fluid, WTURBULENCE *wt, float init_r, float init_g, float init_b)
{
	if (fluid) {
		fluid->initColors(init_r, init_g, init_b);
	}
	if (wt) {
		wt->initColors(init_r, init_g, init_b);
	}
}


/*MantaFlow funcs*/
//extern "C" int smoke_mantaflow_read(struct SmokeDomainSettings *sds, char* name, bool with_wavelets)
//{
//	return read_mantaflow_sim(sds, name, with_wavelets);
//}
//
//
//extern "C" void manta_write_effectors(struct FLUID_3D *fluid)
//{
//	int size_x = fluid->_xRes;
//	int size_y = fluid->_yRes;
//	int size_z = fluid->_zRes;
//	
//	float *force_x = smoke_get_force_x(fluid);
//	float *force_y = smoke_get_force_y(fluid);
//	float *force_z = smoke_get_force_z(fluid);
////	export_force_fields(size_x, size_y, size_z, force_x, force_y, force_z);
//	/*accumulate all force fields in one grid*/	
//	Manta::Vec3 * accumulated_force = NULL;
//	long index(0);
//	if (fluid->manta_resoution == 3){
//		accumulated_force = (Manta::Vec3*)calloc(size_x * size_y * size_z , sizeof(Manta::Vec3));
//			for (int z(0); z < size_z; z++){
//				for (int y(0); y < size_y; y++){
//					for (int x(0); x < size_x; x++){
//					index = smoke_get_index(x, size_x, y, size_y, z);
//					accumulated_force[index] = Manta::Vec3(force_x[index], force_y[index], force_z[index]);
//				}	
//			}		
//		}
//	}
//	else if (fluid->manta_resoution == 2) {
//		accumulated_force = (Manta::Vec3*)malloc(size_x * size_z * sizeof(Manta::Vec3));
//		int step(0);
//		for (int x(0); x < size_x; x++){
//				for (int z(0); z < size_z; z++){
//					index = smoke_get_index(x, size_x, size_y/2, size_y, z);
//					accumulated_force[x + z * size_x] = Manta::Vec3(force_x[index], force_z[index], 0.0);
//				}	
//		}
//	}
//	else {
//		cout << "ERROR: Manta solver resoltion is neither 2 nor 3; Cannot write forces"<<endl;
//		return;
//	}
//	
//	bool is2D = (fluid->manta_resoution == 2);
//	Manta_API::addGrid(accumulated_force, "forces", "Vec3", size_x, size_y, size_z, is2D);
//}
//
//extern "C" void manta_write_emitters(struct SmokeFlowSettings *sfs, bool highRes, int min_x, int min_y, int min_z, int max_x, int max_y, int max_z, int d_x, int d_y, int d_z,float *influence, float *vel)
//{
////	manta_update_effectors(s, smd->domain->manta_obj, smd->domain, 0.1f);
//	bool is2D = (sfs->smd->domain->fluid->manta_resoution == 2);
//	if (! highRes)
//		Manta_API::addAdaptiveGrid(influence, "density", "s", "float",
//										   min_x, min_y, min_z, max_x, max_y, max_z, is2D);
//	else 
//		Manta_API::addAdaptiveGrid(influence, "xl_density", "xl", "float", min_x, min_y, min_z, max_x, max_y, max_z, is2D);
//	//	export_em_fields(Manta_API::instance()->_emission_map,sfs->density, min_x,  min_y,  min_z,  max_x,  max_y,  max_z,  d_x,  d_y,  d_z,  influence,  vel);
//}
//
///*deprecated*/
//extern "C" void manta_export_obstacles(float * influence, int x, int y, int z)
//{
//	
//	cout << "!!!!!!!!!!Deprecated method manta_export_obstacles is being used" << endl;
//	if (influence == NULL){
//		cout<< "ERROR: empty influence object when exporting smoke obstacles" << endl;
//		return;
//	}
////	Manta_API::export_obstacles(influence, x, y, z);
//}
//
//extern "C" void smoke_mantaflow_stop_sim(struct Manta_API * fluid)
//{
//	if (fluid == NULL){
//		cout<< "ERROR: empty manta_API object when stopping smoke simulation" << endl;
//		return;
//	}
//	fluid->stop_manta_sim();
//}
//
//extern "C" int cell_index_3D(int index_2d, int sizex,int sizey, int sizez)
//{
//	
//	return int(sizey * 0.5) * sizex + 
//	(index_2d % (sizex)) + 
//	int(index_2d/(sizex)) * sizex * sizey;
//}

extern "C" float *smoke_get_inflow_grid(FLUID_3D *fluid)
{
	return fluid->_manta_inflow;
}

extern "C" float *smoke_get_fuel_inflow(FLUID_3D *fluid)
{
	return fluid->_fuel_inflow;
}
