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

#include "SMOKE.h"
#include "manta_smoke_API.h"
#include "spectrum.h"

extern "C" int *smoke_get_manta_flags(struct SMOKE *smoke) {
	return smoke->getMantaFlags();
}

extern "C" SMOKE *smoke_init(int *res, struct SmokeModifierData *smd)
{
	SMOKE *smoke = new SMOKE(res, smd);
	return smoke;
}

extern "C" void smoke_free(SMOKE *smoke)
{
	delete smoke;
	smoke = NULL;
}

extern "C" size_t smoke_get_index(int x, int max_x, int y, int max_y, int z /*, int max_z */)
{
	return x + y * max_x + z * max_x*max_y;
}

extern "C" size_t smoke_get_index2d(int x, int max_x, int y /*, int max_y, int z, int max_z */)
{
	return x + y * max_x;
}

extern "C" void smoke_manta_export(SMOKE* smoke, SmokeModifierData *smd)
{
	if (!smoke && !smd) return;
	smoke->exportScript(smd);
	smoke->exportGrids(smd);
}

extern "C" void smoke_step(SMOKE *smoke, SmokeModifierData *smd)
{
	smoke->step(smd);
	smoke->updatePointers(smd);
	if (smoke->usingHighRes())
		smoke->updatePointersHigh(smd);
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

extern "C" void smoke_dissolve(SMOKE *smoke, int speed, int log)
{
	data_dissolve(smoke->getDensity(), smoke->getHeat(), smoke->getColorR(), smoke->getColorG(), smoke->getColorB(), smoke->getTotalCells(), speed, log);
}

extern "C" void smoke_dissolve_wavelet(SMOKE *smoke, int speed, int log)
{
	data_dissolve(smoke->getDensityHigh(), 0, smoke->getColorRHigh(), smoke->getColorGHigh(), smoke->getColorBHigh(), smoke->getTotalCellsHigh(), speed, log);
}

extern "C" void smoke_export(SMOKE *smoke, float *dt, float *dx, float **dens, float **react, float **flame, float **fuel, float **heat, 
							 float **smoke_inflow, float **vx, float **vy, float **vz, float **r, float **g, float **b, unsigned char **obstacles)
{
	*dens = smoke->getDensity();
	if (fuel)
		*fuel = smoke->getFuel();
	if (react)
		*react = smoke->getReact();
	if (flame)
		*flame = smoke->getFlame();
	if( heat)
		*heat = smoke->getHeat();
	if (smoke_inflow)
		*smoke_inflow = smoke->getDensityInflow();
	*vx = smoke->getVelocityX();
	*vy = smoke->getVelocityY();
	*vz = smoke->getVelocityZ();
	if (r)
		*r = smoke->getColorR();
	if (g)
		*g = smoke->getColorG();
	if (b)
		*b = smoke->getColorB();
	*obstacles = smoke->getObstacles();
	*dt = 1; //dummy value, not needed for smoke
	*dx = 1; //dummy value, not needed for smoke
}

extern "C" void smoke_turbulence_export(SMOKE *smoke, float **dens, float **react, float **flame, float **fuel,
                                        float **r, float **g, float **b , float **tcu, float **tcv, float **tcw, float **tcu2, float **tcv2, float **tcw2)
{
	if (!smoke && !smoke->usingHighRes())
		return;

	*dens = smoke->getDensityHigh();
	if (fuel)
		*fuel = smoke->getFuelHigh();
	if (react)
		*react = smoke->getReactHigh();
	if (flame)
		*flame = smoke->getFlameHigh();
	if (r)
		*r = smoke->getColorRHigh();
	if (g)
		*g = smoke->getColorGHigh();
	if (b)
		*b = smoke->getColorBHigh();
	*tcu = smoke->getTextureU();
	*tcv = smoke->getTextureV();
	*tcw = smoke->getTextureW();
	
	*tcu2 = smoke->getTextureU2();
	*tcv2 = smoke->getTextureV2();
	*tcw2 = smoke->getTextureW2();
}

extern "C" float *smoke_get_density(SMOKE *smoke)
{
	return smoke->getDensity();
}

extern "C" float *smoke_get_fuel(SMOKE *smoke)
{
	return smoke->getFuel();
}

extern "C" float *smoke_get_react(SMOKE *smoke)
{
	return smoke->getReact();
}

extern "C" float *smoke_get_heat(SMOKE *smoke)
{
	return smoke->getHeat();
}

extern "C" float *smoke_get_velocity_x(SMOKE *smoke)
{
	return smoke->getVelocityX();
}

extern "C" float *smoke_get_velocity_y(SMOKE *smoke)
{
	return smoke->getVelocityY();
}

extern "C" float *smoke_get_velocity_z(SMOKE *smoke)
{
	return smoke->getVelocityZ();
}

extern "C" float *smoke_get_force_x(SMOKE *smoke)
{
	return smoke->getForceX();
}

extern "C" float *smoke_get_force_y(SMOKE *smoke)
{
	return smoke->getForceY();
}

extern "C" float *smoke_get_force_z(SMOKE *smoke)
{
	return smoke->getForceZ();
}

extern "C" float *smoke_get_flame(SMOKE *smoke)
{
	return smoke->getFlame();
}

extern "C" float *smoke_get_color_r(SMOKE *smoke)
{
	return smoke->getColorR();
}

extern "C" float *smoke_get_color_g(SMOKE *smoke)
{
	return smoke->getColorG();
}

extern "C" float *smoke_get_color_b(SMOKE *smoke)
{
	return smoke->getColorB();
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

extern "C" void smoke_get_rgba(SMOKE *smoke, float *data, int sequential)
{
	get_rgba(smoke->getColorR(), smoke->getColorG(), smoke->getColorB(), smoke->getDensity(), smoke->getTotalCells(), data, sequential);
}

extern "C" void smoke_turbulence_get_rgba(SMOKE *smoke, float *data, int sequential)
{
	get_rgba(smoke->getColorRHigh(), smoke->getColorGHigh(), smoke->getColorBHigh(), smoke->getDensityHigh(), smoke->getTotalCellsHigh(), data, sequential);
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

extern "C" void smoke_get_rgba_from_density(SMOKE *smoke, float color[3], float *data, int sequential)
{
	get_rgba_from_density(color, smoke->getDensity(), smoke->getTotalCells(), data, sequential);
}

extern "C" void smoke_turbulence_get_rgba_from_density(SMOKE *smoke, float color[3], float *data, int sequential)
{
	get_rgba_from_density(color, smoke->getDensityHigh(), smoke->getTotalCellsHigh(), data, sequential);
}

extern "C" float *smoke_turbulence_get_density(SMOKE *smoke)
{
	return (smoke && smoke->usingHighRes()) ? smoke->getDensityHigh() : NULL;
}

extern "C" float *smoke_turbulence_get_fuel(SMOKE *smoke)
{
	return (smoke && smoke->usingHighRes()) ? smoke->getFuelHigh() : NULL;
}

extern "C" float *smoke_turbulence_get_react(SMOKE *smoke)
{
	return (smoke && smoke->usingHighRes()) ? smoke->getReactHigh() : NULL;
}

extern "C" float *smoke_turbulence_get_color_r(SMOKE *smoke)
{
	return (smoke && smoke->usingHighRes()) ? smoke->getColorRHigh() : NULL;
}

extern "C" float *smoke_turbulence_get_color_g(SMOKE *smoke)
{
	return (smoke && smoke->usingHighRes()) ? smoke->getColorGHigh() : NULL;
}

extern "C" float *smoke_turbulence_get_color_b(SMOKE *smoke)
{
	return (smoke && smoke->usingHighRes()) ? smoke->getColorBHigh() : NULL;
}

extern "C" float *smoke_turbulence_get_flame(SMOKE *smoke)
{
	return (smoke && smoke->usingHighRes()) ? smoke->getFlameHigh() : NULL;
}

extern "C" float *liquid_turbulence_get_phi(SMOKE *liquid)
{
	return (liquid && liquid->usingHighRes()) ? liquid->getPhiHigh() : NULL;
}

extern "C" void smoke_turbulence_get_res(SMOKE *smoke, int *res)
{
	if (smoke && smoke->usingHighRes()) {
		res[0] = smoke->getResXHigh();
		res[1] = smoke->getResYHigh();
		res[2] = smoke->getResZHigh();
	}
}

extern "C" int smoke_turbulence_get_cells(SMOKE *smoke)
{
	int total_cells_high = smoke->getResXHigh() * smoke->getResYHigh() * smoke->getResZHigh();
	return (smoke && smoke->usingHighRes()) ? total_cells_high : 0;
}

extern "C" unsigned char *smoke_get_obstacle(SMOKE *smoke)
{
	return smoke->getObstacles();
}

extern "C" void smoke_get_ob_velocity(SMOKE *smoke, float **x, float **y, float **z)
{
	*x = smoke->getObVelocityX();
	*y = smoke->getObVelocityY();
	*z = smoke->getObVelocityZ();
}

#if 0
extern "C" unsigned char *smoke_get_obstacle_anim(SMOKE *smoke)
{
	return smoke->getObstaclesAnim();
}
#endif

extern "C" void flame_get_spectrum(unsigned char *spec, int width, float t1, float t2)
{
	spectrum(t1, t2, width, spec);
}

extern "C" int smoke_has_heat(SMOKE *smoke)
{
	return (smoke->getHeat()) ? 1 : 0;
}

extern "C" int smoke_has_fuel(SMOKE *smoke)
{
	return (smoke->getFuel()) ? 1 : 0;
}

extern "C" int smoke_has_colors(SMOKE *smoke)
{
	return (smoke->getColorR() && smoke->getColorG() && smoke->getColorB()) ? 1 : 0;
}

extern "C" int smoke_turbulence_has_fuel(SMOKE *smoke)
{
	return (smoke->getFuelHigh()) ? 1 : 0;
}

extern "C" int smoke_turbulence_has_colors(SMOKE *smoke)
{
	return (smoke->getColorRHigh() && smoke->getColorGHigh() && smoke->getColorBHigh()) ? 1 : 0;
}

/* additional field initialization */
extern "C" void smoke_ensure_heat(SMOKE *smoke, struct SmokeModifierData *smd)
{
	if (smoke) {
		smoke->initHeat(smd);
		smoke->updatePointers(smd);
	}
}

extern "C" void smoke_ensure_fire(SMOKE *smoke, struct SmokeModifierData *smd)
{
	if (smoke) {
		smoke->initFire(smd);
		smoke->updatePointers(smd);
	}
	if (smoke && smoke->usingHighRes()) {
		smoke->initFireHigh(smd);
		smoke->updatePointersHigh(smd);
	}
}

extern "C" void smoke_ensure_colors(SMOKE *smoke, struct SmokeModifierData *smd)
{
	if (smoke) {
		smoke->initColors(smd);
		smoke->updatePointers(smd);
	}
	if (smoke && smoke->usingHighRes()) {
		smoke->initColorsHigh(smd);
		smoke->updatePointersHigh(smd);
	}
}

extern "C" void liquid_ensure_init(SMOKE *smoke, struct SmokeModifierData *smd)
{
	if (smoke) {
		smoke->initLiquid(smd);
		smoke->updatePointers(smd);
	}
}

extern "C" float *smoke_get_inflow_grid(SMOKE *smoke)
{
	return smoke->getDensityInflow();
}

extern "C" float *smoke_get_fuel_inflow(SMOKE *smoke)
{
	return smoke->getFuelInflow();
}

extern "C" float *liquid_get_phi(SMOKE *liquid)
{
	return liquid->getPhi();
}

extern "C" float *liquid_get_phiinit(SMOKE *liquid)
{
	return liquid->getPhiInit();
}

extern "C" void liquid_save_mesh(SMOKE *liquid, char *filename)
{
	if (liquid) {
		liquid->saveMesh(filename);
	}
}

extern "C" void liquid_save_data(SMOKE *liquid, char *pathname)
{
	if (liquid) {
		liquid->saveLiquidData(pathname);
	}
}

extern "C" void liquid_load_data(SMOKE *liquid, char *pathname)
{
	if (liquid) {
		liquid->loadLiquidData(pathname);
	}
}

extern "C" int liquid_get_num_verts(SMOKE *liquid)
{
	return liquid->getNumVertices();
}

extern "C" int liquid_get_num_normals(SMOKE *liquid)
{
	return liquid->getNumNormals();
}

extern "C" int liquid_get_num_triangles(SMOKE *liquid)
{
	return liquid->getNumTriangles();
}

extern "C" float liquid_get_vertex_x_at(SMOKE *liquid, int i)
{
	return liquid->getVertexXAt(i);
}

extern "C" float liquid_get_vertex_y_at(SMOKE *liquid, int i)
{
	return liquid->getVertexYAt(i);
}

extern "C" float liquid_get_vertex_z_at(SMOKE *liquid, int i)
{
	return liquid->getVertexZAt(i);
}

extern "C" float liquid_get_normal_x_at(SMOKE *liquid, int i)
{
	return liquid->getNormalXAt(i);
}

extern "C" float liquid_get_normal_y_at(SMOKE *liquid, int i)
{
	return liquid->getNormalYAt(i);
}

extern "C" float liquid_get_normal_z_at(SMOKE *liquid, int i)
{
	return liquid->getNormalZAt(i);
}

extern "C" float liquid_get_triangle_x_at(SMOKE *liquid, int i)
{
	return liquid->getTriangleXAt(i);
}

extern "C" float liquid_get_triangle_y_at(SMOKE *liquid, int i)
{
	return liquid->getTriangleYAt(i);
}

extern "C" float liquid_get_triangle_z_at(SMOKE *liquid, int i)
{
	return liquid->getTriangleZAt(i);
}

extern "C" void liquid_update_mesh_data(SMOKE *liquid, char* filename)
{
	liquid->updateMeshData(filename);
}

