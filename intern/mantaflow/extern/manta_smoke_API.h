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

/** \file mantaflow/extern/manta_smoke_API.h
 *  \ingroup mantaflow
 */

#ifndef MANTA_SMOKE_API_H
#define MANTA_SMOKE_API_H

#ifdef __cplusplus
extern "C" {
#endif

struct SMOKE;

int *smoke_get_manta_flags(struct SMOKE *smoke);
struct SMOKE *smoke_init(int *res, struct SmokeModifierData *smd);
void smoke_free(struct SMOKE *smoke);
size_t smoke_get_index(int x, int max_x, int y, int max_y, int z /*, int max_z */);
size_t smoke_get_index2d(int x, int max_x, int y /*, int max_y, int z, int max_z */);
void smoke_manta_export(struct SMOKE* smoke, SmokeModifierData *smd);
void smoke_step(struct SMOKE *smoke, SmokeModifierData *smd);
void smoke_dissolve(struct SMOKE *smoke, int speed, int log);
void smoke_dissolve_wavelet(struct SMOKE *smoke, int speed, int log);
void smoke_export(struct SMOKE *smoke, float *dt, float *dx, float **dens, float **react, float **flame, float **fuel, float **heat, float **smoke_inflow, float **vx, float **vy, float **vz, float **r, float **g, float **b, unsigned char **obstacles);
void smoke_turbulence_export(struct SMOKE *smoke, float **dens, float **react, float **flame, float **fuel, float **r, float **g, float **b , float **tcu, float **tcv, float **tcw, float **tcu2, float **tcv2, float **tcw2);
float *smoke_get_density(struct SMOKE *smoke);
float *smoke_get_fuel(struct SMOKE *smoke);
float *smoke_get_react(struct SMOKE *smoke);
float *smoke_get_heat(struct SMOKE *smoke);
float *smoke_get_velocity_x(struct SMOKE *smoke);
float *smoke_get_velocity_y(struct SMOKE *smoke);
float *smoke_get_velocity_z(struct SMOKE *smoke);
float *smoke_get_force_x(struct SMOKE *smoke);
float *smoke_get_force_y(struct SMOKE *smoke);
float *smoke_get_force_z(struct SMOKE *smoke);
float *smoke_get_flame(struct SMOKE *smoke);
float *smoke_get_color_r(struct SMOKE *smoke);
float *smoke_get_color_g(struct SMOKE *smoke);
float *smoke_get_color_b(struct SMOKE *smoke);
void smoke_get_rgba(struct SMOKE *smoke, float *data, int sequential);
void smoke_turbulence_get_rgba(struct SMOKE *smoke, float *data, int sequential);
void smoke_get_rgba_from_density(struct SMOKE *smoke, float color[3], float *data, int sequential);
void smoke_turbulence_get_rgba_from_density(struct SMOKE *smoke, float color[3], float *data, int sequential);
float *smoke_turbulence_get_density(struct SMOKE *smoke);
float *smoke_turbulence_get_fuel(struct SMOKE *smoke);
float *smoke_turbulence_get_react(struct SMOKE *smoke);
float *smoke_turbulence_get_color_r(struct SMOKE *smoke);
float *smoke_turbulence_get_color_g(struct SMOKE *smoke);
float *smoke_turbulence_get_color_b(struct SMOKE *smoke);
float *smoke_turbulence_get_flame(struct SMOKE *smoke);
void smoke_turbulence_get_res(struct SMOKE *smoke, int *res);
int smoke_turbulence_get_cells(struct SMOKE *smoke);
unsigned char *smoke_get_obstacle(struct SMOKE *smoke);
void smoke_get_ob_velocity(struct SMOKE *smoke, float **x, float **y, float **z);
void flame_get_spectrum(unsigned char *spec, int width, float t1, float t2);
int smoke_has_heat(struct SMOKE *smoke);
int smoke_has_fuel(struct SMOKE *smoke);
int smoke_has_colors(struct SMOKE *smoke);
int smoke_turbulence_has_fuel(struct SMOKE *smoke);
int smoke_turbulence_has_colors(struct SMOKE *smoke);
void smoke_ensure_heat(struct SMOKE *smoke, struct SmokeModifierData *smd);
void smoke_ensure_fire(struct SMOKE *smoke, struct SmokeModifierData *smd);
void smoke_ensure_colors(struct SMOKE *smoke, struct SmokeModifierData *smd);
float *smoke_get_inflow_grid(struct SMOKE *smoke);
float *smoke_get_fuel_inflow(struct SMOKE *smoke);

float *liquid_get_phi(struct SMOKE *liquid);
void liquid_ensure_init(struct SMOKE *smoke, struct SmokeModifierData *smd);


#ifdef __cplusplus
}
#endif

#endif /* MANTA_SMOKE_API_H_ */
