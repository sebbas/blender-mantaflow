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


#ifndef MANTA_SMOKE_API_H_
#define MANTA_SMOKE_API_H_

#ifdef __cplusplus
extern "C" {
#endif

struct MANTA;

int *smoke_get_manta_flags(struct MANTA *manta);
struct MANTA *smoke_init(int *res, struct SmokeModifierData *smd);
void smoke_free(struct MANTA *manta);
size_t smoke_get_index(int x, int max_x, int y, int max_y, int z /*, int max_z */);
size_t smoke_get_index2d(int x, int max_x, int y /*, int max_y, int z, int max_z */);
void smoke_manta_export(struct SmokeModifierData *smd);
void smoke_step(struct MANTA *manta, SmokeModifierData *smd);
void smoke_dissolve(struct MANTA *manta, int speed, int log);
void smoke_dissolve_wavelet(struct MANTA *manta, int speed, int log);
void smoke_export(struct MANTA *manta, float *dt, float *dx, float **dens, float **react, float **flame, float **fuel, float **heat, float **manta_inflow, float **vx, float **vy, float **vz, float **r, float **g, float **b, unsigned char **obstacles);
void smoke_turbulence_export(struct MANTA *manta, float **dens, float **react, float **flame, float **fuel, float **r, float **g, float **b , float **tcu, float **tcv, float **tcw);
float *smoke_get_density(struct MANTA *manta);
float *smoke_get_fuel(struct MANTA *manta);
float *smoke_get_react(struct MANTA *manta);
float *smoke_get_heat(struct MANTA *manta);
float *smoke_get_velocity_x(struct MANTA *manta);
float *smoke_get_velocity_y(struct MANTA *manta);
float *smoke_get_velocity_z(struct MANTA *manta);
float *smoke_get_force_x(struct MANTA *manta);
float *smoke_get_force_y(struct MANTA *manta);
float *smoke_get_force_z(struct MANTA *manta);
float *smoke_get_flame(struct MANTA *manta);
float *smoke_get_color_r(struct MANTA *manta);
float *smoke_get_color_g(struct MANTA *manta);
float *smoke_get_color_b(struct MANTA *manta);
void smoke_get_rgba(struct MANTA *manta, float *data, int sequential);
void smoke_turbulence_get_rgba(struct MANTA *manta, float *data, int sequential);
void smoke_get_rgba_from_density(struct MANTA *manta, float color[3], float *data, int sequential);
void smoke_turbulence_get_rgba_from_density(struct MANTA *manta, float color[3], float *data, int sequential);
float *smoke_turbulence_get_density(struct MANTA *manta);
float *smoke_turbulence_get_fuel(struct MANTA *manta);
float *smoke_turbulence_get_react(struct MANTA *manta);
float *smoke_turbulence_get_color_r(struct MANTA *manta);
float *smoke_turbulence_get_color_g(struct MANTA *manta);
float *smoke_turbulence_get_color_b(struct MANTA *manta);
float *smoke_turbulence_get_flame(struct MANTA *manta);
void smoke_turbulence_get_res(struct MANTA *manta, int *res);
int smoke_turbulence_get_cells(struct MANTA *manta);
unsigned char *smoke_get_obstacle(struct MANTA *manta);
void smoke_get_ob_velocity(struct MANTA *manta, float **x, float **y, float **z);
void flame_get_spectrum(unsigned char *spec, int width, float t1, float t2);
int smoke_has_heat(struct MANTA *manta);
int smoke_has_fuel(struct MANTA *manta);
int smoke_has_colors(struct MANTA *manta);
int smoke_turbulence_has_fuel(struct MANTA *manta);
int smoke_turbulence_has_colors(struct MANTA *manta);
void smoke_ensure_heat(struct MANTA *manta, struct SmokeModifierData *smd);
void smoke_ensure_fire(struct MANTA *manta, struct SmokeModifierData *smd);
void smoke_ensure_colors(struct MANTA *manta, struct SmokeModifierData *smd);
float *smoke_get_inflow_grid(struct MANTA *manta);
float *smoke_get_fuel_inflow(struct MANTA *manta);

#ifdef __cplusplus
}
#endif

#endif /* MANTA_SMOKE_API_H_ */
