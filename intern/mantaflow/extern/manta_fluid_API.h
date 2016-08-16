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

#ifndef MANTA_FLUID_API_H
#define MANTA_FLUID_API_H

#ifdef __cplusplus
extern "C" {
#endif

struct FLUID;

struct FLUID *smoke_init(int *res, struct SmokeModifierData *smd);
void smoke_free(struct FLUID *smoke);
size_t smoke_get_index(int x, int max_x, int y, int max_y, int z /*, int max_z */);
size_t smoke_get_index2d(int x, int max_x, int y /*, int max_y, int z, int max_z */);
void smoke_manta_export(struct FLUID* smoke, SmokeModifierData *smd);
void smoke_step(struct FLUID *smoke, SmokeModifierData *smd);
void smoke_dissolve(struct FLUID *smoke, int speed, int log);
void smoke_dissolve_wavelet(struct FLUID *smoke, int speed, int log);
void smoke_export(struct FLUID *smoke, float *dt, float *dx, float **dens, float **react, float **flame, float **fuel, float **heat, float **smoke_inflow, float **vx, float **vy, float **vz, float **r, float **g, float **b, unsigned char **obstacles);
void smoke_turbulence_export(struct FLUID *smoke, float **dens, float **react, float **flame, float **fuel, float **r, float **g, float **b , float **tcu, float **tcv, float **tcw, float **tcu2, float **tcv2, float **tcw2);
float *smoke_get_density(struct FLUID *smoke);
float *smoke_get_fuel(struct FLUID *smoke);
float *smoke_get_react(struct FLUID *smoke);
float *smoke_get_heat(struct FLUID *smoke);
float *smoke_get_velocity_x(struct FLUID *smoke);
float *smoke_get_velocity_y(struct FLUID *smoke);
float *smoke_get_velocity_z(struct FLUID *smoke);
float *smoke_get_force_x(struct FLUID *smoke);
float *smoke_get_force_y(struct FLUID *smoke);
float *smoke_get_force_z(struct FLUID *smoke);
float *smoke_get_flame(struct FLUID *smoke);
float *smoke_get_color_r(struct FLUID *smoke);
float *smoke_get_color_g(struct FLUID *smoke);
float *smoke_get_color_b(struct FLUID *smoke);
void smoke_get_rgba(struct FLUID *smoke, float *data, int sequential);
void smoke_turbulence_get_rgba(struct FLUID *smoke, float *data, int sequential);
void smoke_get_rgba_from_density(struct FLUID *smoke, float color[3], float *data, int sequential);
void smoke_turbulence_get_rgba_from_density(struct FLUID *smoke, float color[3], float *data, int sequential);
float *smoke_turbulence_get_density(struct FLUID *smoke);
float *smoke_turbulence_get_fuel(struct FLUID *smoke);
float *smoke_turbulence_get_react(struct FLUID *smoke);
float *smoke_turbulence_get_color_r(struct FLUID *smoke);
float *smoke_turbulence_get_color_g(struct FLUID *smoke);
float *smoke_turbulence_get_color_b(struct FLUID *smoke);
float *smoke_turbulence_get_flame(struct FLUID *smoke);
void smoke_turbulence_get_res(struct FLUID *smoke, int *res);
int smoke_turbulence_get_cells(struct FLUID *smoke);
unsigned char *smoke_get_obstacle(struct FLUID *smoke);
void smoke_get_ob_velocity(struct FLUID *smoke, float **x, float **y, float **z);
void flame_get_spectrum(unsigned char *spec, int width, float t1, float t2);
int smoke_has_heat(struct FLUID *smoke);
int smoke_has_fuel(struct FLUID *smoke);
int smoke_has_colors(struct FLUID *smoke);
int smoke_turbulence_has_fuel(struct FLUID *smoke);
int smoke_turbulence_has_colors(struct FLUID *smoke);
void smoke_ensure_heat(struct FLUID *smoke, struct SmokeModifierData *smd);
void smoke_ensure_fire(struct FLUID *smoke, struct SmokeModifierData *smd);
void smoke_ensure_colors(struct FLUID *smoke, struct SmokeModifierData *smd);
float *smoke_get_inflow_grid(struct FLUID *smoke);
float *smoke_get_fuel_inflow(struct FLUID *smoke);
int *smoke_get_flags(struct FLUID *smoke);
int *smoke_turbulence_get_flags(struct FLUID *smoke);

float *liquid_get_phi(struct FLUID *liquid);
float *liquid_get_phiinit(struct FLUID *liquid);
float *liquid_turbulence_get_phi(struct FLUID *liquid);
void liquid_ensure_init(struct FLUID *liquid, struct SmokeModifierData *smd);
void liquid_save_mesh(struct FLUID *liquid, char *filename);
void liquid_save_data(struct FLUID *liquid, char *pathname);
void liquid_load_data(struct FLUID *liquid, char *pathname);
int liquid_get_num_verts(struct FLUID *liquid);
int liquid_get_num_normals(struct FLUID *liquid);
int liquid_get_num_triangles(struct FLUID *liquid);
float liquid_get_vertex_x_at(struct FLUID *liquid, int i);
float liquid_get_vertex_y_at(struct FLUID *liquid, int i);
float liquid_get_vertex_z_at(struct FLUID *liquid, int i);
float liquid_get_normal_x_at(struct FLUID *liquid, int i);
float liquid_get_normal_y_at(struct FLUID *liquid, int i);
float liquid_get_normal_z_at(struct FLUID *liquid, int i);
float liquid_get_triangle_x_at(struct FLUID *liquid, int i);
float liquid_get_triangle_y_at(struct FLUID *liquid, int i);
float liquid_get_triangle_z_at(struct FLUID *liquid, int i);
void liquid_update_mesh_data(struct FLUID *liquid, char *filename);

#ifdef __cplusplus
}
#endif

#endif /* MANTA_FLUID_API_H_ */
