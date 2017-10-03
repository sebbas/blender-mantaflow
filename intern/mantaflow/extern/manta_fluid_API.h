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
void smoke_manta_export(struct FLUID* smoke, struct SmokeModifierData *smd);
void smoke_step(struct FLUID *smoke, int framenr);
void smoke_dissolve(struct FLUID *smoke, int speed, int log);
void smoke_dissolve_wavelet(struct FLUID *smoke, int speed, int log);
void smoke_export(struct FLUID *smoke, float *dt, float *dx, float **dens, float **react, float **flame, float **fuel, float **heat, float **vx, float **vy, float **vz, float **r, float **g, float **b, int **obstacles);
void liquid_export(struct FLUID *liquid, float **phi, float **pp, float **pvel, float **ppSnd, float **pvelSnd, int **ptypeSnd, int **plifeSnd);
void smoke_turbulence_export(struct FLUID *smoke, float **dens, float **react, float **flame, float **fuel, float **r, float **g, float **b , float **tcu, float **tcv, float **tcw, float **tcu2, float **tcv2, float **tcw2);
float *smoke_get_density(struct FLUID *smoke);
float *smoke_get_fuel(struct FLUID *smoke);
float *smoke_get_react(struct FLUID *smoke);
float *smoke_get_heat(struct FLUID *smoke);
float *smoke_get_velocity_x(struct FLUID *smoke);
float *smoke_get_velocity_y(struct FLUID *smoke);
float *smoke_get_velocity_z(struct FLUID *smoke);
float *smoke_get_ob_velocity_x(struct FLUID *fluid);
float *smoke_get_ob_velocity_y(struct FLUID *fluid);
float *smoke_get_ob_velocity_z(struct FLUID *fluid);
float *smoke_get_in_velocity_x(struct FLUID *fluid);
float *smoke_get_in_velocity_y(struct FLUID *fluid);
float *smoke_get_in_velocity_z(struct FLUID *fluid);
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
int *smoke_get_obstacle(struct FLUID *smoke);
void smoke_get_ob_velocity(struct FLUID *smoke, float **x, float **y, float **z);
int smoke_has_heat(struct FLUID *smoke);
int smoke_has_fuel(struct FLUID *smoke);
int smoke_has_colors(struct FLUID *smoke);
int smoke_turbulence_has_fuel(struct FLUID *smoke);
int smoke_turbulence_has_colors(struct FLUID *smoke);
void smoke_ensure_heat(struct FLUID *smoke, struct SmokeModifierData *smd);
void smoke_ensure_fire(struct FLUID *smoke, struct SmokeModifierData *smd);
void smoke_ensure_colors(struct FLUID *smoke, struct SmokeModifierData *smd);
float *smoke_get_guide_velocity_x(struct FLUID *smoke);
float *smoke_get_guide_velocity_y(struct FLUID *smoke);
float *smoke_get_guide_velocity_z(struct FLUID *smoke);

// Liquid grids
float *liquid_get_phiin(struct FLUID *liquid);
float *liquid_get_phiobsin(struct FLUID *liquid);
float *liquid_get_phioutin(struct FLUID *liquid);
void liquid_ensure_init(struct FLUID *liquid, struct SmokeModifierData *smd);
void liquid_manta_export(struct FLUID* smoke, struct SmokeModifierData *smd);

// Liquid Mantaflow IO
void liquid_save_mesh(struct FLUID *liquid, char *filename);
void liquid_save_mesh_high(struct FLUID *liquid, char *filename);
void liquid_save_particles(struct FLUID *liquid, char *filename);
void liquid_save_particle_velocities(struct FLUID *liquid, char *filename);
void liquid_save_data(struct FLUID *liquid, char *pathname);
void liquid_save_data_high(struct FLUID *liquid, char *pathname);
void liquid_load_data(struct FLUID *liquid, char *pathname);
void liquid_load_data_high(struct FLUID *liquid, char *pathname);
void liquid_update_mesh_data(struct FLUID *liquid, char *filename);

// Liquid mesh
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

// Liquids particles
int liquid_get_num_flip_particles(struct FLUID *liquid);
int liquid_get_num_snd_particles(struct FLUID *liquid);

int liquid_get_flip_particle_flag_at(struct FLUID *liquid, int i);
int liquid_get_snd_particle_flag_at(struct FLUID *liquid, int i);

float liquid_get_flip_particle_position_x_at(struct FLUID *liquid, int i);
float liquid_get_flip_particle_position_y_at(struct FLUID *liquid, int i);
float liquid_get_flip_particle_position_z_at(struct FLUID *liquid, int i);
float liquid_get_snd_particle_position_x_at(struct FLUID *liquid, int i);
float liquid_get_snd_particle_position_y_at(struct FLUID *liquid, int i);
float liquid_get_snd_particle_position_z_at(struct FLUID *liquid, int i);

float liquid_get_flip_particle_velocity_x_at(struct FLUID *liquid, int i);
float liquid_get_flip_particle_velocity_y_at(struct FLUID *liquid, int i);
float liquid_get_flip_particle_velocity_z_at(struct FLUID *liquid, int i);
float liquid_get_snd_particle_velocity_x_at(struct FLUID *liquid, int i);
float liquid_get_snd_particle_velocity_y_at(struct FLUID *liquid, int i);
float liquid_get_snd_particle_velocity_z_at(struct FLUID *liquid, int i);
int liquid_get_snd_particle_type_at(struct FLUID *liquid, int i);

void liquid_set_flip_particle_data(struct FLUID* liquid, float* buffer, int numParts);
void liquid_set_snd_particle_data(struct FLUID* liquid, float* buffer, int numParts);

void liquid_set_flip_particle_velocity(struct FLUID* liquid, float* buffer, int numParts);
void liquid_set_snd_particle_velocity(struct FLUID* liquid, float* buffer, int numParts);
void liquid_set_snd_particle_type(struct FLUID* liquid, int* buffer, int numParts);
void liquid_set_snd_particle_life(struct FLUID* liquid, int* buffer, int numParts);

// Fluids in general
int *fluid_get_num_obstacle(struct FLUID *fluid);
float *fluid_get_inflow(struct FLUID *fluid);
float *fluid_get_phiguidein(struct FLUID *fluid);
int *fluid_get_num_guide(struct FLUID *fluid);
int fluid_get_res_x(struct FLUID *fluid);
int fluid_get_res_y(struct FLUID *fluid);
int fluid_get_res_z(struct FLUID *fluid);
void fluid_ensure_obstacle(struct FLUID *fluid, struct SmokeModifierData *smd);
void fluid_ensure_guiding(struct FLUID *fluid, struct SmokeModifierData *smd);
void fluid_ensure_invelocity(struct FLUID *fluid, struct SmokeModifierData *smd);
void fluid_ensure_sndparts(struct FLUID *fluid, struct SmokeModifierData *smd);

#ifdef __cplusplus
}
#endif

#endif /* MANTA_FLUID_API_H_ */
