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

/* Fluid functions */
struct FLUID *fluid_init(int *res, struct SmokeModifierData *smd);
void fluid_free(struct FLUID *fluid);
void fluid_ensure_obstacle(struct FLUID *fluid, struct SmokeModifierData *smd);
void fluid_ensure_guiding(struct FLUID *fluid, struct SmokeModifierData *smd);
void fluid_ensure_invelocity(struct FLUID *fluid, struct SmokeModifierData *smd);
void fluid_ensure_outflow(struct FLUID *fluid, struct SmokeModifierData *smd);
int fluid_write_config(struct FLUID *fluid, struct SmokeModifierData *smd, int framenr);
int fluid_write_data(struct FLUID *fluid, struct SmokeModifierData *smd, int framenr);
int fluid_read_config(struct FLUID *fluid, struct SmokeModifierData *smd, int framenr);
int fluid_read_data(struct FLUID *fluid, struct SmokeModifierData *smd, int framenr);
int fluid_read_noise(struct FLUID *fluid, struct SmokeModifierData *smd, int framenr);
int fluid_read_mesh(struct FLUID *fluid, struct SmokeModifierData *smd, int framenr);
int fluid_read_particles(struct FLUID *fluid, struct SmokeModifierData *smd, int framenr);
int fluid_read_guiding(struct FLUID *fluid,
                       struct SmokeModifierData *smd,
                       int framenr,
                       bool sourceDomain);
int fluid_update_liquid_structures(struct FLUID *fluid,
                                   struct SmokeModifierData *smd,
                                   int framenr);
int fluid_update_mesh_structures(struct FLUID *fluid, struct SmokeModifierData *smd, int framenr);
int fluid_update_particle_structures(struct FLUID *fluid,
                                     struct SmokeModifierData *smd,
                                     int framenr);
int fluid_bake_data(struct FLUID *fluid, struct SmokeModifierData *smd, int framenr);
int fluid_bake_noise(struct FLUID *fluid, struct SmokeModifierData *smd, int framenr);
int fluid_bake_mesh(struct FLUID *fluid, struct SmokeModifierData *smd, int framenr);
int fluid_bake_particles(struct FLUID *fluid, struct SmokeModifierData *smd, int framenr);
int fluid_bake_guiding(struct FLUID *fluid, struct SmokeModifierData *smd, int framenr);
void fluid_update_variables(struct FLUID *fluid, struct SmokeModifierData *smd);
int fluid_get_frame(struct FLUID *fluid);
float fluid_get_timestep(struct FLUID *fluid);
void fluid_adapt_timestep(struct FLUID *fluid);
bool fluid_needs_realloc(struct FLUID *fluid, struct SmokeModifierData *smd);

/* Fluid accessors */
size_t fluid_get_index(int x, int max_x, int y, int max_y, int z /*, int max_z */);
size_t fluid_get_index2d(int x, int max_x, int y /*, int max_y, int z, int max_z */);
float *fluid_get_velocity_x(struct FLUID *fluid);
float *fluid_get_velocity_y(struct FLUID *fluid);
float *fluid_get_velocity_z(struct FLUID *fluid);
float *fluid_get_ob_velocity_x(struct FLUID *fluid);
float *fluid_get_ob_velocity_y(struct FLUID *fluid);
float *fluid_get_ob_velocity_z(struct FLUID *fluid);
float *fluid_get_guide_velocity_x(struct FLUID *fluid);
float *fluid_get_guide_velocity_y(struct FLUID *fluid);
float *fluid_get_guide_velocity_z(struct FLUID *fluid);
float *fluid_get_in_velocity_x(struct FLUID *fluid);
float *fluid_get_in_velocity_y(struct FLUID *fluid);
float *fluid_get_in_velocity_z(struct FLUID *fluid);
float *fluid_get_force_x(struct FLUID *fluid);
float *fluid_get_force_y(struct FLUID *fluid);
float *fluid_get_force_z(struct FLUID *fluid);
float *fluid_get_phiguide_in(struct FLUID *fluid);
int *fluid_get_num_obstacle(struct FLUID *fluid);
int *fluid_get_num_guide(struct FLUID *fluid);
int fluid_get_res_x(struct FLUID *fluid);
int fluid_get_res_y(struct FLUID *fluid);
int fluid_get_res_z(struct FLUID *fluid);
float *fluid_get_phi_in(struct FLUID *fluid);
float *fluid_get_phiobs_in(struct FLUID *fluid);
float *fluid_get_phiout_in(struct FLUID *fluid);

/* Smoke functions */
void smoke_manta_export(struct FLUID *smoke, struct SmokeModifierData *smd);
void smoke_dissolve(struct FLUID *smoke, int speed, int log);
void smoke_dissolve_wavelet(struct FLUID *smoke, int speed, int log);
void smoke_export(struct FLUID *smoke,
                  float *dt,
                  float *dx,
                  float **dens,
                  float **react,
                  float **flame,
                  float **fuel,
                  float **heat,
                  float **vx,
                  float **vy,
                  float **vz,
                  float **r,
                  float **g,
                  float **b,
                  int **obstacles,
                  float **shadow,
                  float **phiin);
void smoke_turbulence_export(struct FLUID *smoke,
                             float **dens,
                             float **react,
                             float **flame,
                             float **fuel,
                             float **r,
                             float **g,
                             float **b,
                             float **tcu,
                             float **tcv,
                             float **tcw,
                             float **tcu2,
                             float **tcv2,
                             float **tcw2);
void smoke_get_rgba(struct FLUID *smoke, float *data, int sequential);
void smoke_turbulence_get_rgba(struct FLUID *smoke, float *data, int sequential);
void smoke_get_rgba_from_density(struct FLUID *smoke, float color[3], float *data, int sequential);
void smoke_turbulence_get_rgba_from_density(struct FLUID *smoke,
                                            float color[3],
                                            float *data,
                                            int sequential);
void smoke_ensure_heat(struct FLUID *smoke, struct SmokeModifierData *smd);
void smoke_ensure_fire(struct FLUID *smoke, struct SmokeModifierData *smd);
void smoke_ensure_colors(struct FLUID *smoke, struct SmokeModifierData *smd);

/* Smoke accessors */
float *smoke_get_density(struct FLUID *smoke);
float *smoke_get_fuel(struct FLUID *smoke);
float *smoke_get_react(struct FLUID *smoke);
float *smoke_get_heat(struct FLUID *smoke);
float *smoke_get_flame(struct FLUID *smoke);
float *smoke_get_shadow(struct FLUID *fluid);
float *smoke_get_color_r(struct FLUID *smoke);
float *smoke_get_color_g(struct FLUID *smoke);
float *smoke_get_color_b(struct FLUID *smoke);
int *smoke_get_obstacle(struct FLUID *smoke);
float *smoke_get_density_in(struct FLUID *smoke);
float *smoke_get_heat_in(struct FLUID *smoke);
float *smoke_get_color_r_in(struct FLUID *smoke);
float *smoke_get_color_g_in(struct FLUID *smoke);
float *smoke_get_color_b_in(struct FLUID *smoke);
float *smoke_get_fuel_in(struct FLUID *smoke);
float *smoke_get_react_in(struct FLUID *smoke);
float *smoke_get_emission_in(struct FLUID *smoke);
int smoke_has_heat(struct FLUID *smoke);
int smoke_has_fuel(struct FLUID *smoke);
int smoke_has_colors(struct FLUID *smoke);
float *smoke_turbulence_get_density(struct FLUID *smoke);
float *smoke_turbulence_get_fuel(struct FLUID *smoke);
float *smoke_turbulence_get_react(struct FLUID *smoke);
float *smoke_turbulence_get_color_r(struct FLUID *smoke);
float *smoke_turbulence_get_color_g(struct FLUID *smoke);
float *smoke_turbulence_get_color_b(struct FLUID *smoke);
float *smoke_turbulence_get_flame(struct FLUID *smoke);
int smoke_turbulence_has_fuel(struct FLUID *smoke);
int smoke_turbulence_has_colors(struct FLUID *smoke);
void smoke_turbulence_get_res(struct FLUID *smoke, int *res);
int smoke_turbulence_get_cells(struct FLUID *smoke);

/* Liquid functions */
void liquid_manta_export(struct FLUID *smoke, struct SmokeModifierData *smd);
void liquid_ensure_sndparts(struct FLUID *fluid, struct SmokeModifierData *smd);

/* Liquid accessors */
int liquid_get_particle_res_x(struct FLUID *fluid);
int liquid_get_particle_res_y(struct FLUID *fluid);
int liquid_get_particle_res_z(struct FLUID *fluid);
int liquid_get_mesh_res_x(struct FLUID *fluid);
int liquid_get_mesh_res_y(struct FLUID *fluid);
int liquid_get_mesh_res_z(struct FLUID *fluid);
int liquid_get_particle_upres(struct FLUID *fluid);
int liquid_get_mesh_upres(struct FLUID *fluid);
int liquid_get_num_verts(struct FLUID *liquid);
int liquid_get_num_normals(struct FLUID *liquid);
int liquid_get_num_triangles(struct FLUID *liquid);
float liquid_get_vertex_x_at(struct FLUID *liquid, int i);
float liquid_get_vertex_y_at(struct FLUID *liquid, int i);
float liquid_get_vertex_z_at(struct FLUID *liquid, int i);
float liquid_get_normal_x_at(struct FLUID *liquid, int i);
float liquid_get_normal_y_at(struct FLUID *liquid, int i);
float liquid_get_normal_z_at(struct FLUID *liquid, int i);
int liquid_get_triangle_x_at(struct FLUID *liquid, int i);
int liquid_get_triangle_y_at(struct FLUID *liquid, int i);
int liquid_get_triangle_z_at(struct FLUID *liquid, int i);
float liquid_get_vertvel_x_at(struct FLUID *liquid, int i);
float liquid_get_vertvel_y_at(struct FLUID *liquid, int i);
float liquid_get_vertvel_z_at(struct FLUID *liquid, int i);
int liquid_get_num_flip_particles(struct FLUID *liquid);
int liquid_get_num_snd_particles(struct FLUID *liquid);
int liquid_get_flip_particle_flag_at(struct FLUID *liquid, int i);
int liquid_get_snd_particle_flag_at(struct FLUID *liquid, int i);
float liquid_get_flip_particle_position_x_at(struct FLUID *liquid, int i);
float liquid_get_flip_particle_position_y_at(struct FLUID *liquid, int i);
float liquid_get_flip_particle_position_z_at(struct FLUID *liquid, int i);
float liquid_get_flip_particle_velocity_x_at(struct FLUID *liquid, int i);
float liquid_get_flip_particle_velocity_y_at(struct FLUID *liquid, int i);
float liquid_get_flip_particle_velocity_z_at(struct FLUID *liquid, int i);
float liquid_get_snd_particle_position_x_at(struct FLUID *liquid, int i);
float liquid_get_snd_particle_position_y_at(struct FLUID *liquid, int i);
float liquid_get_snd_particle_position_z_at(struct FLUID *liquid, int i);
float liquid_get_snd_particle_velocity_x_at(struct FLUID *liquid, int i);
float liquid_get_snd_particle_velocity_y_at(struct FLUID *liquid, int i);
float liquid_get_snd_particle_velocity_z_at(struct FLUID *liquid, int i);

#ifdef __cplusplus
}
#endif

#endif /* MANTA_FLUID_API_H_ */
