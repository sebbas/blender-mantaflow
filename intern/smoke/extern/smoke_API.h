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
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file smoke/extern/smoke_API.h
 *  \ingroup smoke
 */


#ifndef SMOKE_API_H_
#define SMOKE_API_H_

#ifdef __cplusplus
extern "C" {
#endif

	
struct FLUID_3D;

// low res
#ifndef WITH_MANTA
struct FLUID_3D *smoke_init(int *res, float dx, float dtdef, int use_heat, int use_fire, int use_colors);
#else
struct FLUID_3D *smoke_init(int *res, float dx, float dtdef, int use_heat, int use_fire, int use_colors, struct SmokeModifierData *smd);
#endif
void smoke_free(struct FLUID_3D *fluid);

void smoke_initBlenderRNA(struct FLUID_3D *fluid, float *alpha, float *beta, float *dt_factor, float *vorticity, int *border_colli, float *burning_rate,
						  float *flame_smoke, float *flame_smoke_color, float *flame_vorticity, float *flame_ignition_temp, float *flame_max_temp);
void smoke_step(struct FLUID_3D *fluid, float gravity[3], float dtSubdiv);

float *smoke_get_density(struct FLUID_3D *fluid);
float *smoke_get_flame(struct FLUID_3D *fluid);
float *smoke_get_fuel(struct FLUID_3D *fluid);
float *smoke_get_react(struct FLUID_3D *fluid);
float *smoke_get_color_r(struct FLUID_3D *fluid);
float *smoke_get_color_g(struct FLUID_3D *fluid);
float *smoke_get_color_b(struct FLUID_3D *fluid);
void smoke_get_rgba(struct FLUID_3D *fluid, float *data, int sequential);
void smoke_get_rgba_from_density(struct FLUID_3D *fluid, float color[3], float *data, int sequential);
float *smoke_get_heat(struct FLUID_3D *fluid);
float *smoke_get_velocity_x(struct FLUID_3D *fluid);
float *smoke_get_velocity_y(struct FLUID_3D *fluid);
float *smoke_get_velocity_z(struct FLUID_3D *fluid);

/* Moving obstacle velocity provided by blender */
void smoke_get_ob_velocity(struct FLUID_3D *fluid, float **x, float **y, float **z);

float *smoke_get_force_x(struct FLUID_3D *fluid);
float *smoke_get_force_y(struct FLUID_3D *fluid);
float *smoke_get_force_z(struct FLUID_3D *fluid);
	
int *smoke_get_manta_flags(struct FLUID_3D *fluid);
float *smoke_get_inflow_grid(struct FLUID_3D *fluid);
float *smoke_get_fuel_inflow(struct FLUID_3D *fluid);
void smoke_manta_export(struct SmokeModifierData *smd);

unsigned char *smoke_get_obstacle(struct FLUID_3D *fluid);

size_t smoke_get_index(int x, int max_x, int y, int max_y, int z);
size_t smoke_get_index2d(int x, int max_x, int y);
void smoke_dissolve(struct FLUID_3D *fluid, int speed, int log);
// wavelet turbulence functions
#ifndef WITH_MANTA
struct WTURBULENCE *smoke_turbulence_init(int *res, int amplify, int noisetype, const char *noisefile_path, int use_fire, int use_colors);
#else
struct WTURBULENCE *smoke_turbulence_init(int *res, int amplify, int noisetype, const char *noisefile_path, int use_fire, int use_colors, struct SmokeDomainSettings *sds);
#endif

void smoke_turbulence_free(struct WTURBULENCE *wt);
void smoke_turbulence_step(struct WTURBULENCE *wt, struct FLUID_3D *fluid);
float *smoke_turbulence_get_density(struct WTURBULENCE *wt);
float *smoke_turbulence_get_color_r(struct WTURBULENCE *wt);
float *smoke_turbulence_get_color_g(struct WTURBULENCE *wt);
float *smoke_turbulence_get_color_b(struct WTURBULENCE *wt);
void smoke_turbulence_get_rgba(struct WTURBULENCE *wt, float *data, int sequential);
void smoke_turbulence_get_rgba_from_density(struct WTURBULENCE *wt, float color[3], float *data, int sequential);
float *smoke_turbulence_get_flame(struct WTURBULENCE *wt);
float *smoke_turbulence_get_fuel(struct WTURBULENCE *wt);
float *smoke_turbulence_get_react(struct WTURBULENCE *wt);
void smoke_turbulence_get_res(struct WTURBULENCE *wt, int *res);
int smoke_turbulence_get_cells(struct WTURBULENCE *wt);
void smoke_turbulence_set_noise(struct WTURBULENCE *wt, int type, const char *noisefile_path);
void smoke_initWaveletBlenderRNA(struct WTURBULENCE *wt, float *strength);
void smoke_dissolve_wavelet(struct WTURBULENCE *wt, int speed, int log);
/* export */
void smoke_export(struct FLUID_3D *fluid, float *dt, float *dx, float **dens, float **react, float **flame, float **fuel, float **heat, float **heatold,
				  float **vx, float **vy, float **vz, float **r, float **g, float **b, unsigned char **obstacles);
void smoke_turbulence_export(struct WTURBULENCE *wt, float **dens, float **react, float **flame, float **fuel,
							 float **r, float **g, float **b, float **tcu, float **tcv, float **tcw);
/* flame spectrum */
void flame_get_spectrum(unsigned char *spec, int width, float t1, float t2);
/* data fields */
int smoke_has_heat(struct FLUID_3D *fluid);
int smoke_has_fuel(struct FLUID_3D *fluid);
int smoke_has_colors(struct FLUID_3D *fluid);
int smoke_turbulence_has_fuel(struct WTURBULENCE *wt);
int smoke_turbulence_has_colors(struct WTURBULENCE *wt);
void smoke_ensure_heat(struct FLUID_3D *fluid);
void smoke_ensure_fire(struct FLUID_3D *fluid, struct WTURBULENCE *wt);
void smoke_ensure_colors(struct FLUID_3D *fluid, struct WTURBULENCE *wt, float init_r, float init_g, float init_b);


	/*Mantaflow functions*/
int smoke_mantaflow_read(struct SmokeDomainSettings *sds, char* name, bool with_wavelets); //1:success, 0: no file,error
void smoke_mantaflow_stop_sim(struct Manta_API *fluid);
void manta_write_effectors(struct FLUID_3D *fluid);
void manta_update_effectors(struct Scene *scene, struct Object *ob,struct SmokeDomainSettings *sds, float dt);
void manta_write_emitters(struct SmokeFlowSettings *sfs, bool highRes, int min_x, int min_y, int min_z, int max_x, int max_y, int max_z, int d_x, int d_y, int d_z,float *influence, float *vel);	
void manta_export_obstacles(float * influence, int x, int y, int z);
	int cell_index_3D(int index_2d,int sizex,int sizey, int sizez); /*For 2D sims only: transforms 2D manta cell index into 3D Blender index*/
#ifdef __cplusplus
}
#endif

#endif /* SMOKE_API_H_ */
