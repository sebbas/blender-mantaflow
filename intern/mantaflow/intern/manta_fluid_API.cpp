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

#include "MANTA_main.h"
#include "manta_fluid_API.h"

/* Fluid functions */
extern "C" MANTA *manta_init(int *res, struct MantaModifierData *mmd)
{
  return new MANTA(res, mmd);
}
extern "C" void manta_free(MANTA *fluid)
{
  delete fluid;
  fluid = NULL;
}

extern "C" void manta_ensure_obstacle(MANTA *fluid, struct MantaModifierData *mmd)
{
  if (fluid) {
    fluid->initObstacle(mmd);
    fluid->updatePointers();
  }
}
extern "C" void manta_ensure_guiding(MANTA *fluid, struct MantaModifierData *mmd)
{
  if (fluid) {
    fluid->initGuiding(mmd);
    fluid->updatePointers();
  }
}
extern "C" void manta_ensure_invelocity(MANTA *fluid, struct MantaModifierData *mmd)
{
  if (fluid) {
    fluid->initInVelocity(mmd);
    fluid->updatePointers();
  }
}
extern "C" void manta_ensure_outflow(MANTA *fluid, struct MantaModifierData *mmd)
{
  if (fluid) {
    fluid->initOutflow(mmd);
    fluid->updatePointers();
  }
}

extern "C" int manta_write_config(MANTA *fluid, MantaModifierData *mmd, int framenr)
{
  if (!fluid || !mmd)
    return 0;
  return fluid->writeConfiguration(mmd, framenr);
}

extern "C" int manta_write_data(MANTA *fluid, MantaModifierData *mmd, int framenr)
{
  if (!fluid || !mmd)
    return 0;
  return fluid->writeData(mmd, framenr);
}

extern "C" int manta_read_config(MANTA *fluid, MantaModifierData *mmd, int framenr)
{
  if (!fluid || !mmd)
    return 0;
  return fluid->readConfiguration(mmd, framenr);
}

extern "C" int manta_read_data(MANTA *fluid, MantaModifierData *mmd, int framenr)
{
  if (!fluid || !mmd)
    return 0;
  return fluid->readData(mmd, framenr);
}

extern "C" int manta_read_noise(MANTA *fluid, MantaModifierData *mmd, int framenr)
{
  if (!fluid || !mmd)
    return 0;
  return fluid->readNoise(mmd, framenr);
}

extern "C" int manta_read_mesh(MANTA *fluid, MantaModifierData *mmd, int framenr)
{
  if (!fluid || !mmd)
    return 0;
  return fluid->readMesh(mmd, framenr);
}

extern "C" int manta_read_particles(MANTA *fluid, MantaModifierData *mmd, int framenr)
{
  if (!fluid || !mmd)
    return 0;
  return fluid->readParticles(mmd, framenr);
}

extern "C" int manta_read_guiding(MANTA *fluid,
                                  MantaModifierData *mmd,
                                  int framenr,
                                  bool sourceDomain)
{
  if (!fluid || !mmd)
    return 0;
  return fluid->readGuiding(mmd, framenr, sourceDomain);
}

extern "C" int manta_update_liquid_structures(MANTA *fluid, MantaModifierData *mmd, int framenr)
{
  if (!fluid || !mmd)
    return 0;
  return fluid->updateFlipStructures(mmd, framenr);
}

extern "C" int manta_update_mesh_structures(MANTA *fluid, MantaModifierData *mmd, int framenr)
{
  if (!fluid || !mmd)
    return 0;
  return fluid->updateMeshStructures(mmd, framenr);
}

extern "C" int manta_update_particle_structures(MANTA *fluid, MantaModifierData *mmd, int framenr)
{
  if (!fluid || !mmd)
    return 0;
  return fluid->updateParticleStructures(mmd, framenr);
}

extern "C" int manta_bake_data(MANTA *fluid, MantaModifierData *mmd, int framenr)
{
  if (!fluid || !mmd)
    return 0;
  return fluid->bakeData(mmd, framenr);
}

extern "C" int manta_bake_noise(MANTA *fluid, MantaModifierData *mmd, int framenr)
{
  if (!fluid || !mmd)
    return 0;
  return fluid->bakeNoise(mmd, framenr);
}

extern "C" int manta_bake_mesh(MANTA *fluid, MantaModifierData *mmd, int framenr)
{
  if (!fluid || !mmd)
    return 0;
  return fluid->bakeMesh(mmd, framenr);
}

extern "C" int manta_bake_particles(MANTA *fluid, MantaModifierData *mmd, int framenr)
{
  if (!fluid || !mmd)
    return 0;
  return fluid->bakeParticles(mmd, framenr);
}

extern "C" int manta_bake_guiding(MANTA *fluid, MantaModifierData *mmd, int framenr)
{
  if (!fluid || !mmd)
    return 0;
  return fluid->bakeGuiding(mmd, framenr);
}

extern "C" void manta_update_variables(MANTA *fluid, MantaModifierData *mmd)
{
  if (fluid)
    fluid->updateVariables(mmd);
}

extern "C" int manta_get_frame(MANTA *fluid)
{
  if (!fluid)
    return 0;
  return fluid->getFrame();
}

extern "C" float manta_get_timestep(MANTA *fluid)
{
  if (!fluid)
    return 0;
  return fluid->getTimestep();
}

extern "C" void manta_adapt_timestep(MANTA *fluid)
{
  if (fluid)
    fluid->adaptTimestep();
}

extern "C" bool manta_needs_realloc(MANTA *fluid, MantaModifierData *mmd)
{
  if (fluid)
    return fluid->needsRealloc(mmd);
  return false;
}

/* Fluid accessors */
extern "C" size_t manta_get_index(int x, int max_x, int y, int max_y, int z /*, int max_z */)
{
  return x + y * max_x + z * max_x * max_y;
}
extern "C" size_t manta_get_index2d(int x, int max_x, int y /*, int max_y, int z, int max_z */)
{
  return x + y * max_x;
}
extern "C" float *manta_get_velocity_x(MANTA *fluid)
{
  return fluid->getVelocityX();
}
extern "C" float *manta_get_velocity_y(MANTA *fluid)
{
  return fluid->getVelocityY();
}
extern "C" float *manta_get_velocity_z(MANTA *fluid)
{
  return fluid->getVelocityZ();
}

extern "C" float *manta_get_ob_velocity_x(MANTA *fluid)
{
  return fluid->getObVelocityX();
}
extern "C" float *manta_get_ob_velocity_y(MANTA *fluid)
{
  return fluid->getObVelocityY();
}
extern "C" float *manta_get_ob_velocity_z(MANTA *fluid)
{
  return fluid->getObVelocityZ();
}

extern "C" float *manta_get_guide_velocity_x(MANTA *fluid)
{
  return fluid->getGuideVelocityX();
}
extern "C" float *manta_get_guide_velocity_y(MANTA *fluid)
{
  return fluid->getGuideVelocityY();
}
extern "C" float *manta_get_guide_velocity_z(MANTA *fluid)
{
  return fluid->getGuideVelocityZ();
}

extern "C" float *manta_get_in_velocity_x(MANTA *fluid)
{
  return fluid->getInVelocityX();
}
extern "C" float *manta_get_in_velocity_y(MANTA *fluid)
{
  return fluid->getInVelocityY();
}
extern "C" float *manta_get_in_velocity_z(MANTA *fluid)
{
  return fluid->getInVelocityZ();
}

extern "C" float *manta_get_force_x(MANTA *fluid)
{
  return fluid->getForceX();
}
extern "C" float *manta_get_force_y(MANTA *fluid)
{
  return fluid->getForceY();
}
extern "C" float *manta_get_force_z(MANTA *fluid)
{
  return fluid->getForceZ();
}

extern "C" float *manta_get_phiguide_in(MANTA *fluid)
{
  return fluid->getPhiGuideIn();
}

extern "C" int *manta_get_num_obstacle(MANTA *fluid)
{
  return fluid->getNumObstacle();
}
extern "C" int *manta_get_num_guide(MANTA *fluid)
{
  return fluid->getNumGuide();
}

extern "C" int manta_get_res_x(MANTA *fluid)
{
  return fluid->getResX();
}
extern "C" int manta_get_res_y(MANTA *fluid)
{
  return fluid->getResY();
}
extern "C" int manta_get_res_z(MANTA *fluid)
{
  return fluid->getResZ();
}

extern "C" float *manta_get_phi_in(MANTA *fluid)
{
  return fluid->getPhiIn();
}
extern "C" float *manta_get_phiobs_in(MANTA *fluid)
{
  return fluid->getPhiObsIn();
}
extern "C" float *manta_get_phiout_in(MANTA *fluid)
{
  return fluid->getPhiOutIn();
}

/* Smoke functions */
extern "C" void manta_smoke_export_script(MANTA *smoke, MantaModifierData *mmd)
{
  if (!smoke || !mmd)
    return;
  smoke->exportSmokeScript(mmd);
}

extern "C" void manta_smoke_export(MANTA *smoke,
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
                             int **obstacle,
                             float **shadow)
{
  if (dens)
    *dens = smoke->getDensity();
  if (fuel)
    *fuel = smoke->getFuel();
  if (react)
    *react = smoke->getReact();
  if (flame)
    *flame = smoke->getFlame();
  if (heat)
    *heat = smoke->getHeat();
  *vx = smoke->getVelocityX();
  *vy = smoke->getVelocityY();
  *vz = smoke->getVelocityZ();
  if (r)
    *r = smoke->getColorR();
  if (g)
    *g = smoke->getColorG();
  if (b)
    *b = smoke->getColorB();
  *obstacle = smoke->getObstacle();
  *shadow = smoke->getShadow();
  *dt = 1;  //dummy value, not needed for smoke
  *dx = 1;  //dummy value, not needed for smoke
}

extern "C" void manta_smoke_turbulence_export(MANTA *smoke,
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
                                        float **tcw2)
{
  if (!smoke && !(smoke->usingNoise()))
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

static void get_rgba(
    float *r, float *g, float *b, float *a, int total_cells, float *data, int sequential)
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

  for (i = 0; i < total_cells; i++) {
    float alpha = a[i];
    if (alpha) {
      data[i * m] = r[i];
      data[i * m + i_g] = g[i];
      data[i * m + i_b] = b[i];
    }
    else {
      data[i * m] = data[i * m + i_g] = data[i * m + i_b] = 0.0f;
    }
    data[i * m + i_a] = alpha;
  }
}

extern "C" void manta_smoke_get_rgba(MANTA *smoke, float *data, int sequential)
{
  get_rgba(smoke->getColorR(),
           smoke->getColorG(),
           smoke->getColorB(),
           smoke->getDensity(),
           smoke->getTotalCells(),
           data,
           sequential);
}

extern "C" void manta_smoke_turbulence_get_rgba(MANTA *smoke, float *data, int sequential)
{
  get_rgba(smoke->getColorRHigh(),
           smoke->getColorGHigh(),
           smoke->getColorBHigh(),
           smoke->getDensityHigh(),
           smoke->getTotalCellsHigh(),
           data,
           sequential);
}

/* get a single color premultiplied voxel grid */
static void get_rgba_from_density(
    float color[3], float *a, int total_cells, float *data, int sequential)
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

  for (i = 0; i < total_cells; i++) {
    float alpha = a[i];
    if (alpha) {
      data[i * m] = color[0] * alpha;
      data[i * m + i_g] = color[1] * alpha;
      data[i * m + i_b] = color[2] * alpha;
    }
    else {
      data[i * m] = data[i * m + i_g] = data[i * m + i_b] = 0.0f;
    }
    data[i * m + i_a] = alpha;
  }
}

extern "C" void manta_smoke_get_rgba_from_density(MANTA *smoke,
                                            float color[3],
                                            float *data,
                                            int sequential)
{
  get_rgba_from_density(color, smoke->getDensity(), smoke->getTotalCells(), data, sequential);
}

extern "C" void manta_smoke_turbulence_get_rgba_from_density(MANTA *smoke,
                                                       float color[3],
                                                       float *data,
                                                       int sequential)
{
  get_rgba_from_density(
      color, smoke->getDensityHigh(), smoke->getTotalCellsHigh(), data, sequential);
}

extern "C" void manta_smoke_ensure_heat(MANTA *smoke, struct MantaModifierData *mmd)
{
  if (smoke) {
    smoke->initHeat(mmd);
    smoke->updatePointers();
  }
}

extern "C" void manta_smoke_ensure_fire(MANTA *smoke, struct MantaModifierData *mmd)
{
  if (smoke) {
    smoke->initFire(mmd);
    smoke->updatePointers();
  }
  if (smoke && smoke->usingNoise()) {
    smoke->initFireHigh(mmd);
    smoke->updatePointersNoise();
  }
}

extern "C" void manta_smoke_ensure_colors(MANTA *smoke, struct MantaModifierData *mmd)
{
  if (smoke) {
    smoke->initColors(mmd);
    smoke->updatePointers();
  }
  if (smoke && smoke->usingNoise()) {
    smoke->initColorsHigh(mmd);
    smoke->updatePointersNoise();
  }
}

/* Smoke accessors */
extern "C" float *manta_smoke_get_density(MANTA *smoke)
{
  return smoke->getDensity();
}
extern "C" float *manta_smoke_get_fuel(MANTA *smoke)
{
  return smoke->getFuel();
}
extern "C" float *manta_smoke_get_react(MANTA *smoke)
{
  return smoke->getReact();
}
extern "C" float *manta_smoke_get_heat(MANTA *smoke)
{
  return smoke->getHeat();
}
extern "C" float *manta_smoke_get_flame(MANTA *smoke)
{
  return smoke->getFlame();
}
extern "C" float *manta_smoke_get_shadow(MANTA *fluid)
{
  return fluid->getShadow();
}

extern "C" float *manta_smoke_get_color_r(MANTA *smoke)
{
  return smoke->getColorR();
}
extern "C" float *manta_smoke_get_color_g(MANTA *smoke)
{
  return smoke->getColorG();
}
extern "C" float *manta_smoke_get_color_b(MANTA *smoke)
{
  return smoke->getColorB();
}

extern "C" int *manta_smoke_get_obstacle(MANTA *smoke)
{
  return smoke->getObstacle();
}

extern "C" float *manta_smoke_get_density_in(MANTA *smoke)
{
  return smoke->getDensityIn();
}
extern "C" float *manta_smoke_get_heat_in(MANTA *smoke)
{
  return smoke->getHeatIn();
}
extern "C" float *manta_smoke_get_color_r_in(MANTA *smoke)
{
  return smoke->getColorRIn();
}
extern "C" float *manta_smoke_get_color_g_in(MANTA *smoke)
{
  return smoke->getColorGIn();
}
extern "C" float *manta_smoke_get_color_b_in(MANTA *smoke)
{
  return smoke->getColorBIn();
}
extern "C" float *manta_smoke_get_fuel_in(MANTA *smoke)
{
  return smoke->getFuelIn();
}
extern "C" float *manta_smoke_get_react_in(MANTA *smoke)
{
  return smoke->getReactIn();
}
extern "C" float *manta_smoke_get_emission_in(MANTA *smoke)
{
  return smoke->getEmissionIn();
}

extern "C" int manta_smoke_has_heat(MANTA *smoke)
{
  return (smoke->getHeat()) ? 1 : 0;
}
extern "C" int manta_smoke_has_fuel(MANTA *smoke)
{
  return (smoke->getFuel()) ? 1 : 0;
}
extern "C" int manta_smoke_has_colors(MANTA *smoke)
{
  return (smoke->getColorR() && smoke->getColorG() && smoke->getColorB()) ? 1 : 0;
}

extern "C" float *manta_smoke_turbulence_get_density(MANTA *smoke)
{
  return (smoke && smoke->usingNoise()) ? smoke->getDensityHigh() : NULL;
}
extern "C" float *manta_smoke_turbulence_get_fuel(MANTA *smoke)
{
  return (smoke && smoke->usingNoise()) ? smoke->getFuelHigh() : NULL;
}
extern "C" float *manta_smoke_turbulence_get_react(MANTA *smoke)
{
  return (smoke && smoke->usingNoise()) ? smoke->getReactHigh() : NULL;
}
extern "C" float *manta_smoke_turbulence_get_color_r(MANTA *smoke)
{
  return (smoke && smoke->usingNoise()) ? smoke->getColorRHigh() : NULL;
}
extern "C" float *manta_smoke_turbulence_get_color_g(MANTA *smoke)
{
  return (smoke && smoke->usingNoise()) ? smoke->getColorGHigh() : NULL;
}
extern "C" float *manta_smoke_turbulence_get_color_b(MANTA *smoke)
{
  return (smoke && smoke->usingNoise()) ? smoke->getColorBHigh() : NULL;
}
extern "C" float *manta_smoke_turbulence_get_flame(MANTA *smoke)
{
  return (smoke && smoke->usingNoise()) ? smoke->getFlameHigh() : NULL;
}

extern "C" int manta_smoke_turbulence_has_fuel(MANTA *smoke)
{
  return (smoke->getFuelHigh()) ? 1 : 0;
}
extern "C" int manta_smoke_turbulence_has_colors(MANTA *smoke)
{
  return (smoke->getColorRHigh() && smoke->getColorGHigh() && smoke->getColorBHigh()) ? 1 : 0;
}

extern "C" void manta_smoke_turbulence_get_res(MANTA *smoke, int *res)
{
  if (smoke && smoke->usingNoise()) {
    res[0] = smoke->getResXHigh();
    res[1] = smoke->getResYHigh();
    res[2] = smoke->getResZHigh();
  }
}
extern "C" int manta_smoke_turbulence_get_cells(MANTA *smoke)
{
  int total_cells_high = smoke->getResXHigh() * smoke->getResYHigh() * smoke->getResZHigh();
  return (smoke && smoke->usingNoise()) ? total_cells_high : 0;
}

/* Liquid functions */
extern "C" void manta_liquid_export_script(MANTA *liquid, MantaModifierData *mmd)
{
  if (!liquid || !mmd)
    return;
  liquid->exportLiquidScript(mmd);
}

extern "C" void manta_liquid_ensure_sndparts(MANTA *liquid, struct MantaModifierData *mmd)
{
  if (liquid) {
    liquid->initLiquidSndParts(mmd);
    liquid->updatePointers();
  }
}

/* Liquid accessors */
extern "C" int manta_liquid_get_particle_res_x(MANTA *liquid)
{
  return liquid->getParticleResX();
}
extern "C" int manta_liquid_get_particle_res_y(MANTA *liquid)
{
  return liquid->getParticleResY();
}
extern "C" int manta_liquid_get_particle_res_z(MANTA *liquid)
{
  return liquid->getParticleResZ();
}

extern "C" int manta_liquid_get_mesh_res_x(MANTA *liquid)
{
  return liquid->getMeshResX();
}
extern "C" int manta_liquid_get_mesh_res_y(MANTA *liquid)
{
  return liquid->getMeshResY();
}
extern "C" int manta_liquid_get_mesh_res_z(MANTA *liquid)
{
  return liquid->getMeshResZ();
}

extern "C" int manta_liquid_get_particle_upres(MANTA *liquid)
{
  return liquid->getParticleUpres();
}
extern "C" int manta_liquid_get_mesh_upres(MANTA *liquid)
{
  return liquid->getMeshUpres();
}

extern "C" int manta_liquid_get_num_verts(MANTA *liquid)
{
  return liquid->getNumVertices();
}
extern "C" int manta_liquid_get_num_normals(MANTA *liquid)
{
  return liquid->getNumNormals();
}
extern "C" int manta_liquid_get_num_triangles(MANTA *liquid)
{
  return liquid->getNumTriangles();
}

extern "C" float manta_liquid_get_vertex_x_at(MANTA *liquid, int i)
{
  return liquid->getVertexXAt(i);
}
extern "C" float manta_liquid_get_vertex_y_at(MANTA *liquid, int i)
{
  return liquid->getVertexYAt(i);
}
extern "C" float manta_liquid_get_vertex_z_at(MANTA *liquid, int i)
{
  return liquid->getVertexZAt(i);
}

extern "C" float manta_liquid_get_normal_x_at(MANTA *liquid, int i)
{
  return liquid->getNormalXAt(i);
}
extern "C" float manta_liquid_get_normal_y_at(MANTA *liquid, int i)
{
  return liquid->getNormalYAt(i);
}
extern "C" float manta_liquid_get_normal_z_at(MANTA *liquid, int i)
{
  return liquid->getNormalZAt(i);
}

extern "C" int manta_liquid_get_triangle_x_at(MANTA *liquid, int i)
{
  return liquid->getTriangleXAt(i);
}
extern "C" int manta_liquid_get_triangle_y_at(MANTA *liquid, int i)
{
  return liquid->getTriangleYAt(i);
}
extern "C" int manta_liquid_get_triangle_z_at(MANTA *liquid, int i)
{
  return liquid->getTriangleZAt(i);
}

extern "C" float manta_liquid_get_vertvel_x_at(MANTA *liquid, int i)
{
  return liquid->getVertVelXAt(i);
}
extern "C" float manta_liquid_get_vertvel_y_at(MANTA *liquid, int i)
{
  return liquid->getVertVelYAt(i);
}
extern "C" float manta_liquid_get_vertvel_z_at(MANTA *liquid, int i)
{
  return liquid->getVertVelZAt(i);
}

extern "C" int manta_liquid_get_num_flip_particles(MANTA *liquid)
{
  return liquid->getNumFlipParticles();
}
extern "C" int manta_liquid_get_num_snd_particles(MANTA *liquid)
{
  return liquid->getNumSndParticles();
}

extern "C" int manta_liquid_get_flip_particle_flag_at(MANTA *liquid, int i)
{
  return liquid->getFlipParticleFlagAt(i);
}
extern "C" int manta_liquid_get_snd_particle_flag_at(MANTA *liquid, int i)
{
  return liquid->getSndParticleFlagAt(i);
}

extern "C" float manta_liquid_get_flip_particle_position_x_at(MANTA *liquid, int i)
{
  return liquid->getFlipParticlePositionXAt(i);
}
extern "C" float manta_liquid_get_flip_particle_position_y_at(MANTA *liquid, int i)
{
  return liquid->getFlipParticlePositionYAt(i);
}
extern "C" float manta_liquid_get_flip_particle_position_z_at(MANTA *liquid, int i)
{
  return liquid->getFlipParticlePositionZAt(i);
}

extern "C" float manta_liquid_get_flip_particle_velocity_x_at(MANTA *liquid, int i)
{
  return liquid->getFlipParticleVelocityXAt(i);
}
extern "C" float manta_liquid_get_flip_particle_velocity_y_at(MANTA *liquid, int i)
{
  return liquid->getFlipParticleVelocityYAt(i);
}
extern "C" float manta_liquid_get_flip_particle_velocity_z_at(MANTA *liquid, int i)
{
  return liquid->getFlipParticleVelocityZAt(i);
}

extern "C" float manta_liquid_get_snd_particle_position_x_at(MANTA *liquid, int i)
{
  return liquid->getSndParticlePositionXAt(i);
}
extern "C" float manta_liquid_get_snd_particle_position_y_at(MANTA *liquid, int i)
{
  return liquid->getSndParticlePositionYAt(i);
}
extern "C" float manta_liquid_get_snd_particle_position_z_at(MANTA *liquid, int i)
{
  return liquid->getSndParticlePositionZAt(i);
}

extern "C" float manta_liquid_get_snd_particle_velocity_x_at(MANTA *liquid, int i)
{
  return liquid->getSndParticleVelocityXAt(i);
}
extern "C" float manta_liquid_get_snd_particle_velocity_y_at(MANTA *liquid, int i)
{
  return liquid->getSndParticleVelocityYAt(i);
}
extern "C" float manta_liquid_get_snd_particle_velocity_z_at(MANTA *liquid, int i)
{
  return liquid->getSndParticleVelocityZAt(i);
}
