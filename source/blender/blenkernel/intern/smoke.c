/*
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
 * The Original Code is Copyright (C) Blender Foundation.
 * All rights reserved.
 */

/** \file
 * \ingroup bke
 */

/* Part of the code copied from elbeem fluid library, copyright by Nils Thuerey */

#include "MEM_guardedalloc.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h> /* memset */

#include "BLI_blenlib.h"
#include "BLI_math.h"
#include "BLI_kdopbvh.h"
#include "BLI_threads.h"
#include "BLI_utildefines.h"

#include "DNA_anim_types.h"
#include "DNA_armature_types.h"
#include "DNA_constraint_types.h"
#include "DNA_customdata_types.h"
#include "DNA_light_types.h"
#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_modifier_types.h"
#include "DNA_object_types.h"
#include "DNA_particle_types.h"
#include "DNA_scene_types.h"
#include "DNA_smoke_types.h"

#include "BKE_appdir.h"
#include "BKE_animsys.h"
#include "BKE_armature.h"
#include "BKE_bvhutils.h"
#include "BKE_collision.h"
#include "BKE_colortools.h"
#include "BKE_constraint.h"
#include "BKE_customdata.h"
#include "BKE_deform.h"
#include "BKE_effect.h"
#include "BKE_library.h"
#include "BKE_mesh.h"
#include "BKE_mesh_runtime.h"
#include "BKE_modifier.h"
#include "BKE_object.h"
#include "BKE_particle.h"
#include "BKE_pointcache.h"
#include "BKE_scene.h"
#include "BKE_smoke.h"
#include "BKE_texture.h"

#include "DEG_depsgraph.h"
#include "DEG_depsgraph_query.h"

#include "RE_shader_ext.h"

#include "GPU_glew.h"

/* UNUSED so far, may be enabled later */
/* #define USE_SMOKE_COLLISION_DM */

//#define DEBUG_TIME

#include "manta_fluid_API.h"

#ifdef DEBUG_TIME
#  include "PIL_time.h"
#endif

#ifdef WITH_MANTA

#  include "BLI_task.h"
#  include "BLI_kdtree.h"
#  include "BLI_voxel.h"

static ThreadMutex object_update_lock = BLI_MUTEX_INITIALIZER;

struct Mesh;
struct Object;
struct Scene;
struct SmokeModifierData;

// timestep default value for nice appearance 0.1f
#  define DT_DEFAULT 0.1f

#  define ADD_IF_LOWER_POS(a, b) (min_ff((a) + (b), max_ff((a), (b))))
#  define ADD_IF_LOWER_NEG(a, b) (max_ff((a) + (b), min_ff((a), (b))))
#  define ADD_IF_LOWER(a, b) (((b) > 0) ? ADD_IF_LOWER_POS((a), (b)) : ADD_IF_LOWER_NEG((a), (b)))

#else /* WITH_MANTA */

/* Stubs to use when smoke is disabled */
void fluid_free(struct FLUID *UNUSED(fluid));
float *smoke_get_density(struct FLUID *UNUSED(fluid));

void fluid_free(struct FLUID *UNUSED(fluid))
{
}
float *smoke_get_density(struct FLUID *UNUSED(fluid))
{
  return NULL;
}
struct Mesh *smokeModifier_do(SmokeModifierData *UNUSED(smd),
                              Depsgraph *UNUSED(depsgraph),
                              Scene *UNUSED(scene),
                              Object *UNUSED(ob),
                              Mesh *UNUSED(me))
{
  return NULL;
}
float BKE_smoke_get_velocity_at(struct Object *UNUSED(ob),
                                float UNUSED(position[3]),
                                float UNUSED(velocity[3]))
{
  return 0.0f;
}

#endif /* WITH_MANTA */

#ifdef WITH_MANTA

void BKE_smoke_reallocate_fluid(SmokeDomainSettings *sds, int res[3], int free_old)
{
  if (free_old && sds->fluid) {
    fluid_free(sds->fluid);
  }
  if (!min_iii(res[0], res[1], res[2])) {
    sds->fluid = NULL;
    return;
  }

  sds->fluid = fluid_init(res, sds->smd);

  if (sds->flags & FLUID_DOMAIN_USE_NOISE) {
    sds->res_noise[0] = res[0] * sds->noise_scale;
    sds->res_noise[1] = res[1] * sds->noise_scale;
    sds->res_noise[2] = res[2] * sds->noise_scale;
  }
}

/* convert global position to domain cell space */
static void smoke_pos_to_cell(SmokeDomainSettings *sds, float pos[3])
{
  mul_m4_v3(sds->imat, pos);
  sub_v3_v3(pos, sds->p0);
  pos[0] *= 1.0f / sds->cell_size[0];
  pos[1] *= 1.0f / sds->cell_size[1];
  pos[2] *= 1.0f / sds->cell_size[2];
}

/* set domain transformations and base resolution from object mesh */
static void smoke_set_domain_from_mesh(SmokeDomainSettings *sds,
                                       Object *ob,
                                       Mesh *me,
                                       bool init_resolution)
{
  size_t i;
  float min[3] = {FLT_MAX, FLT_MAX, FLT_MAX}, max[3] = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
  float size[3];
  MVert *verts = me->mvert;
  float scale = 0.0;
  int res;

  res = sds->maxres;

  // get BB of domain
  for (i = 0; i < me->totvert; i++) {
    // min BB
    min[0] = MIN2(min[0], verts[i].co[0]);
    min[1] = MIN2(min[1], verts[i].co[1]);
    min[2] = MIN2(min[2], verts[i].co[2]);

    // max BB
    max[0] = MAX2(max[0], verts[i].co[0]);
    max[1] = MAX2(max[1], verts[i].co[1]);
    max[2] = MAX2(max[2], verts[i].co[2]);
  }

  /* set domain bounds */
  copy_v3_v3(sds->p0, min);
  copy_v3_v3(sds->p1, max);
  sds->dx = 1.0f / res;

  /* calculate domain dimensions */
  sub_v3_v3v3(size, max, min);
  if (init_resolution) {
    zero_v3_int(sds->base_res);
    copy_v3_v3(sds->cell_size, size);
  }
  /* apply object scale */
  for (i = 0; i < 3; i++) {
    size[i] = fabsf(size[i] * ob->scale[i]);
  }
  copy_v3_v3(sds->global_size, size);
  copy_v3_v3(sds->dp0, min);

  invert_m4_m4(sds->imat, ob->obmat);

  // prevent crash when initializing a plane as domain
  if (!init_resolution || (size[0] < FLT_EPSILON) || (size[1] < FLT_EPSILON) ||
      (size[2] < FLT_EPSILON)) {
    return;
  }

  /* define grid resolutions from longest domain side */
  if (size[0] >= MAX2(size[1], size[2])) {
    scale = res / size[0];
    sds->scale = size[0] / fabsf(ob->scale[0]);
    sds->base_res[0] = res;
    sds->base_res[1] = max_ii((int)(size[1] * scale + 0.5f), 4);
    sds->base_res[2] = max_ii((int)(size[2] * scale + 0.5f), 4);
  }
  else if (size[1] >= MAX2(size[0], size[2])) {
    scale = res / size[1];
    sds->scale = size[1] / fabsf(ob->scale[1]);
    sds->base_res[0] = max_ii((int)(size[0] * scale + 0.5f), 4);
    sds->base_res[1] = res;
    sds->base_res[2] = max_ii((int)(size[2] * scale + 0.5f), 4);
  }
  else {
    scale = res / size[2];
    sds->scale = size[2] / fabsf(ob->scale[2]);
    sds->base_res[0] = max_ii((int)(size[0] * scale + 0.5f), 4);
    sds->base_res[1] = max_ii((int)(size[1] * scale + 0.5f), 4);
    sds->base_res[2] = res;
  }

  /* set cell size */
  sds->cell_size[0] /= (float)sds->base_res[0];
  sds->cell_size[1] /= (float)sds->base_res[1];
  sds->cell_size[2] /= (float)sds->base_res[2];
}

static void smoke_set_domain_gravity(Scene *scene, SmokeDomainSettings *sds)
{
  float gravity[3] = {0.0f, 0.0f, -1.0f};
  float gravity_mag;

  /* use global gravity if enabled */
  if (scene->physics_settings.flag & PHYS_GLOBAL_GRAVITY) {
    copy_v3_v3(gravity, scene->physics_settings.gravity);
    /* map default value to 1.0 */
    mul_v3_fl(gravity, 1.0f / 9.810f);

    /* convert gravity to domain space */
    gravity_mag = len_v3(gravity);
    mul_mat3_m4_v3(sds->imat, gravity);
    normalize_v3(gravity);
    mul_v3_fl(gravity, gravity_mag);

    sds->gravity[0] = gravity[0];
    sds->gravity[1] = gravity[1];
    sds->gravity[2] = gravity[2];
  }
}

static int smokeModifier_init(
    SmokeModifierData *smd, Depsgraph *depsgraph, Object *ob, Scene *scene, Mesh *me)
{
  int scene_framenr = (int)DEG_get_ctime(depsgraph);

  if ((smd->type & MOD_SMOKE_TYPE_DOMAIN) && smd->domain && !smd->domain->fluid) {
    SmokeDomainSettings *sds = smd->domain;
    int res[3];
    /* set domain dimensions from mesh */
    smoke_set_domain_from_mesh(sds, ob, me, true);
    /* set domain gravity */
    smoke_set_domain_gravity(scene, sds);
    /* reset domain values */
    zero_v3_int(sds->shift);
    zero_v3(sds->shift_f);
    add_v3_fl(sds->shift_f, 0.5f);
    zero_v3(sds->prev_loc);
    mul_m4_v3(ob->obmat, sds->prev_loc);
    copy_m4_m4(sds->obmat, ob->obmat);

    /* set resolutions */
    if (smd->domain->flags & FLUID_DOMAIN_USE_ADAPTIVE_DOMAIN) {
      res[0] = res[1] = res[2] = 1; /* use minimum res for adaptive init */
    }
    else {
      copy_v3_v3_int(res, sds->base_res);
    }
    copy_v3_v3_int(sds->res, res);
    sds->total_cells = sds->res[0] * sds->res[1] * sds->res[2];
    sds->res_min[0] = sds->res_min[1] = sds->res_min[2] = 0;
    copy_v3_v3_int(sds->res_max, res);

    /* set time, dt = 0.1 is at 25fps */
    float fps = scene->r.frs_sec / scene->r.frs_sec_base;
    sds->dt = DT_DEFAULT * (25.0f / fps);

    /* allocate fluid */
    BKE_smoke_reallocate_fluid(sds, sds->res, 0);

    smd->time = scene_framenr;

    return 1;
  }
  else if (smd->type & MOD_SMOKE_TYPE_FLOW) {
    if (!smd->flow) {
      smokeModifier_createType(smd);
    }
    smd->time = scene_framenr;
    return 1;
  }
  else if (smd->type & MOD_SMOKE_TYPE_EFFEC) {
    if (!smd->effec) {
      smokeModifier_createType(smd);
    }
    smd->time = scene_framenr;
    return 1;
  }
  return 2;
}

#endif /* WITH_MANTA */

static void smokeModifier_freeDomain(SmokeModifierData *smd)
{
  if (smd->domain) {
    if (smd->domain->fluid) {
      fluid_free(smd->domain->fluid);
    }

    if (smd->domain->fluid_mutex) {
      BLI_rw_mutex_free(smd->domain->fluid_mutex);
    }

    if (smd->domain->effector_weights) {
      MEM_freeN(smd->domain->effector_weights);
    }
    smd->domain->effector_weights = NULL;

    if (!(smd->modifier.flag & eModifierFlag_SharedCaches)) {
      BKE_ptcache_free_list(&(smd->domain->ptcaches[0]));
      smd->domain->point_cache[0] = NULL;
    }

    if (smd->domain->mesh_velocities) {
      MEM_freeN(smd->domain->mesh_velocities);
    }
    smd->domain->mesh_velocities = NULL;

    if (smd->domain->coba) {
      MEM_freeN(smd->domain->coba);
    }

    MEM_freeN(smd->domain);
    smd->domain = NULL;
  }
}

static void smokeModifier_freeFlow(SmokeModifierData *smd)
{
  if (smd->flow) {
    if (smd->flow->mesh) {
      BKE_id_free(NULL, smd->flow->mesh);
    }
    smd->flow->mesh = NULL;

    if (smd->flow->verts_old) {
      MEM_freeN(smd->flow->verts_old);
    }
    smd->flow->verts_old = NULL;
    smd->flow->numverts = 0;

    MEM_freeN(smd->flow);
    smd->flow = NULL;
  }
}

static void smokeModifier_freeCollision(SmokeModifierData *smd)
{
  if (smd->effec) {
    if (smd->effec->mesh) {
      BKE_id_free(NULL, smd->effec->mesh);
    }
    smd->effec->mesh = NULL;

    if (smd->effec->verts_old) {
      MEM_freeN(smd->effec->verts_old);
    }
    smd->effec->verts_old = NULL;
    smd->effec->numverts = 0;

    MEM_freeN(smd->effec);
    smd->effec = NULL;
  }
}

static void smokeModifier_reset_ex(struct SmokeModifierData *smd, bool need_lock)
{
  if (smd) {
    if (smd->domain) {
      if (smd->domain->fluid) {
        if (need_lock) {
          BLI_rw_mutex_lock(smd->domain->fluid_mutex, THREAD_LOCK_WRITE);
        }

        fluid_free(smd->domain->fluid);
        smd->domain->fluid = NULL;

        if (need_lock) {
          BLI_rw_mutex_unlock(smd->domain->fluid_mutex);
        }
      }

      smd->time = -1;
      smd->domain->total_cells = 0;
      smd->domain->active_fields = 0;
    }
    else if (smd->flow) {
      if (smd->flow->verts_old) {
        MEM_freeN(smd->flow->verts_old);
      }
      smd->flow->verts_old = NULL;
      smd->flow->numverts = 0;
    }
    else if (smd->effec) {
      if (smd->effec->verts_old) {
        MEM_freeN(smd->effec->verts_old);
      }
      smd->effec->verts_old = NULL;
      smd->effec->numverts = 0;
    }
  }
}

void smokeModifier_reset(struct SmokeModifierData *smd)
{
  smokeModifier_reset_ex(smd, true);
}

void smokeModifier_free(SmokeModifierData *smd)
{
  if (smd) {
    smokeModifier_freeDomain(smd);
    smokeModifier_freeFlow(smd);
    smokeModifier_freeCollision(smd);
  }
}

void smokeModifier_createType(struct SmokeModifierData *smd)
{
  if (smd) {
    if (smd->type & MOD_SMOKE_TYPE_DOMAIN) {
      if (smd->domain) {
        smokeModifier_freeDomain(smd);
      }

      /* domain object data */
      smd->domain = MEM_callocN(sizeof(SmokeDomainSettings), "SmokeDomain");
      smd->domain->smd = smd;
      smd->domain->effector_weights = BKE_effector_add_weights(NULL);
      smd->domain->fluid = NULL;
      smd->domain->fluid_mutex = BLI_rw_mutex_alloc();
      smd->domain->eff_group = NULL;
      smd->domain->fluid_group = NULL;
      smd->domain->coll_group = NULL;

      /* adaptive domain options */
      smd->domain->adapt_margin = 4;
      smd->domain->adapt_res = 0;
      smd->domain->adapt_threshold = 0.02f;

      /* fluid domain options */
      smd->domain->maxres = 64;
      smd->domain->solver_res = 3;
      smd->domain->border_collisions = 0;  // open domain
      smd->domain->flags = FLUID_DOMAIN_USE_DISSOLVE_LOG | FLUID_DOMAIN_USE_ADAPTIVE_TIME;
      smd->domain->gravity[0] = 0.0f;
      smd->domain->gravity[1] = 0.0f;
      smd->domain->gravity[2] = -1.0f;
      smd->domain->active_fields = 0;
      smd->domain->type = FLUID_DOMAIN_TYPE_GAS;

      /* smoke domain options */
      smd->domain->alpha = -0.001;
      smd->domain->beta = 0.3;
      smd->domain->diss_speed = 5;
      smd->domain->vorticity = 0;
      smd->domain->active_color[0] = 0.0f;
      smd->domain->active_color[1] = 0.0f;
      smd->domain->active_color[2] = 0.0f;
      smd->domain->highres_sampling = SM_HRES_FULLSAMPLE;

      /* flame options */
      smd->domain->burning_rate = 0.75f;
      smd->domain->flame_smoke = 1.0f;
      smd->domain->flame_vorticity = 0.5f;
      smd->domain->flame_ignition = 1.5f;
      smd->domain->flame_max_temp = 3.0f;
      smd->domain->flame_smoke_color[0] = 0.7f;
      smd->domain->flame_smoke_color[1] = 0.7f;
      smd->domain->flame_smoke_color[2] = 0.7f;

      /* noise options */
      smd->domain->noise_strength = 1.0;
      smd->domain->noise_pos_scale = 2.0f;
      smd->domain->noise_time_anim = 0.1f;
      smd->domain->noise_scale = 2;
      smd->domain->noise_type = FLUID_NOISE_TYPE_WAVELET;

      /* liquid domain options */
      smd->domain->particle_randomness = 0.1f;
      smd->domain->particle_number = 2;
      smd->domain->particle_minimum = 8;
      smd->domain->particle_maximum = 16;
      smd->domain->particle_radius = 1.8f;
      smd->domain->particle_band_width = 3.0f;

      /* diffusion options*/
      smd->domain->surface_tension = 0.0f;
      smd->domain->viscosity_base = 1.0f;
      smd->domain->viscosity_exponent = 6.0f;
      smd->domain->domain_size = 0.5f;

      /* mesh options */
      smd->domain->mesh_velocities = NULL;
      smd->domain->mesh_concave_upper = 3.5f;
      smd->domain->mesh_concave_lower = 0.4f;
      smd->domain->mesh_smoothen_pos = 1;
      smd->domain->mesh_smoothen_neg = 1;
      smd->domain->mesh_scale = 2;
      smd->domain->totvert = 0;
      smd->domain->mesh_generator = FLUID_DOMAIN_MESH_IMPROVED;

      /* secondary particle options */
      smd->domain->sndparticle_tau_min_wc = 2.0;
      smd->domain->sndparticle_tau_max_wc = 8.0;
      smd->domain->sndparticle_tau_min_ta = 5.0;
      smd->domain->sndparticle_tau_max_ta = 20.0;
      smd->domain->sndparticle_tau_min_k = 1.0;
      smd->domain->sndparticle_tau_max_k = 5.0;
      smd->domain->sndparticle_k_wc = 200;
      smd->domain->sndparticle_k_ta = 40;
      smd->domain->sndparticle_k_b = 0.5;
      smd->domain->sndparticle_k_d = 0.6;
      smd->domain->sndparticle_l_min = 10.0;
      smd->domain->sndparticle_l_max = 25.0;
      smd->domain->sndparticle_boundary = SNDPARTICLE_BOUNDARY_DELETE;
      smd->domain->sndparticle_combined_export = SNDPARTICLE_COMBINED_EXPORT_OFF;
      smd->domain->sndparticle_potential_radius = 2;
      smd->domain->sndparticle_update_radius = 2;
      smd->domain->particle_type = 0;
      smd->domain->particle_scale = 1;

      /* fluid guiding options */
      smd->domain->guiding_parent = NULL;
      smd->domain->guiding_alpha = 2.0f;
      smd->domain->guiding_beta = 5;
      smd->domain->guiding_vel_factor = 2.0f;
      smd->domain->guide_res = NULL;
      smd->domain->guiding_source = FLUID_DOMAIN_GUIDING_SRC_DOMAIN;

      /* cache options */
      smd->domain->cache_frame_start = 1;
      smd->domain->cache_frame_end = 50;
      smd->domain->cache_frame_pause_data = 0;
      smd->domain->cache_frame_pause_noise = 0;
      smd->domain->cache_frame_pause_mesh = 0;
      smd->domain->cache_frame_pause_particles = 0;
      smd->domain->cache_frame_pause_guiding = 0;
      smd->domain->cache_flag = 0;
      smd->domain->cache_mesh_format = FLUID_DOMAIN_FILE_BIN_OBJECT;
      smd->domain->cache_data_format = FLUID_DOMAIN_FILE_UNI;
      smd->domain->cache_particle_format = FLUID_DOMAIN_FILE_UNI;
      smd->domain->cache_noise_format = FLUID_DOMAIN_FILE_UNI;
      modifier_path_init(smd->domain->cache_directory,
                         sizeof(smd->domain->cache_directory),
                         FLUID_DOMAIN_DIR_DEFAULT);

      /* time options */
      smd->domain->time_scale = 1.0;
      smd->domain->cfl_condition = 4.0;

      /* display options */
      smd->domain->slice_method = FLUID_DOMAIN_SLICE_VIEW_ALIGNED;
      smd->domain->axis_slice_method = AXIS_SLICE_FULL;
      smd->domain->slice_axis = 0;
      smd->domain->interp_method = 0;
      smd->domain->draw_velocity = false;
      smd->domain->slice_per_voxel = 5.0f;
      smd->domain->slice_depth = 0.5f;
      smd->domain->display_thickness = 1.0f;
      smd->domain->coba = NULL;
      smd->domain->vector_scale = 1.0f;
      smd->domain->vector_draw_type = VECTOR_DRAW_NEEDLE;
      smd->domain->use_coba = false;
      smd->domain->coba_field = FLUID_DOMAIN_FIELD_DENSITY;

      /* -- Deprecated / unsed options (below)-- */

      /* pointcache options */
      BLI_listbase_clear(&smd->domain->ptcaches[1]);
      smd->domain->point_cache[0] = BKE_ptcache_add(&(smd->domain->ptcaches[0]));
      smd->domain->point_cache[0]->flag |= PTCACHE_DISK_CACHE;
      smd->domain->point_cache[0]->step = 1;
      smd->domain->point_cache[1] = NULL; /* Deprecated */
      smd->domain->cache_comp = SM_CACHE_LIGHT;
      smd->domain->cache_high_comp = SM_CACHE_LIGHT;

      /* OpenVDB cache options */
#ifdef WITH_OPENVDB_BLOSC
      smd->domain->openvdb_comp = VDB_COMPRESSION_BLOSC;
#else
      smd->domain->openvdb_comp = VDB_COMPRESSION_ZIP;
#endif
      smd->domain->clipping = 1e-3f;
      smd->domain->data_depth = 0;
    }
    else if (smd->type & MOD_SMOKE_TYPE_FLOW) {
      if (smd->flow) {
        smokeModifier_freeFlow(smd);
      }

      /* flow object data */
      smd->flow = MEM_callocN(sizeof(SmokeFlowSettings), "SmokeFlow");
      smd->flow->smd = smd;
      smd->flow->mesh = NULL;
      smd->flow->psys = NULL;
      smd->flow->noise_texture = NULL;

      /* initial velocity */
      smd->flow->verts_old = NULL;
      smd->flow->numverts = 0;
      smd->flow->vel_multi = 1.0f;
      smd->flow->vel_normal = 0.0f;
      smd->flow->vel_random = 0.0f;
      smd->flow->vel_coord[0] = 0.0f;
      smd->flow->vel_coord[1] = 0.0f;
      smd->flow->vel_coord[2] = 0.0f;

      /* emission */
      smd->flow->density = 1.0f;
      smd->flow->color[0] = 0.7f;
      smd->flow->color[1] = 0.7f;
      smd->flow->color[2] = 0.7f;
      smd->flow->fuel_amount = 1.0f;
      smd->flow->temp = 1.0f;
      smd->flow->volume_density = 0.0f;
      smd->flow->surface_distance = 1.5f;
      smd->flow->particle_size = 1.0f;
      smd->flow->subframes = 0;

      /* texture control */
      smd->flow->source = FLUID_FLOW_SOURCE_MESH;
      smd->flow->texture_size = 1.0f;

      smd->flow->type = FLUID_FLOW_TYPE_SMOKE;
      smd->flow->behavior = FLUID_FLOW_BEHAVIOR_GEOMETRY;
      smd->flow->type = FLUID_FLOW_TYPE_SMOKE;
      smd->flow->flags = FLUID_FLOW_ABSOLUTE | FLUID_FLOW_USE_PART_SIZE | FLUID_FLOW_USE_INFLOW;
    }
    else if (smd->type & MOD_SMOKE_TYPE_EFFEC) {
      if (smd->effec) {
        smokeModifier_freeCollision(smd);
      }

      /* effector object data */
      smd->effec = MEM_callocN(sizeof(SmokeCollSettings), "SmokeColl");
      smd->effec->smd = smd;
      smd->effec->mesh = NULL;
      smd->effec->verts_old = NULL;
      smd->effec->numverts = 0;
      smd->effec->surface_distance = 0.5f;
      smd->effec->type = FLUID_EFFECTOR_TYPE_COLLISION;

      /* guiding options */
      smd->effec->guiding_mode = FLUID_EFFECTOR_GUIDING_MAXIMUM;
      smd->effec->vel_multi = 1.0f;
    }
  }
}

void smokeModifier_copy(const struct SmokeModifierData *smd,
                        struct SmokeModifierData *tsmd,
                        const int flag)
{
  tsmd->type = smd->type;
  tsmd->time = smd->time;

  smokeModifier_createType(tsmd);

  if (tsmd->domain) {
    SmokeDomainSettings *tsds = tsmd->domain;
    SmokeDomainSettings *sds = smd->domain;

    /* domain object data */
    tsds->fluid_group = sds->fluid_group;
    tsds->eff_group = sds->eff_group;
    tsds->coll_group = sds->coll_group;
    MEM_freeN(tsds->effector_weights);
    tsds->effector_weights = MEM_dupallocN(sds->effector_weights);

    /* adaptive domain options */
    tsds->adapt_margin = sds->adapt_margin;
    tsds->adapt_res = sds->adapt_res;
    tsds->adapt_threshold = sds->adapt_threshold;

    /* fluid domain options */
    tsds->maxres = sds->maxres;
    tsds->solver_res = sds->solver_res;
    tsds->border_collisions = sds->border_collisions;
    tsds->flags = sds->flags;
    tsds->gravity[0] = sds->gravity[0];
    tsds->gravity[1] = sds->gravity[1];
    tsds->gravity[2] = sds->gravity[2];
    tsds->active_fields = sds->active_fields;
    tsds->type = sds->type;

    /* smoke domain options */
    tsds->alpha = sds->alpha;
    tsds->beta = sds->beta;
    tsds->diss_speed = sds->diss_speed;
    tsds->vorticity = sds->vorticity;
    tsds->highres_sampling = sds->highres_sampling;

    /* flame options */
    tsds->burning_rate = sds->burning_rate;
    tsds->flame_smoke = sds->flame_smoke;
    tsds->flame_vorticity = sds->flame_vorticity;
    tsds->flame_ignition = sds->flame_ignition;
    tsds->flame_max_temp = sds->flame_max_temp;
    copy_v3_v3(tsds->flame_smoke_color, sds->flame_smoke_color);

    /* noise options */
    tsds->noise_strength = sds->noise_strength;
    tsds->noise_pos_scale = sds->noise_pos_scale;
    tsds->noise_time_anim = sds->noise_time_anim;
    tsds->noise_scale = sds->noise_scale;
    tsds->noise_type = sds->noise_type;

    /* liquid domain options */
    tsds->particle_randomness = sds->particle_randomness;
    tsds->particle_number = sds->particle_number;
    tsds->particle_minimum = sds->particle_minimum;
    tsds->particle_maximum = sds->particle_maximum;
    tsds->particle_radius = sds->particle_radius;
    tsds->particle_band_width = sds->particle_band_width;

    /* diffusion options*/
    tsds->surface_tension = sds->surface_tension;
    tsds->viscosity_base = sds->viscosity_base;
    tsds->viscosity_exponent = sds->viscosity_exponent;
    tsds->domain_size = sds->domain_size;

    /* mesh options */
    if (sds->mesh_velocities) {
      tsds->mesh_velocities = MEM_dupallocN(sds->mesh_velocities);
    }
    tsds->mesh_concave_upper = sds->mesh_concave_upper;
    tsds->mesh_concave_lower = sds->mesh_concave_lower;
    tsds->mesh_smoothen_pos = sds->mesh_smoothen_pos;
    tsds->mesh_smoothen_neg = sds->mesh_smoothen_neg;
    tsds->mesh_scale = sds->mesh_scale;
    tsds->totvert = sds->totvert;
    tsds->mesh_generator = sds->mesh_generator;

    /* secondary particle options */
    tsds->sndparticle_k_b = sds->sndparticle_k_b;
    tsds->sndparticle_k_d = sds->sndparticle_k_d;
    tsds->sndparticle_k_ta = sds->sndparticle_k_ta;
    tsds->sndparticle_k_wc = sds->sndparticle_k_wc;
    tsds->sndparticle_l_max = sds->sndparticle_l_max;
    tsds->sndparticle_l_min = sds->sndparticle_l_min;
    tsds->sndparticle_tau_max_k = sds->sndparticle_tau_max_k;
    tsds->sndparticle_tau_max_ta = sds->sndparticle_tau_max_ta;
    tsds->sndparticle_tau_max_wc = sds->sndparticle_tau_max_wc;
    tsds->sndparticle_tau_min_k = sds->sndparticle_tau_min_k;
    tsds->sndparticle_tau_min_ta = sds->sndparticle_tau_min_ta;
    tsds->sndparticle_tau_min_wc = sds->sndparticle_tau_min_wc;
    tsds->sndparticle_boundary = sds->sndparticle_boundary;
    tsds->sndparticle_combined_export = sds->sndparticle_combined_export;
    tsds->sndparticle_potential_radius = sds->sndparticle_potential_radius;
    tsds->sndparticle_update_radius = sds->sndparticle_update_radius;
    tsds->particle_type = sds->particle_type;
    tsds->particle_scale = sds->particle_scale;

    /* fluid guiding options */
    tsds->guiding_parent = sds->guiding_parent;
    tsds->guiding_alpha = sds->guiding_alpha;
    tsds->guiding_beta = sds->guiding_beta;
    tsds->guiding_vel_factor = sds->guiding_vel_factor;
    tsds->guide_res = sds->guide_res;
    tsds->guiding_source = sds->guiding_source;

    /* cache options */
    tsds->cache_frame_start = sds->cache_frame_start;
    tsds->cache_frame_end = sds->cache_frame_end;
    tsds->cache_frame_pause_data = sds->cache_frame_pause_data;
    tsds->cache_frame_pause_noise = sds->cache_frame_pause_noise;
    tsds->cache_frame_pause_mesh = sds->cache_frame_pause_mesh;
    tsds->cache_frame_pause_particles = sds->cache_frame_pause_particles;
    tsds->cache_frame_pause_guiding = sds->cache_frame_pause_guiding;
    tsds->cache_flag = sds->cache_flag;
    tsds->cache_mesh_format = sds->cache_mesh_format;
    tsds->cache_data_format = sds->cache_data_format;
    tsds->cache_particle_format = sds->cache_particle_format;
    tsds->cache_noise_format = sds->cache_noise_format;
    BLI_strncpy(tsds->cache_directory, sds->cache_directory, sizeof(tsds->cache_directory));

    /* time options */
    tsds->time_scale = sds->time_scale;
    tsds->cfl_condition = sds->cfl_condition;

    /* display options */
    tsds->slice_method = sds->slice_method;
    tsds->axis_slice_method = sds->axis_slice_method;
    tsds->slice_axis = sds->slice_axis;
    tsds->interp_method = sds->interp_method;
    tsds->draw_velocity = sds->draw_velocity;
    tsds->slice_per_voxel = sds->slice_per_voxel;
    tsds->slice_depth = sds->slice_depth;
    tsds->display_thickness = sds->display_thickness;
    if (sds->coba) {
      tsds->coba = MEM_dupallocN(sds->coba);
    }
    tsds->vector_scale = sds->vector_scale;
    tsds->vector_draw_type = sds->vector_draw_type;
    tsds->use_coba = sds->use_coba;
    tsds->coba_field = sds->coba_field;

    /* -- Deprecated / unsed options (below)-- */

    /* pointcache options */
    BKE_ptcache_free_list(&(tsds->ptcaches[0]));
    if (flag & LIB_ID_CREATE_NO_MAIN) {
      /* Share the cache with the original object's modifier. */
      tsmd->modifier.flag |= eModifierFlag_SharedCaches;
      tsds->point_cache[0] = sds->point_cache[0];
      tsds->ptcaches[0] = sds->ptcaches[0];
    }
    else {
      tsds->point_cache[0] = BKE_ptcache_copy_list(
          &(tsds->ptcaches[0]), &(sds->ptcaches[0]), flag);
    }

    /* OpenVDB cache options */
    tsds->openvdb_comp = sds->openvdb_comp;
    tsds->clipping = sds->clipping;
    tsds->data_depth = sds->data_depth;
  }
  else if (tsmd->flow) {
    SmokeFlowSettings *tsfs = tsmd->flow;
    SmokeFlowSettings *sfs = smd->flow;

    tsfs->psys = sfs->psys;
    tsfs->noise_texture = sfs->noise_texture;

    /* initial velocity */
    tsfs->vel_multi = sfs->vel_multi;
    tsfs->vel_normal = sfs->vel_normal;
    tsfs->vel_random = sfs->vel_random;
    tsfs->vel_coord[0] = sfs->vel_coord[0];
    tsfs->vel_coord[1] = sfs->vel_coord[1];
    tsfs->vel_coord[2] = sfs->vel_coord[2];

    /* emission */
    tsfs->density = sfs->density;
    copy_v3_v3(tsfs->color, sfs->color);
    tsfs->fuel_amount = sfs->fuel_amount;
    tsfs->temp = sfs->temp;
    tsfs->volume_density = sfs->volume_density;
    tsfs->surface_distance = sfs->surface_distance;
    tsfs->particle_size = sfs->particle_size;
    tsfs->subframes = sfs->subframes;

    /* texture control */
    tsfs->texture_size = sfs->texture_size;
    tsfs->texture_offset = sfs->texture_offset;
    BLI_strncpy(tsfs->uvlayer_name, sfs->uvlayer_name, sizeof(tsfs->uvlayer_name));
    tsfs->vgroup_density = sfs->vgroup_density;

    tsfs->type = sfs->type;
    tsfs->behavior = sfs->behavior;
    tsfs->source = sfs->source;
    tsfs->texture_type = sfs->texture_type;
    tsfs->flags = sfs->flags;
  }
  else if (tsmd->effec) {
    SmokeCollSettings *tscs = tsmd->effec;
    SmokeCollSettings *scs = smd->effec;

    tscs->surface_distance = scs->surface_distance;
    tscs->type = scs->type;

    /* guiding options */
    tscs->guiding_mode = scs->guiding_mode;
    tscs->vel_multi = scs->vel_multi;
  }
}

#ifdef WITH_MANTA

// forward declaration
static void smoke_calc_transparency(SmokeDomainSettings *sds, ViewLayer *view_layer);
static float calc_voxel_transp(
    float *result, float *input, int res[3], int *pixel, float *tRay, float correct);
static void update_mesh_distances(int index,
                                  float *mesh_distances,
                                  BVHTreeFromMesh *treeData,
                                  const float ray_start[3],
                                  float surface_thickness);

static int get_light(ViewLayer *view_layer, float *light)
{
  Base *base_tmp = NULL;
  int found_light = 0;

  // try to find a lamp, preferably local
  for (base_tmp = FIRSTBASE(view_layer); base_tmp; base_tmp = base_tmp->next) {
    if (base_tmp->object->type == OB_LAMP) {
      Light *la = base_tmp->object->data;

      if (la->type == LA_LOCAL) {
        copy_v3_v3(light, base_tmp->object->obmat[3]);
        return 1;
      }
      else if (!found_light) {
        copy_v3_v3(light, base_tmp->object->obmat[3]);
        found_light = 1;
      }
    }
  }

  return found_light;
}

/**********************************************************
 * Obstacles
 **********************************************************/

typedef struct ObstaclesFromDMData {
  SmokeDomainSettings *sds;
  SmokeCollSettings *scs;
  const MVert *mvert;
  const MLoop *mloop;
  const MLoopTri *looptri;
  BVHTreeFromMesh *tree;

  bool has_velocity;
  float *vert_vel;
  float *velocityX, *velocityY, *velocityZ;
  int *num_objects;
  float *distances_map;
} ObstaclesFromDMData;

static void obstacles_from_mesh_task_cb(void *__restrict userdata,
                                        const int z,
                                        const ParallelRangeTLS *__restrict UNUSED(tls))
{
  ObstaclesFromDMData *data = userdata;
  SmokeDomainSettings *sds = data->sds;

  /* slightly rounded-up sqrt(3 * (0.5)^2) == max. distance of cell boundary along the diagonal */
  const float surface_distance = 2.0f;  //0.867f;
  /* Note: Use larger surface distance to cover larger area with obvel. Manta will use these obvels and extrapolate them (inside and outside obstacle) */

  for (int x = sds->res_min[0]; x < sds->res_max[0]; x++) {
    for (int y = sds->res_min[1]; y < sds->res_max[1]; y++) {
      const int index = fluid_get_index(
          x - sds->res_min[0], sds->res[0], y - sds->res_min[1], sds->res[1], z - sds->res_min[2]);

      float ray_start[3] = {(float)x + 0.5f, (float)y + 0.5f, (float)z + 0.5f};
      BVHTreeNearest nearest = {0};
      nearest.index = -1;
      nearest.dist_sq = surface_distance *
                        surface_distance; /* find_nearest uses squared distance */
      bool hasIncObj = false;

      /* find the nearest point on the mesh */
      if (BLI_bvhtree_find_nearest(
              data->tree->tree, ray_start, &nearest, data->tree->nearest_callback, data->tree) !=
          -1) {
        const MLoopTri *lt = &data->looptri[nearest.index];
        float weights[3];
        int v1, v2, v3;

        /* calculate barycentric weights for nearest point */
        v1 = data->mloop[lt->tri[0]].v;
        v2 = data->mloop[lt->tri[1]].v;
        v3 = data->mloop[lt->tri[2]].v;
        interp_weights_tri_v3(
            weights, data->mvert[v1].co, data->mvert[v2].co, data->mvert[v3].co, nearest.co);

        if (data->has_velocity) {
          /* increase object count */
          data->num_objects[index]++;
          hasIncObj = true;

          /* apply object velocity */
          float hit_vel[3];
          interp_v3_v3v3v3(hit_vel,
                           &data->vert_vel[v1 * 3],
                           &data->vert_vel[v2 * 3],
                           &data->vert_vel[v3 * 3],
                           weights);

          /* Guiding has additional velocity multiplier */
          if (data->scs->type == FLUID_EFFECTOR_TYPE_GUIDE) {
            mul_v3_fl(hit_vel, data->scs->vel_multi);

            switch (data->scs->guiding_mode) {
              case FLUID_EFFECTOR_GUIDING_AVERAGED:
                data->velocityX[index] = (data->velocityX[index] + hit_vel[0]) * 0.5f;
                data->velocityY[index] = (data->velocityY[index] + hit_vel[1]) * 0.5f;
                data->velocityZ[index] = (data->velocityZ[index] + hit_vel[2]) * 0.5f;
                break;
              case FLUID_EFFECTOR_GUIDING_OVERRIDE:
                data->velocityX[index] = hit_vel[0];
                data->velocityY[index] = hit_vel[1];
                data->velocityZ[index] = hit_vel[2];
                break;
              case FLUID_EFFECTOR_GUIDING_MINIMUM:
                data->velocityX[index] = MIN2(fabsf(hit_vel[0]), fabsf(data->velocityX[index]));
                data->velocityY[index] = MIN2(fabsf(hit_vel[1]), fabsf(data->velocityY[index]));
                data->velocityZ[index] = MIN2(fabsf(hit_vel[2]), fabsf(data->velocityZ[index]));
                break;
              case FLUID_EFFECTOR_GUIDING_MAXIMUM:
              default:
                data->velocityX[index] = MAX2(fabsf(hit_vel[0]), fabsf(data->velocityX[index]));
                data->velocityY[index] = MAX2(fabsf(hit_vel[1]), fabsf(data->velocityY[index]));
                data->velocityZ[index] = MAX2(fabsf(hit_vel[2]), fabsf(data->velocityZ[index]));
                break;
            }
          }
          else {
            /* Apply (i.e. add) effector object velocity */
            data->velocityX[index] += (data->scs->type == FLUID_EFFECTOR_TYPE_GUIDE) ?
                                          hit_vel[0] * data->scs->vel_multi :
                                          hit_vel[0];
            data->velocityY[index] += (data->scs->type == FLUID_EFFECTOR_TYPE_GUIDE) ?
                                          hit_vel[1] * data->scs->vel_multi :
                                          hit_vel[1];
            data->velocityZ[index] += (data->scs->type == FLUID_EFFECTOR_TYPE_GUIDE) ?
                                          hit_vel[2] * data->scs->vel_multi :
                                          hit_vel[2];
            //printf("adding effector object vel: [%f, %f, %f], dx is: %f\n", hit_vel[0], hit_vel[1], hit_vel[2], sds->dx);
          }
        }
      }

      /* Get distance to mesh surface from both within and outside grid (mantaflow phi grid) */
      if (data->distances_map) {
        update_mesh_distances(
            index, data->distances_map, data->tree, ray_start, data->scs->surface_distance);

        /* Ensure that num objects are also counted inside object. But dont count twice (see object inc for nearest point) */
        if (data->distances_map[index] < 0 && !hasIncObj) {
          data->num_objects[index]++;
        }
      }
    }
  }
}

static void obstacles_from_mesh(Object *coll_ob,
                                SmokeDomainSettings *sds,
                                SmokeCollSettings *scs,
                                float *distances_map,
                                float *velocityX,
                                float *velocityY,
                                float *velocityZ,
                                int *num_objects,
                                float dt)
{
  if (!scs->mesh) {
    return;
  }
  {
    Mesh *me = NULL;
    MVert *mvert = NULL;
    const MLoopTri *looptri;
    const MLoop *mloop;
    BVHTreeFromMesh treeData = {NULL};
    int numverts, i;

    float *vert_vel = NULL;
    bool has_velocity = false;

    me = BKE_mesh_copy_for_eval(scs->mesh, true);

    /* Duplicate vertices to modify. */
    if (me->mvert) {
      me->mvert = MEM_dupallocN(me->mvert);
      CustomData_set_layer(&me->vdata, CD_MVERT, me->mvert);
    }

    BKE_mesh_ensure_normals(me);
    mvert = me->mvert;
    mloop = me->mloop;
    looptri = BKE_mesh_runtime_looptri_ensure(me);
    numverts = me->totvert;

    /* TODO (sebbas):
     * Make vert_vel init optional?
     * code is in trouble if the object moves but is declared as "does not move" */
    {
      vert_vel = MEM_callocN(sizeof(float) * numverts * 3, "smoke_obs_velocity");

      if (scs->numverts != numverts || !scs->verts_old) {
        if (scs->verts_old) {
          MEM_freeN(scs->verts_old);
        }

        scs->verts_old = MEM_callocN(sizeof(float) * numverts * 3, "smoke_obs_verts_old");
        scs->numverts = numverts;
      }
      else {
        has_velocity = true;
      }
    }

    /*  Transform collider vertices to
     *   domain grid space for fast lookups */
    for (i = 0; i < numverts; i++) {
      float n[3];
      float co[3];

      /* vert pos */
      mul_m4_v3(coll_ob->obmat, mvert[i].co);
      smoke_pos_to_cell(sds, mvert[i].co);

      /* vert normal */
      normal_short_to_float_v3(n, mvert[i].no);
      mul_mat3_m4_v3(coll_ob->obmat, n);
      mul_mat3_m4_v3(sds->imat, n);
      normalize_v3(n);
      normal_float_to_short_v3(mvert[i].no, n);

      /* vert velocity */
      add_v3fl_v3fl_v3i(co, mvert[i].co, sds->shift);
      if (has_velocity) {
        sub_v3_v3v3(&vert_vel[i * 3], co, &scs->verts_old[i * 3]);
        mul_v3_fl(&vert_vel[i * 3], sds->dx / dt);
      }
      copy_v3_v3(&scs->verts_old[i * 3], co);
    }

    if (BKE_bvhtree_from_mesh_get(&treeData, me, BVHTREE_FROM_LOOPTRI, 4)) {
      ObstaclesFromDMData data = {.sds = sds,
                                  .scs = scs,
                                  .mvert = mvert,
                                  .mloop = mloop,
                                  .looptri = looptri,
                                  .tree = &treeData,
                                  .has_velocity = has_velocity,
                                  .vert_vel = vert_vel,
                                  .velocityX = velocityX,
                                  .velocityY = velocityY,
                                  .velocityZ = velocityZ,
                                  .num_objects = num_objects,
                                  .distances_map = distances_map};
      ParallelRangeSettings settings;
      BLI_parallel_range_settings_defaults(&settings);
      settings.scheduling_mode = TASK_SCHEDULING_DYNAMIC;
      BLI_task_parallel_range(
          sds->res_min[2], sds->res_max[2], &data, obstacles_from_mesh_task_cb, &settings);
    }
    /* free bvh tree */
    free_bvhtree_from_mesh(&treeData);
    BKE_id_free(NULL, me);

    if (vert_vel) {
      MEM_freeN(vert_vel);
    }
    if (me->mvert) {
      MEM_freeN(me->mvert);
    }
  }
}

static void update_obstacleflags(SmokeDomainSettings *sds, Object **collobjs, int numcollobj)
{
  int active_fields = sds->active_fields;
  unsigned int collIndex;

  /* Monitor active fields based on flow settings */
  for (collIndex = 0; collIndex < numcollobj; collIndex++) {
    Object *collob = collobjs[collIndex];
    SmokeModifierData *smd2 = (SmokeModifierData *)modifiers_findByType(collob,
                                                                        eModifierType_Smoke);

    if ((smd2->type & MOD_SMOKE_TYPE_EFFEC) && smd2->effec) {
      SmokeCollSettings *scs = smd2->effec;
      if (!scs) {
        break;
      }
      if (scs->type == FLUID_EFFECTOR_TYPE_COLLISION) {
        active_fields |= FLUID_DOMAIN_ACTIVE_OBSTACLE;
      }
      if (scs->type == FLUID_EFFECTOR_TYPE_GUIDE) {
        active_fields |= FLUID_DOMAIN_ACTIVE_GUIDING;
      }
    }
  }
  /* Finally, initialize new data fields if any */
  if (active_fields & FLUID_DOMAIN_ACTIVE_OBSTACLE) {
    fluid_ensure_obstacle(sds->fluid, sds->smd);
  }
  if (active_fields & FLUID_DOMAIN_ACTIVE_GUIDING) {
    fluid_ensure_guiding(sds->fluid, sds->smd);
  }
  sds->active_fields = active_fields;
}

/* Animated obstacles: dx_step = ((x_new - x_old) / totalsteps) * substep */
static void update_obstacles(Depsgraph *depsgraph,
                             Scene *scene,
                             Object *ob,
                             SmokeDomainSettings *sds,
                             float time_per_frame,
                             float frame_length,
                             int frame,
                             float dt)
{
  Object **collobjs = NULL;
  unsigned int numcollobj = 0, collIndex = 0;

  collobjs = BKE_collision_objects_create(
      depsgraph, ob, sds->coll_group, &numcollobj, eModifierType_Smoke);

  /* Update all flow related flags and ensure that corresponding grids get initialized */
  update_obstacleflags(sds, collobjs, numcollobj);

  float *velx = fluid_get_ob_velocity_x(sds->fluid);
  float *vely = fluid_get_ob_velocity_y(sds->fluid);
  float *velz = fluid_get_ob_velocity_z(sds->fluid);
  float *velxGuide = fluid_get_guide_velocity_x(sds->fluid);
  float *velyGuide = fluid_get_guide_velocity_y(sds->fluid);
  float *velzGuide = fluid_get_guide_velocity_z(sds->fluid);
  float *velxOrig = fluid_get_velocity_x(sds->fluid);
  float *velyOrig = fluid_get_velocity_y(sds->fluid);
  float *velzOrig = fluid_get_velocity_z(sds->fluid);
  float *density = smoke_get_density(sds->fluid);
  float *fuel = smoke_get_fuel(sds->fluid);
  float *flame = smoke_get_flame(sds->fluid);
  float *r = smoke_get_color_r(sds->fluid);
  float *g = smoke_get_color_g(sds->fluid);
  float *b = smoke_get_color_b(sds->fluid);
  float *phiObsIn = fluid_get_phiobs_in(sds->fluid);
  float *phiGuideIn = fluid_get_phiguide_in(sds->fluid);
  int *obstacles = smoke_get_obstacle(sds->fluid);
  int *num_obstacles = fluid_get_num_obstacle(sds->fluid);
  int *num_guides = fluid_get_num_guide(sds->fluid);
  unsigned int z;

  /* Grid reset before writing again */
  for (z = 0; z < sds->res[0] * sds->res[1] * sds->res[2]; z++) {
    if (phiObsIn) {
      phiObsIn[z] = 9999;
    }
    if (phiGuideIn) {
      phiGuideIn[z] = 9999;
    }
    if (num_obstacles) {
      num_obstacles[z] = 0;
    }
    if (num_guides) {
      num_guides[z] = 0;
    }

    if (velx && vely && velz) {
      velx[z] = 0.0f;
      vely[z] = 0.0f;
      velz[z] = 0.0f;
    }
    if (velxGuide && velyGuide && velzGuide) {
      velxGuide[z] = 0.0f;
      velyGuide[z] = 0.0f;
      velzGuide[z] = 0.0f;
    }
  }

  /* Prepare grids from effector objects */
  for (collIndex = 0; collIndex < numcollobj; collIndex++) {
    Object *collob = collobjs[collIndex];
    SmokeModifierData *smd2 = (SmokeModifierData *)modifiers_findByType(collob,
                                                                        eModifierType_Smoke);

    // DG TODO: check if modifier is active?
    if ((smd2->type & MOD_SMOKE_TYPE_EFFEC) && smd2->effec) {
      SmokeCollSettings *scs = smd2->effec;

      /* Length of one frame. If using adaptive stepping, length is smaller than actual frame length */
      float adaptframe_length = time_per_frame / frame_length;

      /* Handle adaptive subframe (ie has subframe fraction). Need to set according scene subframe parameter */
      if (time_per_frame < frame_length) {
        scene->r.subframe = adaptframe_length;
        scene->r.cfra = frame - 1;
      }
      /* Handle absolute endframe (ie no subframe fraction). Need to set the scene subframe parameter to 0 and advance current scene frame */
      else {
        scene->r.subframe = 0.0f;
        scene->r.cfra = frame;
      }
      //printf("effector: frame: %d // scene current frame: %d // scene current subframe: %f\n", frame, scene->r.cfra, scene->r.subframe);

      /* TODO (sebbas): Using BKE_scene_frame_get(scene) instead of new DEG_get_ctime(depsgraph) as subframes dont work with the latter yet */
      BKE_object_modifier_update_subframe(
          depsgraph, scene, collob, true, 5, BKE_scene_frame_get(scene), eModifierType_Smoke);

      if (scs && (scs->type == FLUID_EFFECTOR_TYPE_COLLISION)) {
        obstacles_from_mesh(collob, sds, scs, phiObsIn, velx, vely, velz, num_obstacles, dt);
      }
      if (scs && (scs->type == FLUID_EFFECTOR_TYPE_GUIDE)) {
        obstacles_from_mesh(
            collob, sds, scs, phiGuideIn, velxGuide, velyGuide, velzGuide, num_guides, dt);
      }
    }
  }

  BKE_collision_objects_free(collobjs);

  /* obstacle cells should not contain any velocity from the smoke simulation */
  for (z = 0; z < sds->res[0] * sds->res[1] * sds->res[2]; z++) {
    if (obstacles[z] & 2)  // mantaflow convention: FlagObstacle
    {
      if (velxOrig && velyOrig && velzOrig) {
        velxOrig[z] = 0;
        velyOrig[z] = 0;
        velzOrig[z] = 0;
      }
      if (density) {
        density[z] = 0;
      }
      if (fuel) {
        fuel[z] = 0;
        flame[z] = 0;
      }
      if (r) {
        r[z] = 0;
        g[z] = 0;
        b[z] = 0;
      }
    }
    /* average velocities from multiple obstacles in one cell */
    if (num_obstacles && num_obstacles[z]) {
      velx[z] /= num_obstacles[z];
      vely[z] /= num_obstacles[z];
      velz[z] /= num_obstacles[z];
    }
    /* average velocities from multiple guides in one cell */
    if (num_guides && num_guides[z]) {
      velxGuide[z] /= num_guides[z];
      velyGuide[z] /= num_guides[z];
      velzGuide[z] /= num_guides[z];
    }
  }
}

/**********************************************************
 * Flow emission code
 **********************************************************/

typedef struct EmissionMap {
  float *influence;
  float *influence_high;
  float *velocity;
  float *distances;
  float *distances_high;
  int min[3], max[3], res[3];
  int hmin[3], hmax[3], hres[3];
  int total_cells, valid;
} EmissionMap;

static void em_boundInsert(EmissionMap *em, float point[3])
{
  int i = 0;
  if (!em->valid) {
    for (; i < 3; i++) {
      em->min[i] = (int)floor(point[i]);
      em->max[i] = (int)ceil(point[i]);
    }
    em->valid = 1;
  }
  else {
    for (; i < 3; i++) {
      if (point[i] < em->min[i]) {
        em->min[i] = (int)floor(point[i]);
      }
      if (point[i] > em->max[i]) {
        em->max[i] = (int)ceil(point[i]);
      }
    }
  }
}

static void clampBoundsInDomain(SmokeDomainSettings *sds,
                                int min[3],
                                int max[3],
                                float *min_vel,
                                float *max_vel,
                                int margin,
                                float dt)
{
  int i;
  for (i = 0; i < 3; i++) {
    int adapt = (sds->flags & FLUID_DOMAIN_USE_ADAPTIVE_DOMAIN) ? sds->adapt_res : 0;
    /* add margin */
    min[i] -= margin;
    max[i] += margin;

    /* adapt to velocity */
    if (min_vel && min_vel[i] < 0.0f) {
      min[i] += (int)floor(min_vel[i] * dt);
    }
    if (max_vel && max_vel[i] > 0.0f) {
      max[i] += (int)ceil(max_vel[i] * dt);
    }

    /* clamp within domain max size */
    CLAMP(min[i], -adapt, sds->base_res[i] + adapt);
    CLAMP(max[i], -adapt, sds->base_res[i] + adapt);
  }
}

static void em_allocateData(EmissionMap *em, bool use_velocity, int hires_mul)
{
  int i, res[3];

  for (i = 0; i < 3; i++) {
    res[i] = em->max[i] - em->min[i];
    if (res[i] <= 0) {
      return;
    }
  }
  em->total_cells = res[0] * res[1] * res[2];
  copy_v3_v3_int(em->res, res);

  em->influence = MEM_callocN(sizeof(float) * em->total_cells, "smoke_flow_influence");
  if (use_velocity) {
    em->velocity = MEM_callocN(sizeof(float) * em->total_cells * 3, "smoke_flow_velocity");
  }

  em->distances = MEM_callocN(sizeof(float) * em->total_cells, "fluid_flow_distances");
  memset(em->distances, 0x7f7f7f7f, sizeof(float) * em->total_cells);  // init to inf

  /* allocate high resolution map if required */
  if (hires_mul > 1) {
    int total_cells_high = em->total_cells * (hires_mul * hires_mul * hires_mul);

    for (i = 0; i < 3; i++) {
      em->hmin[i] = em->min[i] * hires_mul;
      em->hmax[i] = em->max[i] * hires_mul;
      em->hres[i] = em->res[i] * hires_mul;
    }

    em->influence_high = MEM_callocN(sizeof(float) * total_cells_high,
                                     "smoke_flow_influence_high");
    em->distances_high = MEM_callocN(sizeof(float) * total_cells_high,
                                     "fluid_flow_distances_high");
    memset(em->distances_high, 0x7f7f7f7f, sizeof(float) * total_cells_high);  // init to inf
  }
  em->valid = 1;
}

static void em_freeData(EmissionMap *em)
{
  if (em->influence) {
    MEM_freeN(em->influence);
  }
  if (em->influence_high) {
    MEM_freeN(em->influence_high);
  }
  if (em->velocity) {
    MEM_freeN(em->velocity);
  }
  if (em->distances) {
    MEM_freeN(em->distances);
  }
  if (em->distances_high) {
    MEM_freeN(em->distances_high);
  }
}

static void em_combineMaps(
    EmissionMap *output, EmissionMap *em2, int hires_multiplier, int additive, float sample_size)
{
  int i, x, y, z;

  /* copyfill input 1 struct and clear output for new allocation */
  EmissionMap em1;
  memcpy(&em1, output, sizeof(EmissionMap));
  memset(output, 0, sizeof(EmissionMap));

  for (i = 0; i < 3; i++) {
    if (em1.valid) {
      output->min[i] = MIN2(em1.min[i], em2->min[i]);
      output->max[i] = MAX2(em1.max[i], em2->max[i]);
    }
    else {
      output->min[i] = em2->min[i];
      output->max[i] = em2->max[i];
    }
  }
  /* allocate output map */
  em_allocateData(output, (em1.velocity || em2->velocity), hires_multiplier);

  /* base resolution inputs */
  for (x = output->min[0]; x < output->max[0]; x++) {
    for (y = output->min[1]; y < output->max[1]; y++) {
      for (z = output->min[2]; z < output->max[2]; z++) {
        int index_out = fluid_get_index(x - output->min[0],
                                        output->res[0],
                                        y - output->min[1],
                                        output->res[1],
                                        z - output->min[2]);

        /* initialize with first input if in range */
        if (x >= em1.min[0] && x < em1.max[0] && y >= em1.min[1] && y < em1.max[1] &&
            z >= em1.min[2] && z < em1.max[2]) {
          int index_in = fluid_get_index(
              x - em1.min[0], em1.res[0], y - em1.min[1], em1.res[1], z - em1.min[2]);

          /* values */
          output->influence[index_out] = em1.influence[index_in];
          output->distances[index_out] = em1.distances[index_in];
          if (output->velocity && em1.velocity) {
            copy_v3_v3(&output->velocity[index_out * 3], &em1.velocity[index_in * 3]);
          }
        }

        /* apply second input if in range */
        if (x >= em2->min[0] && x < em2->max[0] && y >= em2->min[1] && y < em2->max[1] &&
            z >= em2->min[2] && z < em2->max[2]) {
          int index_in = fluid_get_index(
              x - em2->min[0], em2->res[0], y - em2->min[1], em2->res[1], z - em2->min[2]);

          /* values */
          if (additive) {
            output->influence[index_out] += em2->influence[index_in] * sample_size;
          }
          else {
            output->influence[index_out] = MAX2(em2->influence[index_in],
                                                output->influence[index_out]);
          }
          output->distances[index_out] = MIN2(em2->distances[index_in],
                                              output->distances[index_out]);
          if (output->velocity && em2->velocity) {
            /* last sample replaces the velocity */
            output->velocity[index_out * 3] = ADD_IF_LOWER(output->velocity[index_out * 3],
                                                           em2->velocity[index_in * 3]);
            output->velocity[index_out * 3 + 1] = ADD_IF_LOWER(output->velocity[index_out * 3 + 1],
                                                               em2->velocity[index_in * 3 + 1]);
            output->velocity[index_out * 3 + 2] = ADD_IF_LOWER(output->velocity[index_out * 3 + 2],
                                                               em2->velocity[index_in * 3 + 2]);
          }
        }
      }  // low res loop
    }
  }

  /* initialize high resolution input if available */
  if (output->influence_high) {
    for (x = output->hmin[0]; x < output->hmax[0]; x++) {
      for (y = output->hmin[1]; y < output->hmax[1]; y++) {
        for (z = output->hmin[2]; z < output->hmax[2]; z++) {
          int index_out = fluid_get_index(x - output->hmin[0],
                                          output->hres[0],
                                          y - output->hmin[1],
                                          output->hres[1],
                                          z - output->hmin[2]);

          /* initialize with first input if in range */
          if (x >= em1.hmin[0] && x < em1.hmax[0] && y >= em1.hmin[1] && y < em1.hmax[1] &&
              z >= em1.hmin[2] && z < em1.hmax[2]) {
            int index_in = fluid_get_index(
                x - em1.hmin[0], em1.hres[0], y - em1.hmin[1], em1.hres[1], z - em1.hmin[2]);
            /* values */
            output->influence_high[index_out] = em1.influence_high[index_in];
          }

          /* apply second input if in range */
          if (x >= em2->hmin[0] && x < em2->hmax[0] && y >= em2->hmin[1] && y < em2->hmax[1] &&
              z >= em2->hmin[2] && z < em2->hmax[2]) {
            int index_in = fluid_get_index(
                x - em2->hmin[0], em2->hres[0], y - em2->hmin[1], em2->hres[1], z - em2->hmin[2]);

            /* values */
            if (additive) {
              output->influence_high[index_out] += em2->distances_high[index_in] * sample_size;
            }
            else {
              output->distances_high[index_out] = MAX2(em2->distances_high[index_in],
                                                       output->distances_high[index_out]);
            }
            output->distances_high[index_out] = MIN2(em2->distances_high[index_in],
                                                     output->distances_high[index_out]);
          }
        }  // high res loop
      }
    }
  }

  /* free original data */
  em_freeData(&em1);
}

typedef struct EmitFromParticlesData {
  SmokeFlowSettings *sfs;
  KDTree_3d *tree;
  int hires_multiplier;

  EmissionMap *em;
  float *particle_vel;
  float hr;

  int *min, *max, *res;

  float solid;
  float smooth;
  float hr_smooth;
} EmitFromParticlesData;

static void emit_from_particles_task_cb(void *__restrict userdata,
                                        const int z,
                                        const ParallelRangeTLS *__restrict UNUSED(tls))
{
  EmitFromParticlesData *data = userdata;
  SmokeFlowSettings *sfs = data->sfs;
  EmissionMap *em = data->em;
  const int hires_multiplier = data->hires_multiplier;

  for (int x = data->min[0]; x < data->max[0]; x++) {
    for (int y = data->min[1]; y < data->max[1]; y++) {
      /* take low res samples where possible */
      if (hires_multiplier <= 1 ||
          !(x % hires_multiplier || y % hires_multiplier || z % hires_multiplier)) {
        /* get low res space coordinates */
        const int lx = x / hires_multiplier;
        const int ly = y / hires_multiplier;
        const int lz = z / hires_multiplier;

        const int index = fluid_get_index(
            lx - em->min[0], em->res[0], ly - em->min[1], em->res[1], lz - em->min[2]);
        const float ray_start[3] = {((float)lx) + 0.5f, ((float)ly) + 0.5f, ((float)lz) + 0.5f};

        /* find particle distance from the kdtree */
        KDTreeNearest_3d nearest;
        const float range = data->solid + data->smooth;
        BLI_kdtree_3d_find_nearest(data->tree, ray_start, &nearest);

        if (nearest.dist < range) {
          em->influence[index] = (nearest.dist < data->solid) ?
                                     1.0f :
                                     (1.0f - (nearest.dist - data->solid) / data->smooth);
          /* Uses particle velocity as initial velocity for smoke */
          if (sfs->flags & FLUID_FLOW_INITVELOCITY &&
              (sfs->psys->part->phystype != PART_PHYS_NO)) {
            madd_v3_v3fl(
                &em->velocity[index * 3], &data->particle_vel[nearest.index * 3], sfs->vel_multi);
          }
        }
      }

      /* take high res samples if required */
      if (hires_multiplier > 1) {
        /* get low res space coordinates */
        const float lx = ((float)x) * data->hr;
        const float ly = ((float)y) * data->hr;
        const float lz = ((float)z) * data->hr;

        const int index = fluid_get_index(
            x - data->min[0], data->res[0], y - data->min[1], data->res[1], z - data->min[2]);
        const float ray_start[3] = {
            lx + 0.5f * data->hr, ly + 0.5f * data->hr, lz + 0.5f * data->hr};

        /* find particle distance from the kdtree */
        KDTreeNearest_3d nearest;
        const float range = data->solid + data->hr_smooth;
        BLI_kdtree_3d_find_nearest(data->tree, ray_start, &nearest);

        if (nearest.dist < range) {
          em->influence_high[index] = (nearest.dist < data->solid) ?
                                          1.0f :
                                          (1.0f - (nearest.dist - data->solid) / data->smooth);
        }
      }
    }
  }
}

static void emit_from_particles(Object *flow_ob,
                                SmokeDomainSettings *sds,
                                SmokeFlowSettings *sfs,
                                EmissionMap *em,
                                Depsgraph *depsgraph,
                                Scene *scene,
                                float dt)
{
  if (sfs && sfs->psys && sfs->psys->part &&
      ELEM(sfs->psys->part->type, PART_EMITTER, PART_FLUID))  // is particle system selected
  {
    ParticleSimulationData sim;
    ParticleSystem *psys = sfs->psys;
    float *particle_pos;
    float *particle_vel;
    int totpart = psys->totpart, totchild;
    int p = 0;
    int valid_particles = 0;
    int bounds_margin = 1;

    /* radius based flow */
    const float solid = sfs->particle_size * 0.5f;
    const float smooth = 0.5f; /* add 0.5 cells of linear falloff to reduce aliasing */
    int hires_multiplier = 1;
    KDTree_3d *tree = NULL;

    sim.depsgraph = depsgraph;
    sim.scene = scene;
    sim.ob = flow_ob;
    sim.psys = psys;
    sim.psys->lattice_deform_data = psys_create_lattice_deform_data(&sim);

    /* prepare curvemapping tables */
    if ((psys->part->child_flag & PART_CHILD_USE_CLUMP_CURVE) && psys->part->clumpcurve) {
      curvemapping_changed_all(psys->part->clumpcurve);
    }
    if ((psys->part->child_flag & PART_CHILD_USE_ROUGH_CURVE) && psys->part->roughcurve) {
      curvemapping_changed_all(psys->part->roughcurve);
    }
    if ((psys->part->child_flag & PART_CHILD_USE_TWIST_CURVE) && psys->part->twistcurve) {
      curvemapping_changed_all(psys->part->twistcurve);
    }

    /* initialize particle cache */
    if (psys->part->type == PART_HAIR) {
      // TODO: PART_HAIR not supported whatsoever
      totchild = 0;
    }
    else {
      totchild = psys->totchild * psys->part->disp / 100;
    }

    particle_pos = MEM_callocN(sizeof(float) * (totpart + totchild) * 3, "smoke_flow_particles");
    particle_vel = MEM_callocN(sizeof(float) * (totpart + totchild) * 3, "smoke_flow_particles");

    /* setup particle radius emission if enabled */
    if (sfs->flags & FLUID_FLOW_USE_PART_SIZE) {
      tree = BLI_kdtree_3d_new(psys->totpart + psys->totchild);

      /* check need for high resolution map */
      if ((sds->flags & FLUID_DOMAIN_USE_NOISE) && (sds->highres_sampling == SM_HRES_FULLSAMPLE)) {
        hires_multiplier = sds->noise_scale;
      }

      bounds_margin = (int)ceil(solid + smooth);
    }

    /* calculate local position for each particle */
    for (p = 0; p < totpart + totchild; p++) {
      ParticleKey state;
      float *pos;
      if (p < totpart) {
        if (psys->particles[p].flag & (PARS_NO_DISP | PARS_UNEXIST)) {
          continue;
        }
      }
      else {
        /* handle child particle */
        ChildParticle *cpa = &psys->child[p - totpart];
        if (psys->particles[cpa->parent].flag & (PARS_NO_DISP | PARS_UNEXIST)) {
          continue;
        }
      }

      state.time = DEG_get_ctime(depsgraph); /* use depsgraph time */
      if (psys_get_particle_state(&sim, p, &state, 0) == 0) {
        continue;
      }

      /* location */
      pos = &particle_pos[valid_particles * 3];
      copy_v3_v3(pos, state.co);
      smoke_pos_to_cell(sds, pos);

      /* velocity */
      copy_v3_v3(&particle_vel[valid_particles * 3], state.vel);
      mul_mat3_m4_v3(sds->imat, &particle_vel[valid_particles * 3]);

      if (sfs->flags & FLUID_FLOW_USE_PART_SIZE) {
        BLI_kdtree_3d_insert(tree, valid_particles, pos);
      }

      /* calculate emission map bounds */
      em_boundInsert(em, pos);
      valid_particles++;
    }

    /* set emission map */
    clampBoundsInDomain(sds, em->min, em->max, NULL, NULL, bounds_margin, dt);
    em_allocateData(em, sfs->flags & FLUID_FLOW_INITVELOCITY, hires_multiplier);

    if (!(sfs->flags & FLUID_FLOW_USE_PART_SIZE)) {
      for (p = 0; p < valid_particles; p++) {
        int cell[3];
        size_t i = 0;
        size_t index = 0;
        int badcell = 0;

        /* 1. get corresponding cell */
        cell[0] = floor(particle_pos[p * 3]) - em->min[0];
        cell[1] = floor(particle_pos[p * 3 + 1]) - em->min[1];
        cell[2] = floor(particle_pos[p * 3 + 2]) - em->min[2];
        /* check if cell is valid (in the domain boundary) */
        for (i = 0; i < 3; i++) {
          if ((cell[i] > em->res[i] - 1) || (cell[i] < 0)) {
            badcell = 1;
            break;
          }
        }
        if (badcell) {
          continue;
        }
        /* get cell index */
        index = fluid_get_index(cell[0], em->res[0], cell[1], em->res[1], cell[2]);
        /* Add influence to emission map */
        em->influence[index] = 1.0f;
        /* Uses particle velocity as initial velocity for smoke */
        if (sfs->flags & FLUID_FLOW_INITVELOCITY && (psys->part->phystype != PART_PHYS_NO)) {
          madd_v3_v3fl(&em->velocity[index * 3], &particle_vel[p * 3], sfs->vel_multi);
        }
      }  // particles loop
    }
    else if (valid_particles > 0) {  // FLUID_FLOW_USE_PART_SIZE
      int min[3], max[3], res[3];
      const float hr = 1.0f / ((float)hires_multiplier);
      /* slightly adjust high res antialias smoothness based on number of divisions
       * to allow smaller details but yet not differing too much from the low res size */
      const float hr_smooth = smooth * powf(hr, 1.0f / 3.0f);

      /* setup loop bounds */
      for (int i = 0; i < 3; i++) {
        min[i] = em->min[i] * hires_multiplier;
        max[i] = em->max[i] * hires_multiplier;
        res[i] = em->res[i] * hires_multiplier;
      }

      BLI_kdtree_3d_balance(tree);

      EmitFromParticlesData data = {
          .sfs = sfs,
          .tree = tree,
          .hires_multiplier = hires_multiplier,
          .hr = hr,
          .em = em,
          .particle_vel = particle_vel,
          .min = min,
          .max = max,
          .res = res,
          .solid = solid,
          .smooth = smooth,
          .hr_smooth = hr_smooth,
      };

      ParallelRangeSettings settings;
      BLI_parallel_range_settings_defaults(&settings);
      settings.scheduling_mode = TASK_SCHEDULING_DYNAMIC;
      BLI_task_parallel_range(min[2], max[2], &data, emit_from_particles_task_cb, &settings);
    }

    if (sfs->flags & FLUID_FLOW_USE_PART_SIZE) {
      BLI_kdtree_3d_free(tree);
    }

    /* free data */
    if (particle_pos) {
      MEM_freeN(particle_pos);
    }
    if (particle_vel) {
      MEM_freeN(particle_vel);
    }
  }
}

/* Calculate map of (minimum) distances to flow/obstacle surface. Distances outside mesh are positive, inside negative */
static void update_mesh_distances(int index,
                                  float *mesh_distances,
                                  BVHTreeFromMesh *treeData,
                                  const float ray_start[3],
                                  float surface_thickness)
{

  /* First pass: Raycasts in 26 directions (6 main axis + 12 quadrant diagonals (2D) + 8 octant diagonals (3D)) */
  float ray_dirs[26][3] = {
      {1.0f, 0.0f, 0.0f},   {0.0f, 1.0f, 0.0f},   {0.0f, 0.0f, 1.0f},  {-1.0f, 0.0f, 0.0f},
      {0.0f, -1.0f, 0.0f},  {0.0f, 0.0f, -1.0f},  {1.0f, 1.0f, 0.0f},  {1.0f, -1.0f, 0.0f},
      {-1.0f, 1.0f, 0.0f},  {-1.0f, -1.0f, 0.0f}, {1.0f, 0.0f, 1.0f},  {1.0f, 0.0f, -1.0f},
      {-1.0f, 0.0f, 1.0f},  {-1.0f, 0.0f, -1.0f}, {0.0f, 1.0f, 1.0f},  {0.0f, 1.0f, -1.0f},
      {0.0f, -1.0f, 1.0f},  {0.0f, -1.0f, -1.0f}, {1.0f, 1.0f, 1.0f},  {1.0f, -1.0f, 1.0f},
      {-1.0f, 1.0f, 1.0f},  {-1.0f, -1.0f, 1.0f}, {1.0f, 1.0f, -1.0f}, {1.0f, -1.0f, -1.0f},
      {-1.0f, 1.0f, -1.0f}, {-1.0f, -1.0f, -1.0f}};
  size_t ray_cnt = sizeof ray_dirs / sizeof ray_dirs[0];

  /* Count for ray misses (no face hit) and cases where ray direction matches face normal direction. */
  int miss_cnt = 0;
  int dir_cnt = 0;

  for (int i = 0; i < ray_cnt; i++) {
    BVHTreeRayHit hit_tree = {0};
    hit_tree.index = -1;
    hit_tree.dist = 9999;

    normalize_v3(ray_dirs[i]);
    BLI_bvhtree_ray_cast(treeData->tree,
                         ray_start,
                         ray_dirs[i],
                         0.0f,
                         &hit_tree,
                         treeData->raycast_callback,
                         treeData);

    /* Ray did not hit mesh. Current point definitely not inside mesh. Inside mesh all rays have to hit. */
    if (hit_tree.index == -1) {
      miss_cnt++;
      continue;
    }

    /* Ray and normal are in pointing opposite directions. */
    if (dot_v3v3(ray_dirs[i], hit_tree.no) <= 0) {
      dir_cnt++;
    }
  }

  /* Initialize grid points to -0.5 inside and 0.5 outside mesh.
   * Inside mesh: All rays have to hit (no misses) or all face normals have to match ray direction */
  if (mesh_distances[index] != -0.5f) {
    mesh_distances[index] = (miss_cnt > 0 || dir_cnt == ray_cnt) ? 0.5f : -0.5f;
  }

  /* Second pass: Ensure that single planes get initialized. */
  BVHTreeNearest nearest = {0};
  nearest.index = -1;
  nearest.dist_sq = surface_thickness * surface_thickness; /* find_nearest uses squared distance */

  if (BLI_bvhtree_find_nearest(
          treeData->tree, ray_start, &nearest, treeData->nearest_callback, treeData) != -1) {
    if (mesh_distances[index] != -0.5f) {
      mesh_distances[index] = -0.5f;
    }
  }
}

static void sample_mesh(SmokeFlowSettings *sfs,
                        const MVert *mvert,
                        const MLoop *mloop,
                        const MLoopTri *mlooptri,
                        const MLoopUV *mloopuv,
                        float *influence_map,
                        float *velocity_map,
                        int index,
                        const int base_res[3],
                        float flow_center[3],
                        BVHTreeFromMesh *treeData,
                        const float ray_start[3],
                        const float *vert_vel,
                        bool has_velocity,
                        int defgrp_index,
                        MDeformVert *dvert,
                        float x,
                        float y,
                        float z)
{
  float ray_dir[3] = {1.0f, 0.0f, 0.0f};
  BVHTreeRayHit hit = {0};
  BVHTreeNearest nearest = {0};

  float volume_factor = 0.0f;
  float sample_str = 0.0f;

  hit.index = -1;
  hit.dist = 9999;
  nearest.index = -1;
  nearest.dist_sq = sfs->surface_distance *
                    sfs->surface_distance; /* find_nearest uses squared distance */

  /* Check volume collision */
  if (sfs->volume_density) {
    if (BLI_bvhtree_ray_cast(treeData->tree,
                             ray_start,
                             ray_dir,
                             0.0f,
                             &hit,
                             treeData->raycast_callback,
                             treeData) != -1) {
      float dot = ray_dir[0] * hit.no[0] + ray_dir[1] * hit.no[1] + ray_dir[2] * hit.no[2];
      /* If ray and hit face normal are facing same direction
       * hit point is inside a closed mesh. */
      if (dot >= 0) {
        /* Also cast a ray in opposite direction to make sure
         * point is at least surrounded by two faces */
        negate_v3(ray_dir);
        hit.index = -1;
        hit.dist = 9999;

        BLI_bvhtree_ray_cast(
            treeData->tree, ray_start, ray_dir, 0.0f, &hit, treeData->raycast_callback, treeData);
        if (hit.index != -1) {
          volume_factor = sfs->volume_density;
        }
      }
    }
  }

  /* find the nearest point on the mesh */
  if (BLI_bvhtree_find_nearest(
          treeData->tree, ray_start, &nearest, treeData->nearest_callback, treeData) != -1) {
    float weights[3];
    int v1, v2, v3, f_index = nearest.index;
    float n1[3], n2[3], n3[3], hit_normal[3];

    /* emit from surface based on distance */
    if (sfs->surface_distance) {
      sample_str = sqrtf(nearest.dist_sq) / sfs->surface_distance;
      CLAMP(sample_str, 0.0f, 1.0f);
      sample_str = pow(1.0f - sample_str, 0.5f);
    }
    else {
      sample_str = 0.0f;
    }

    /* calculate barycentric weights for nearest point */
    v1 = mloop[mlooptri[f_index].tri[0]].v;
    v2 = mloop[mlooptri[f_index].tri[1]].v;
    v3 = mloop[mlooptri[f_index].tri[2]].v;
    interp_weights_tri_v3(weights, mvert[v1].co, mvert[v2].co, mvert[v3].co, nearest.co);

    if (sfs->flags & FLUID_FLOW_INITVELOCITY && velocity_map) {
      /* apply normal directional velocity */
      if (sfs->vel_normal) {
        /* interpolate vertex normal vectors to get nearest point normal */
        normal_short_to_float_v3(n1, mvert[v1].no);
        normal_short_to_float_v3(n2, mvert[v2].no);
        normal_short_to_float_v3(n3, mvert[v3].no);
        interp_v3_v3v3v3(hit_normal, n1, n2, n3, weights);
        normalize_v3(hit_normal);
        /* apply normal directional and random velocity
         * - TODO: random disabled for now since it doesn't really work well
         *   as pressure calc smoothens it out. */
        velocity_map[index * 3] += hit_normal[0] * sfs->vel_normal * 0.25f;
        velocity_map[index * 3 + 1] += hit_normal[1] * sfs->vel_normal * 0.25f;
        velocity_map[index * 3 + 2] += hit_normal[2] * sfs->vel_normal * 0.25f;
        /* TODO: for fire emitted from mesh surface we can use
         * Vf = Vs + (Ps/Pf - 1)*S to model gaseous expansion from solid to fuel */
      }
      /* apply object velocity */
      if (has_velocity && sfs->vel_multi) {
        float hit_vel[3];
        interp_v3_v3v3v3(
            hit_vel, &vert_vel[v1 * 3], &vert_vel[v2 * 3], &vert_vel[v3 * 3], weights);
        velocity_map[index * 3] += hit_vel[0] * sfs->vel_multi;
        velocity_map[index * 3 + 1] += hit_vel[1] * sfs->vel_multi;
        velocity_map[index * 3 + 2] += hit_vel[2] * sfs->vel_multi;
        //printf("adding flow object vel: [%f, %f, %f]\n", hit_vel[0], hit_vel[1], hit_vel[2]);
      }
      velocity_map[index * 3] += sfs->vel_coord[0];
      velocity_map[index * 3 + 1] += sfs->vel_coord[1];
      velocity_map[index * 3 + 2] += sfs->vel_coord[2];
    }

    /* apply vertex group influence if used */
    if (defgrp_index != -1 && dvert) {
      float weight_mask = defvert_find_weight(&dvert[v1], defgrp_index) * weights[0] +
                          defvert_find_weight(&dvert[v2], defgrp_index) * weights[1] +
                          defvert_find_weight(&dvert[v3], defgrp_index) * weights[2];
      sample_str *= weight_mask;
    }

    /* apply emission texture */
    if ((sfs->flags & FLUID_FLOW_TEXTUREEMIT) && sfs->noise_texture) {
      float tex_co[3] = {0};
      TexResult texres;

      if (sfs->texture_type == FLUID_FLOW_TEXTURE_MAP_AUTO) {
        tex_co[0] = ((x - flow_center[0]) / base_res[0]) / sfs->texture_size;
        tex_co[1] = ((y - flow_center[1]) / base_res[1]) / sfs->texture_size;
        tex_co[2] = ((z - flow_center[2]) / base_res[2] - sfs->texture_offset) / sfs->texture_size;
      }
      else if (mloopuv) {
        const float *uv[3];
        uv[0] = mloopuv[mlooptri[f_index].tri[0]].uv;
        uv[1] = mloopuv[mlooptri[f_index].tri[1]].uv;
        uv[2] = mloopuv[mlooptri[f_index].tri[2]].uv;

        interp_v2_v2v2v2(tex_co, UNPACK3(uv), weights);

        /* map between -1.0f and 1.0f */
        tex_co[0] = tex_co[0] * 2.0f - 1.0f;
        tex_co[1] = tex_co[1] * 2.0f - 1.0f;
        tex_co[2] = sfs->texture_offset;
      }
      texres.nor = NULL;
      BKE_texture_get_value(NULL, sfs->noise_texture, tex_co, &texres, false);
      sample_str *= texres.tin;
    }
  }

  /* multiply initial velocity by emitter influence */
  if (sfs->flags & FLUID_FLOW_INITVELOCITY && velocity_map) {
    mul_v3_fl(&velocity_map[index * 3], sample_str);
  }

  /* apply final influence based on volume factor */
  influence_map[index] = MAX2(volume_factor, sample_str);
}

typedef struct EmitFromDMData {
  SmokeDomainSettings *sds;
  SmokeFlowSettings *sfs;
  const MVert *mvert;
  const MLoop *mloop;
  const MLoopTri *mlooptri;
  const MLoopUV *mloopuv;
  MDeformVert *dvert;
  int defgrp_index;

  BVHTreeFromMesh *tree;
  int hires_multiplier;
  float hr;

  EmissionMap *em;
  bool has_velocity;
  float *vert_vel;

  float *flow_center;
  int *min, *max, *res;
} EmitFromDMData;

static void emit_from_mesh_task_cb(void *__restrict userdata,
                                   const int z,
                                   const ParallelRangeTLS *__restrict UNUSED(tls))
{
  EmitFromDMData *data = userdata;
  EmissionMap *em = data->em;
  const int hires_multiplier = data->hires_multiplier;

  for (int x = data->min[0]; x < data->max[0]; x++) {
    for (int y = data->min[1]; y < data->max[1]; y++) {
      /* take low res samples where possible */
      if (hires_multiplier <= 1 ||
          !(x % hires_multiplier || y % hires_multiplier || z % hires_multiplier)) {
        /* get low res space coordinates */
        const int lx = x / hires_multiplier;
        const int ly = y / hires_multiplier;
        const int lz = z / hires_multiplier;

        const int index = fluid_get_index(
            lx - em->min[0], em->res[0], ly - em->min[1], em->res[1], lz - em->min[2]);
        const float ray_start[3] = {((float)lx) + 0.5f, ((float)ly) + 0.5f, ((float)lz) + 0.5f};

        /* Emission for smoke and fire. Result in em->influence. Also, calculate invels */
        sample_mesh(data->sfs,
                    data->mvert,
                    data->mloop,
                    data->mlooptri,
                    data->mloopuv,
                    em->influence,
                    em->velocity,
                    index,
                    data->sds->base_res,
                    data->flow_center,
                    data->tree,
                    ray_start,
                    data->vert_vel,
                    data->has_velocity,
                    data->defgrp_index,
                    data->dvert,
                    (float)lx,
                    (float)ly,
                    (float)lz);

        /* Calculate levelset from meshes. Result in em->distances */
        update_mesh_distances(
            index, em->distances, data->tree, ray_start, data->sfs->surface_distance);
      }

      /* take high res samples if required */
      if (hires_multiplier > 1) {
        /* get low res space coordinates */
        const float lx = ((float)x) * data->hr;
        const float ly = ((float)y) * data->hr;
        const float lz = ((float)z) * data->hr;

        const int index = fluid_get_index(
            x - data->min[0], data->res[0], y - data->min[1], data->res[1], z - data->min[2]);
        const float ray_start[3] = {
            lx + 0.5f * data->hr,
            ly + 0.5f * data->hr,
            lz + 0.5f * data->hr,
        };

        /* Emission for smoke and fire high. Result in em->influence_high */
        if (data->sfs->type == FLUID_FLOW_TYPE_SMOKE || data->sfs->type == FLUID_FLOW_TYPE_FIRE ||
            data->sfs->type == FLUID_FLOW_TYPE_SMOKEFIRE) {
          sample_mesh(data->sfs,
                      data->mvert,
                      data->mloop,
                      data->mlooptri,
                      data->mloopuv,
                      em->influence_high,
                      NULL,
                      index,
                      data->sds->base_res,
                      data->flow_center,
                      data->tree,
                      ray_start,
                      data->vert_vel,
                      data->has_velocity,
                      data->defgrp_index,
                      data->dvert,
                      /* x,y,z needs to be always lowres */
                      lx,
                      ly,
                      lz);
        }
      }
    }
  }
}

static void emit_from_mesh(
    Object *flow_ob, SmokeDomainSettings *sds, SmokeFlowSettings *sfs, EmissionMap *em, float dt)
{
  if (sfs->mesh) {
    Mesh *me = NULL;
    MVert *mvert = NULL;
    const MLoopTri *mlooptri = NULL;
    const MLoop *mloop = NULL;
    const MLoopUV *mloopuv = NULL;
    MDeformVert *dvert = NULL;
    BVHTreeFromMesh treeData = {NULL};
    int numverts, i;

    float *vert_vel = NULL;
    bool has_velocity = false;

    int defgrp_index = sfs->vgroup_density - 1;
    float flow_center[3] = {0};
    int min[3], max[3], res[3];
    int hires_multiplier = 1;

    /* copy mesh for thread safety because we modify it,
     * main issue is its VertArray being modified, then replaced and freed
     */
    me = BKE_mesh_copy_for_eval(sfs->mesh, true);

    /* Duplicate vertices to modify. */
    if (me->mvert) {
      me->mvert = MEM_dupallocN(me->mvert);
      CustomData_set_layer(&me->vdata, CD_MVERT, me->mvert);
    }

    BKE_mesh_ensure_normals(me);
    mvert = me->mvert;
    mloop = me->mloop;
    mlooptri = BKE_mesh_runtime_looptri_ensure(me);
    numverts = me->totvert;
    dvert = CustomData_get_layer(&me->vdata, CD_MDEFORMVERT);
    mloopuv = CustomData_get_layer_named(&me->ldata, CD_MLOOPUV, sfs->uvlayer_name);

    if (sfs->flags & FLUID_FLOW_INITVELOCITY) {
      vert_vel = MEM_callocN(sizeof(float) * numverts * 3, "smoke_flow_velocity");

      if (sfs->numverts != numverts || !sfs->verts_old) {
        if (sfs->verts_old) {
          MEM_freeN(sfs->verts_old);
        }
        sfs->verts_old = MEM_callocN(sizeof(float) * numverts * 3, "smoke_flow_verts_old");
        sfs->numverts = numverts;
      }
      else {
        has_velocity = true;
      }
    }

    /*  Transform mesh vertices to
     *   domain grid space for fast lookups */
    for (i = 0; i < numverts; i++) {
      float n[3];

      /* vert pos */
      mul_m4_v3(flow_ob->obmat, mvert[i].co);
      smoke_pos_to_cell(sds, mvert[i].co);

      /* vert normal */
      normal_short_to_float_v3(n, mvert[i].no);
      mul_mat3_m4_v3(flow_ob->obmat, n);
      mul_mat3_m4_v3(sds->imat, n);
      normalize_v3(n);
      normal_float_to_short_v3(mvert[i].no, n);

      /* vert velocity */
      if (sfs->flags & FLUID_FLOW_INITVELOCITY) {
        float co[3];
        add_v3fl_v3fl_v3i(co, mvert[i].co, sds->shift);
        if (has_velocity) {
          sub_v3_v3v3(&vert_vel[i * 3], co, &sfs->verts_old[i * 3]);
          mul_v3_fl(&vert_vel[i * 3], sds->dx / dt);
        }
        copy_v3_v3(&sfs->verts_old[i * 3], co);
      }

      /* calculate emission map bounds */
      em_boundInsert(em, mvert[i].co);
    }
    mul_m4_v3(flow_ob->obmat, flow_center);
    smoke_pos_to_cell(sds, flow_center);

    /* check need for high resolution map */
    if ((sds->flags & FLUID_DOMAIN_USE_NOISE) && (sds->highres_sampling == SM_HRES_FULLSAMPLE)) {
      hires_multiplier = sds->noise_scale;
    }

    /* set emission map */
    clampBoundsInDomain(sds, em->min, em->max, NULL, NULL, (int)ceil(sfs->surface_distance), dt);
    em_allocateData(em, sfs->flags & FLUID_FLOW_INITVELOCITY, hires_multiplier);

    /* setup loop bounds */
    for (i = 0; i < 3; i++) {
      min[i] = em->min[i] * hires_multiplier;
      max[i] = em->max[i] * hires_multiplier;
      res[i] = em->res[i] * hires_multiplier;
    }

    if (BKE_bvhtree_from_mesh_get(&treeData, me, BVHTREE_FROM_LOOPTRI, 4)) {
      const float hr = 1.0f / ((float)hires_multiplier);

      EmitFromDMData data = {
          .sds = sds,
          .sfs = sfs,
          .mvert = mvert,
          .mloop = mloop,
          .mlooptri = mlooptri,
          .mloopuv = mloopuv,
          .dvert = dvert,
          .defgrp_index = defgrp_index,
          .tree = &treeData,
          .hires_multiplier = hires_multiplier,
          .hr = hr,
          .em = em,
          .has_velocity = has_velocity,
          .vert_vel = vert_vel,
          .flow_center = flow_center,
          .min = min,
          .max = max,
          .res = res,
      };

      ParallelRangeSettings settings;
      BLI_parallel_range_settings_defaults(&settings);
      settings.scheduling_mode = TASK_SCHEDULING_DYNAMIC;
      BLI_task_parallel_range(min[2], max[2], &data, emit_from_mesh_task_cb, &settings);
    }
    /* free bvh tree */
    free_bvhtree_from_mesh(&treeData);

    if (vert_vel) {
      MEM_freeN(vert_vel);
    }

    if (me->mvert) {
      MEM_freeN(me->mvert);
    }
    BKE_id_free(NULL, me);
  }
}

/**********************************************************
 *  Smoke step
 **********************************************************/

static void adaptiveDomainCopy(SmokeDomainSettings *sds,
                               int o_res[3],
                               int n_res[3],
                               int o_min[3],
                               int n_min[3],
                               int o_max[3],
                               int o_shift[3],
                               int n_shift[3],
                               int isNoise)
{
  int x, y, z;
  struct FLUID *fluid_old = sds->fluid;
  const int block_size = sds->noise_scale;
  int new_shift[3] = {0};
  sub_v3_v3v3_int(new_shift, n_shift, o_shift);

  /* allocate new fluid data */
  BKE_smoke_reallocate_fluid(sds, n_res, 0);

  int o_total_cells = o_res[0] * o_res[1] * o_res[2];
  int n_total_cells = n_res[0] * n_res[1] * n_res[2];

  /* copy values from old fluid to new */
  if (o_total_cells > 1 && n_total_cells > 1) {
    /* base smoke */
    float *o_dens, *o_react, *o_flame, *o_fuel, *o_heat, *o_vx, *o_vy, *o_vz,
      *o_r, *o_g, *o_b, *o_phiin;
    float *n_dens, *n_react, *n_flame, *n_fuel, *n_heat, *n_vx, *n_vy, *n_vz,
      *n_r, *n_g, *n_b, *n_phiin;
    float dummy, *dummy_s;
    int *dummy_p;
    /* noise smoke */
    int wt_res_old[3];
    float *o_wt_dens, *o_wt_react, *o_wt_flame, *o_wt_fuel, *o_wt_tcu, *o_wt_tcv, *o_wt_tcw,
      *o_wt_tcu2, *o_wt_tcv2, *o_wt_tcw2, *o_wt_r, *o_wt_g, *o_wt_b;
    float *n_wt_dens, *n_wt_react, *n_wt_flame, *n_wt_fuel, *n_wt_tcu, *n_wt_tcv, *n_wt_tcw,
      *n_wt_tcu2, *n_wt_tcv2, *n_wt_tcw2, *n_wt_r, *n_wt_g, *n_wt_b;

    if (isNoise && sds->flags & FLUID_DOMAIN_USE_NOISE) {
      smoke_turbulence_export(fluid_old,
                  &o_wt_dens,
                  &o_wt_react,
                  &o_wt_flame,
                  &o_wt_fuel,
                  &o_wt_r,
                  &o_wt_g,
                  &o_wt_b,
                  &o_wt_tcu,
                  &o_wt_tcv,
                  &o_wt_tcw,
                  &o_wt_tcu2,
                  &o_wt_tcv2,
                  &o_wt_tcw2);
      smoke_turbulence_get_res(fluid_old, wt_res_old);
      smoke_turbulence_export(sds->fluid,
                  &n_wt_dens,
                  &n_wt_react,
                  &n_wt_flame,
                  &n_wt_fuel,
                  &n_wt_r,
                  &n_wt_g,
                  &n_wt_b,
                  &n_wt_tcu,
                  &n_wt_tcv,
                  &n_wt_tcw,
                  &n_wt_tcu2,
                  &n_wt_tcv2,
                  &n_wt_tcw2);
    }
    else {
      smoke_export(fluid_old,
           &dummy,
           &dummy,
           &o_dens,
           &o_react,
           &o_flame,
           &o_fuel,
           &o_heat,
           &o_vx,
           &o_vy,
           &o_vz,
           &o_r,
           &o_g,
           &o_b,
           &dummy_p,
           &dummy_s,
           &o_phiin);
      smoke_export(sds->fluid,
           &dummy,
           &dummy,
           &n_dens,
           &n_react,
           &n_flame,
           &n_fuel,
           &n_heat,
           &n_vx,
           &n_vy,
           &n_vz,
           &n_r,
           &n_g,
           &n_b,
           &dummy_p,
           &dummy_s,
           &n_phiin);
    }

    for (x = o_min[0]; x < o_max[0]; x++) {
      for (y = o_min[1]; y < o_max[1]; y++) {
        for (z = o_min[2]; z < o_max[2]; z++) {
          /* old grid index */
          int xo = x - o_min[0];
          int yo = y - o_min[1];
          int zo = z - o_min[2];
          int index_old = fluid_get_index(xo, o_res[0], yo, o_res[1], zo);
          /* new grid index */
          int xn = x - n_min[0] - new_shift[0];
          int yn = y - n_min[1] - new_shift[1];
          int zn = z - n_min[2] - new_shift[2];
          int index_new = fluid_get_index(xn, n_res[0], yn, n_res[1], zn);

          /* skip if outside new domain */
          if (xn < 0 || xn >= n_res[0] || yn < 0 || yn >= n_res[1] || zn < 0 || zn >= n_res[2]) {
            continue;
          }

          /* copy data */
          if (isNoise && sds->flags & FLUID_DOMAIN_USE_NOISE) {
            int i, j, k;
            /* old grid index */
            int xx_o = xo * block_size;
            int yy_o = yo * block_size;
            int zz_o = zo * block_size;
            /* new grid index */
            int xx_n = xn * block_size;
            int yy_n = yn * block_size;
            int zz_n = zn * block_size;

            n_wt_tcu[index_new] = o_wt_tcu[index_old];
            n_wt_tcv[index_new] = o_wt_tcv[index_old];
            n_wt_tcw[index_new] = o_wt_tcw[index_old];

            n_wt_tcu2[index_new] = o_wt_tcu2[index_old];
            n_wt_tcv2[index_new] = o_wt_tcv2[index_old];
            n_wt_tcw2[index_new] = o_wt_tcw2[index_old];

            for (i = 0; i < block_size; i++) {
              for (j = 0; j < block_size; j++) {
                for (k = 0; k < block_size; k++) {
                  int big_index_old = fluid_get_index(
                    xx_o + i, wt_res_old[0], yy_o + j, wt_res_old[1], zz_o + k);
                  int big_index_new = fluid_get_index(
                    xx_n + i, sds->res_noise[0], yy_n + j, sds->res_noise[1], zz_n + k);
                  /* copy data */
                  n_wt_dens[big_index_new] = o_wt_dens[big_index_old];
                  if (n_wt_flame && o_wt_flame) {
                    n_wt_flame[big_index_new] = o_wt_flame[big_index_old];
                    n_wt_fuel[big_index_new] = o_wt_fuel[big_index_old];
                    n_wt_react[big_index_new] = o_wt_react[big_index_old];
                  }
                  if (n_wt_r && o_wt_r) {
                    n_wt_r[big_index_new] = o_wt_r[big_index_old];
                    n_wt_g[big_index_new] = o_wt_g[big_index_old];
                    n_wt_b[big_index_new] = o_wt_b[big_index_old];
                  }
                }
              }
            }
          }
          else {
            n_dens[index_new] = o_dens[index_old];
            /* heat */
            if (n_heat && o_heat) {
              n_heat[index_new] = o_heat[index_old];
            }
            /* fuel */
            if (n_fuel && o_fuel) {
              n_flame[index_new] = o_flame[index_old];
              n_fuel[index_new] = o_fuel[index_old];
              n_react[index_new] = o_react[index_old];
            }
            /* color */
            if (o_r && n_r) {
              n_r[index_new] = o_r[index_old];
              n_g[index_new] = o_g[index_old];
              n_b[index_new] = o_b[index_old];
            }
            n_vx[index_new] = o_vx[index_old];
            n_vy[index_new] = o_vy[index_old];
            n_vz[index_new] = o_vz[index_old];
            /* levelset */
            if (n_phiin && o_phiin) {
              n_phiin[index_new] = o_phiin[index_old];
            }
          }
        }
      }
    }
  }
  fluid_free(fluid_old);
}

static void adaptiveDomainAdjust(SmokeDomainSettings *sds,
                                 Object *ob,
                                 EmissionMap *emaps,
                                 unsigned int numflowobj,
                                 float dt)
{
  /* calculate domain shift for current frame */
  int new_shift[3] = {0};
  int total_shift[3];
  float frame_shift_f[3];
  float ob_loc[3] = {0};

  mul_m4_v3(ob->obmat, ob_loc);

  sub_v3_v3v3(frame_shift_f, ob_loc, sds->prev_loc);
  copy_v3_v3(sds->prev_loc, ob_loc);
  /* convert global space shift to local "cell" space */
  mul_mat3_m4_v3(sds->imat, frame_shift_f);
  frame_shift_f[0] = frame_shift_f[0] / sds->cell_size[0];
  frame_shift_f[1] = frame_shift_f[1] / sds->cell_size[1];
  frame_shift_f[2] = frame_shift_f[2] / sds->cell_size[2];
  /* add to total shift */
  add_v3_v3(sds->shift_f, frame_shift_f);
  /* convert to integer */
  total_shift[0] = (int)(floorf(sds->shift_f[0]));
  total_shift[1] = (int)(floorf(sds->shift_f[1]));
  total_shift[2] = (int)(floorf(sds->shift_f[2]));
  int tmpShift[3];
  copy_v3_v3_int(tmpShift, sds->shift);
  sub_v3_v3v3_int(new_shift, total_shift, sds->shift);
  copy_v3_v3_int(sds->shift, total_shift);

  /* calculate new domain boundary points so that smoke doesn't slide on sub-cell movement */
  sds->p0[0] = sds->dp0[0] - sds->cell_size[0] * (sds->shift_f[0] - total_shift[0] - 0.5f);
  sds->p0[1] = sds->dp0[1] - sds->cell_size[1] * (sds->shift_f[1] - total_shift[1] - 0.5f);
  sds->p0[2] = sds->dp0[2] - sds->cell_size[2] * (sds->shift_f[2] - total_shift[2] - 0.5f);
  sds->p1[0] = sds->p0[0] + sds->cell_size[0] * sds->base_res[0];
  sds->p1[1] = sds->p0[1] + sds->cell_size[1] * sds->base_res[1];
  sds->p1[2] = sds->p0[2] + sds->cell_size[2] * sds->base_res[2];

  /* adjust domain resolution */
  const int block_size = sds->noise_scale;
  int min[3] = {32767, 32767, 32767}, max[3] = {-32767, -32767, -32767}, res[3];
  int total_cells = 1, res_changed = 0, shift_changed = 0;
  float min_vel[3], max_vel[3];
  int x, y, z;
  float *density = smoke_get_density(sds->fluid);
  float *fuel = smoke_get_fuel(sds->fluid);
  float *bigdensity = smoke_turbulence_get_density(sds->fluid);
  float *bigfuel = smoke_turbulence_get_fuel(sds->fluid);
  float *vx = fluid_get_velocity_x(sds->fluid);
  float *vy = fluid_get_velocity_y(sds->fluid);
  float *vz = fluid_get_velocity_z(sds->fluid);
  int wt_res[3];

  if (sds->flags & FLUID_DOMAIN_USE_NOISE && sds->fluid) {
    smoke_turbulence_get_res(sds->fluid, wt_res);
  }

  INIT_MINMAX(min_vel, max_vel);

  /* Calculate bounds for current domain content */
  for (x = sds->res_min[0]; x < sds->res_max[0]; x++) {
    for (y = sds->res_min[1]; y < sds->res_max[1]; y++) {
      for (z = sds->res_min[2]; z < sds->res_max[2]; z++) {
        int xn = x - new_shift[0];
        int yn = y - new_shift[1];
        int zn = z - new_shift[2];
        int index;
        float max_den;

        /* skip if cell already belongs to new area */
        if (xn >= min[0] && xn <= max[0] && yn >= min[1] && yn <= max[1] && zn >= min[2] &&
            zn <= max[2]) {
          continue;
        }

        index = fluid_get_index(x - sds->res_min[0],
                                sds->res[0],
                                y - sds->res_min[1],
                                sds->res[1],
                                z - sds->res_min[2]);
        max_den = (fuel) ? MAX2(density[index], fuel[index]) : density[index];

        /* check high resolution bounds if max density isnt already high enough */
        if (max_den < sds->adapt_threshold && sds->flags & FLUID_DOMAIN_USE_NOISE && sds->fluid) {
          int i, j, k;
          /* high res grid index */
          int xx = (x - sds->res_min[0]) * block_size;
          int yy = (y - sds->res_min[1]) * block_size;
          int zz = (z - sds->res_min[2]) * block_size;

          for (i = 0; i < block_size; i++) {
            for (j = 0; j < block_size; j++) {
              for (k = 0; k < block_size; k++) {
                int big_index = fluid_get_index(xx + i, wt_res[0], yy + j, wt_res[1], zz + k);
                float den = (bigfuel) ? MAX2(bigdensity[big_index], bigfuel[big_index]) :
                                        bigdensity[big_index];
                if (den > max_den) {
                  max_den = den;
                }
              }
            }
          }
        }

        /* content bounds (use shifted coordinates) */
        if (max_den >= sds->adapt_threshold) {
          if (min[0] > xn) {
            min[0] = xn;
          }
          if (min[1] > yn) {
            min[1] = yn;
          }
          if (min[2] > zn) {
            min[2] = zn;
          }
          if (max[0] < xn) {
            max[0] = xn;
          }
          if (max[1] < yn) {
            max[1] = yn;
          }
          if (max[2] < zn) {
            max[2] = zn;
          }
        }

        /* velocity bounds */
        if (min_vel[0] > vx[index]) {
          min_vel[0] = vx[index];
        }
        if (min_vel[1] > vy[index]) {
          min_vel[1] = vy[index];
        }
        if (min_vel[2] > vz[index]) {
          min_vel[2] = vz[index];
        }
        if (max_vel[0] < vx[index]) {
          max_vel[0] = vx[index];
        }
        if (max_vel[1] < vy[index]) {
          max_vel[1] = vy[index];
        }
        if (max_vel[2] < vz[index]) {
          max_vel[2] = vz[index];
        }
      }
    }
  }

  /* also apply emission maps */
  for (int i = 0; i < numflowobj; i++) {
    EmissionMap *em = &emaps[i];

    for (x = em->min[0]; x < em->max[0]; x++) {
      for (y = em->min[1]; y < em->max[1]; y++) {
        for (z = em->min[2]; z < em->max[2]; z++) {
          int index = fluid_get_index(
              x - em->min[0], em->res[0], y - em->min[1], em->res[1], z - em->min[2]);
          float max_den = em->influence[index];

          /* density bounds */
          if (max_den >= sds->adapt_threshold) {
            if (min[0] > x) {
              min[0] = x;
            }
            if (min[1] > y) {
              min[1] = y;
            }
            if (min[2] > z) {
              min[2] = z;
            }
            if (max[0] < x) {
              max[0] = x;
            }
            if (max[1] < y) {
              max[1] = y;
            }
            if (max[2] < z) {
              max[2] = z;
            }
          }
        }
      }
    }
  }

  /* calculate new bounds based on these values */
  clampBoundsInDomain(sds, min, max, min_vel, max_vel, sds->adapt_margin + 1, dt);

  for (int i = 0; i < 3; i++) {
    /* calculate new resolution */
    res[i] = max[i] - min[i];
    total_cells *= res[i];

    if (new_shift[i]) {
      shift_changed = 1;
    }

    /* if no content set minimum dimensions */
    if (res[i] <= 0) {
      int j;
      for (j = 0; j < 3; j++) {
        min[j] = 0;
        max[j] = 1;
        res[j] = 1;
      }
      res_changed = 1;
      total_cells = 1;
      break;
    }
    if (min[i] != sds->res_min[i] || max[i] != sds->res_max[i]) {
      res_changed = 1;
    }
  }

  if (res_changed || shift_changed) {
    adaptiveDomainCopy(sds, sds->res, res, sds->res_min, min, sds->res_max, tmpShift, total_shift, 0);

    /* set new domain dimensions */
    copy_v3_v3_int(sds->res_min, min);
    copy_v3_v3_int(sds->res_max, max);
    copy_v3_v3_int(sds->res, res);
    sds->total_cells = total_cells;
  }
}

BLI_INLINE void apply_outflow_fields(int index,
                                     float distance_value,
                                     float *density,
                                     float *heat,
                                     float *fuel,
                                     float *react,
                                     float *color_r,
                                     float *color_g,
                                     float *color_b,
                                     float *phiout)
{
  /* determine outflow cells - phiout used in smoke and liquids */
  if (phiout) {
    phiout[index] = distance_value;
  }

  /* set smoke outflow */
  if (density) {
    density[index] = 0.f;
  }
  if (heat) {
    heat[index] = 0.f;
  }
  if (fuel) {
    fuel[index] = 0.f;
    react[index] = 0.f;
  }
  if (color_r) {
    color_r[index] = 0.f;
    color_g[index] = 0.f;
    color_b[index] = 0.f;
  }
}

BLI_INLINE void apply_inflow_fields(SmokeFlowSettings *sfs,
                                    float emission_value,
                                    float distance_value,
                                    int index,
                                    float *density,
                                    float *heat,
                                    float *fuel,
                                    float *react,
                                    float *color_r,
                                    float *color_g,
                                    float *color_b,
                                    float *phi,
                                    float *emission)
{
  /* add inflow */
  if (phi) {
    phi[index] = distance_value;
  }

  /* save emission value for manta inflow */
  if (emission) {
    emission[index] = emission_value;
  }

  /* add smoke inflow */
  int absolute_flow = (sfs->flags & FLUID_FLOW_ABSOLUTE);
  float dens_old = (density) ? density[index] : 0.0;
  // float fuel_old = (fuel) ? fuel[index] : 0.0f;  /* UNUSED */
  float dens_flow = (sfs->type == FLUID_FLOW_TYPE_FIRE) ? 0.0f : emission_value * sfs->density;
  float fuel_flow = (fuel) ? emission_value * sfs->fuel_amount : 0.0f;
  /* add heat */
  if (heat && emission_value > 0.0f) {
    heat[index] = ADD_IF_LOWER(heat[index], sfs->temp);
  }

  /* set density and fuel - absolute mode */
  if (absolute_flow) {
    if (density && sfs->type != FLUID_FLOW_TYPE_FIRE) {
      if (dens_flow > density[index]) {
        density[index] = dens_flow;
      }
    }
    if (fuel && sfs->type != FLUID_FLOW_TYPE_SMOKE && fuel_flow) {
      if (fuel_flow > fuel[index]) {
        fuel[index] = fuel_flow;
      }
    }
  }
  /* set density and fuel - additive mode */
  else {
    if (density && sfs->type != FLUID_FLOW_TYPE_FIRE) {
      density[index] += dens_flow;
      CLAMP(density[index], 0.0f, 1.0f);
    }
    if (fuel && sfs->type != FLUID_FLOW_TYPE_SMOKE && sfs->fuel_amount) {
      fuel[index] += fuel_flow;
      CLAMP(fuel[index], 0.0f, 10.0f);
    }
  }

  /* set color */
  if (color_r && dens_flow) {
    float total_dens = density[index] / (dens_old + dens_flow);
    color_r[index] = (color_r[index] + sfs->color[0] * dens_flow) * total_dens;
    color_g[index] = (color_g[index] + sfs->color[1] * dens_flow) * total_dens;
    color_b[index] = (color_b[index] + sfs->color[2] * dens_flow) * total_dens;
  }

  /* set fire reaction coordinate */
  if (fuel && fuel[index] > FLT_EPSILON) {
    /* instead of using 1.0 for all new fuel add slight falloff
     * to reduce flow blockiness */
    float value = 1.0f - pow2f(1.0f - emission_value);

    if (value > react[index]) {
      float f = fuel_flow / fuel[index];
      react[index] = value * f + (1.0f - f) * react[index];
      CLAMP(react[index], 0.0f, value);
    }
  }
}

static void update_flowsflags(SmokeDomainSettings *sds, Object **flowobjs, int numflowobj)
{
  int active_fields = sds->active_fields;
  unsigned int flowIndex;

  /* Monitor active fields based on flow settings */
  for (flowIndex = 0; flowIndex < numflowobj; flowIndex++) {
    Object *collob = flowobjs[flowIndex];
    SmokeModifierData *smd2 = (SmokeModifierData *)modifiers_findByType(collob,
                                                                        eModifierType_Smoke);

    // Sanity check
    if (!smd2) {
      continue;
    }

    if ((smd2->type & MOD_SMOKE_TYPE_FLOW) && smd2->flow) {
      SmokeFlowSettings *sfs = smd2->flow;
      if (!sfs) {
        break;
      }
      if (sfs->flags & FLUID_FLOW_INITVELOCITY) {
        active_fields |= FLUID_DOMAIN_ACTIVE_INVEL;
      }
      if (sfs->behavior == FLUID_FLOW_BEHAVIOR_OUTFLOW) {
        active_fields |= FLUID_DOMAIN_ACTIVE_OUTFLOW;
      }
      /* liquids done from here */
      if (sds->type == FLUID_DOMAIN_TYPE_LIQUID) {
        continue;
      }

      /* activate heat field if flow produces any heat */
      if (sfs->temp) {
        active_fields |= FLUID_DOMAIN_ACTIVE_HEAT;
      }
      /* activate fuel field if flow adds any fuel */
      if (sfs->fuel_amount &&
          (sfs->type == FLUID_FLOW_TYPE_FIRE || sfs->type == FLUID_FLOW_TYPE_SMOKEFIRE)) {
        active_fields |= FLUID_DOMAIN_ACTIVE_FIRE;
      }
      /* activate color field if flows add smoke with varying colors */
      if (sfs->density &&
          (sfs->type == FLUID_FLOW_TYPE_SMOKE || sfs->type == FLUID_FLOW_TYPE_SMOKEFIRE)) {
        if (!(active_fields & FLUID_DOMAIN_ACTIVE_COLOR_SET)) {
          copy_v3_v3(sds->active_color, sfs->color);
          active_fields |= FLUID_DOMAIN_ACTIVE_COLOR_SET;
        }
        else if (!equals_v3v3(sds->active_color, sfs->color)) {
          copy_v3_v3(sds->active_color, sfs->color);
          active_fields |= FLUID_DOMAIN_ACTIVE_COLORS;
        }
      }
    }
  }
  /* Monitor active fields based on domain settings */
  if (sds->type == FLUID_DOMAIN_TYPE_GAS && active_fields & FLUID_DOMAIN_ACTIVE_FIRE) {
    /* heat is always needed for fire */
    active_fields |= FLUID_DOMAIN_ACTIVE_HEAT;
    /* also activate colors if domain smoke color differs from active color */
    if (!(active_fields & FLUID_DOMAIN_ACTIVE_COLOR_SET)) {
      copy_v3_v3(sds->active_color, sds->flame_smoke_color);
      active_fields |= FLUID_DOMAIN_ACTIVE_COLOR_SET;
    }
    else if (!equals_v3v3(sds->active_color, sds->flame_smoke_color)) {
      copy_v3_v3(sds->active_color, sds->flame_smoke_color);
      active_fields |= FLUID_DOMAIN_ACTIVE_COLORS;
    }
  }
  /* Finally, initialize new data fields if any */
  if (active_fields & FLUID_DOMAIN_ACTIVE_INVEL) {
    fluid_ensure_invelocity(sds->fluid, sds->smd);
  }
  if (active_fields & FLUID_DOMAIN_ACTIVE_OUTFLOW) {
    fluid_ensure_outflow(sds->fluid, sds->smd);
  }
  if (active_fields & FLUID_DOMAIN_ACTIVE_HEAT) {
    smoke_ensure_heat(sds->fluid, sds->smd);
  }
  if (active_fields & FLUID_DOMAIN_ACTIVE_FIRE) {
    smoke_ensure_fire(sds->fluid, sds->smd);
  }
  if (active_fields & FLUID_DOMAIN_ACTIVE_COLORS) {
    /* initialize all smoke with "active_color" */
    smoke_ensure_colors(sds->fluid, sds->smd);
  }
  if (sds->type == FLUID_DOMAIN_TYPE_LIQUID &&
      (sds->particle_type & FLUID_DOMAIN_PARTICLE_SPRAY ||
       sds->particle_type & FLUID_DOMAIN_PARTICLE_FOAM ||
       sds->particle_type & FLUID_DOMAIN_PARTICLE_TRACER)) {
    liquid_ensure_sndparts(sds->fluid, sds->smd);
  }
  sds->active_fields = active_fields;
}

static void update_flowsfluids(struct Depsgraph *depsgraph,
                               Scene *scene,
                               Object *ob,
                               SmokeDomainSettings *sds,
                               float time_per_frame,
                               float frame_length,
                               int frame,
                               bool is_first_frame)
{
  EmissionMap *emaps = NULL;
  Object **flowobjs = NULL;
  unsigned int numflowobj = 0, flowIndex = 0;

  flowobjs = BKE_collision_objects_create(
      depsgraph, ob, sds->fluid_group, &numflowobj, eModifierType_Smoke);

  /* Update all flow related flags and ensure that corresponding grids get initialized */
  update_flowsflags(sds, flowobjs, numflowobj);

  /* init emission maps for each flow */
  emaps = MEM_callocN(sizeof(struct EmissionMap) * numflowobj, "smoke_flow_maps");

  /* Prepare flow emission maps */
  for (flowIndex = 0; flowIndex < numflowobj; flowIndex++) {
    Object *flowobj = flowobjs[flowIndex];
    SmokeModifierData *smd2 = (SmokeModifierData *)modifiers_findByType(flowobj,
                                                                        eModifierType_Smoke);

    /* Check for initialized smoke object */
    if ((smd2->type & MOD_SMOKE_TYPE_FLOW) && smd2->flow) {
      SmokeFlowSettings *sfs = smd2->flow;
      int subframes = sfs->subframes;
      EmissionMap *em = &emaps[flowIndex];

      /* Length of one frame. If using adaptive stepping, length is smaller than actual frame length */
      float adaptframe_length = time_per_frame / frame_length;

      /* Further splitting because of emission subframe: If no subframes present, sample_size is 1 */
      float sample_size = 1.0f / (float)(subframes + 1);
      float sdt = adaptframe_length * sample_size;
      int hires_multiplier = 1;

      /* First frame cannot have any subframes because there is (obviously) no previous frame from where subframes could come from */
      if (is_first_frame) {
        subframes = 0;
      }

      int subframe;
      float prev_frame_pos;

      /* Emission loop. When not using subframes this will loop only once. */
      for (subframe = 0; subframe <= subframes; subframe++) {

        /* Temporary emission map used when subframes are enabled, i.e. at least one subframe */
        EmissionMap em_temp = {NULL};

        /* Set scene time */
        /* Handle emission subframe */
        if (subframe < subframes && !is_first_frame) {
          prev_frame_pos = sdt * (float)(subframe + 1);
          scene->r.subframe = prev_frame_pos;
          scene->r.cfra = frame - 1;
        }
        /* Last frame in this loop (subframe == suframes). Can be real end frame or in between frames (adaptive frame) */
        else {
          /* Handle adaptive subframe (ie has subframe fraction). Need to set according scene subframe parameter */
          if (time_per_frame < frame_length) {
            scene->r.subframe = adaptframe_length;
            scene->r.cfra = frame - 1;
          }
          /* Handle absolute endframe (ie no subframe fraction). Need to set the scene subframe parameter to 0 and advance current scene frame */
          else {
            scene->r.subframe = 0.0f;
            scene->r.cfra = frame;
          }
        }
        //printf("flow: frame (is first: %d): %d // scene current frame: %d // scene current subframe: %f\n", is_first_frame, frame, scene->r.cfra, scene->r.subframe);

        /* Emission from particles */
        if (sfs->source == FLUID_FLOW_SOURCE_PARTICLES) {
          /* emit_from_particles() updates timestep internally */
          if (subframes) {
            emit_from_particles(flowobj, sds, sfs, &em_temp, depsgraph, scene, sdt);
          }
          else {
            emit_from_particles(flowobj, sds, sfs, em, depsgraph, scene, sdt);
          }

          if (!(sfs->flags & FLUID_FLOW_USE_PART_SIZE)) {
            hires_multiplier = 1;
          }
        }
        /* Emission from mesh */
        else if (sfs->source == FLUID_FLOW_SOURCE_MESH) {
          /* Update flow object frame */
          // BLI_mutex_lock() called in smoke_step(), so safe to update subframe here

          /* TODO (sebbas): Using BKE_scene_frame_get(scene) instead of new DEG_get_ctime(depsgraph) as subframes dont work with the latter yet */
          BKE_object_modifier_update_subframe(
              depsgraph, scene, flowobj, true, 5, BKE_scene_frame_get(scene), eModifierType_Smoke);

          /* Apply flow */
          if (subframes) {
            emit_from_mesh(flowobj, sds, sfs, &em_temp, sdt);
          }
          else {
            emit_from_mesh(flowobj, sds, sfs, em, sdt);
          }
        }
        else {
          printf("Error: unknown flow emission source\n");
        }

        /* If this we emitted with temp emission map in this loop (subframe emission), we combine the temp map with the original emission map */
        if (subframes) {
          /* Combine emission maps */
          em_combineMaps(
              em, &em_temp, hires_multiplier, !(sfs->flags & FLUID_FLOW_ABSOLUTE), sample_size);
          em_freeData(&em_temp);
        }
      }
    }
  }

  /* Adjust domain size if needed */
  if (sds->flags & FLUID_DOMAIN_USE_ADAPTIVE_DOMAIN) {
    adaptiveDomainAdjust(sds, ob, emaps, numflowobj, time_per_frame);
  }

  float *phi_in = fluid_get_phi_in(sds->fluid);
  float *phiout_in = fluid_get_phiout_in(sds->fluid);
  float *density = smoke_get_density(sds->fluid);
  float *color_r = smoke_get_color_r(sds->fluid);
  float *color_g = smoke_get_color_g(sds->fluid);
  float *color_b = smoke_get_color_b(sds->fluid);
  float *fuel = smoke_get_fuel(sds->fluid);
  float *heat = smoke_get_heat(sds->fluid);
  float *react = smoke_get_react(sds->fluid);

  float *density_in = smoke_get_density_in(sds->fluid);
  float *heat_in = smoke_get_heat_in(sds->fluid);
  float *color_r_in = smoke_get_color_r_in(sds->fluid);
  float *color_g_in = smoke_get_color_g_in(sds->fluid);
  float *color_b_in = smoke_get_color_b_in(sds->fluid);
  float *fuel_in = smoke_get_fuel_in(sds->fluid);
  float *react_in = smoke_get_react_in(sds->fluid);
  float *emission_in = smoke_get_emission_in(sds->fluid);

  float *velx_initial = fluid_get_in_velocity_x(sds->fluid);
  float *vely_initial = fluid_get_in_velocity_y(sds->fluid);
  float *velz_initial = fluid_get_in_velocity_z(sds->fluid);
  unsigned int z;

  /* Grid reset before writing again */
  for (z = 0; z < sds->res[0] * sds->res[1] * sds->res[2]; z++) {
    if (phi_in) {
      phi_in[z] = 9999;
    }
    if (phiout_in) {
      phiout_in[z] = 9999;
    }
  }

  /* Apply emission data */
  for (flowIndex = 0; flowIndex < numflowobj; flowIndex++) {
    Object *flowobj = flowobjs[flowIndex];
    SmokeModifierData *smd2 = (SmokeModifierData *)modifiers_findByType(flowobj,
                                                                        eModifierType_Smoke);

    // check for initialized smoke object
    if ((smd2->type & MOD_SMOKE_TYPE_FLOW) && smd2->flow) {
      SmokeFlowSettings *sfs = smd2->flow;
      EmissionMap *em = &emaps[flowIndex];
      float *velocity_map = em->velocity;
      float *emission_map = em->influence;
      float *distance_map = em->distances;

      int gx, gy, gz, ex, ey, ez, dx, dy, dz;
      size_t e_index, d_index;

      // loop through every emission map cell
      for (gx = em->min[0]; gx < em->max[0]; gx++) {
        for (gy = em->min[1]; gy < em->max[1]; gy++) {
          for (gz = em->min[2]; gz < em->max[2]; gz++) {
            /* get emission map index */
            ex = gx - em->min[0];
            ey = gy - em->min[1];
            ez = gz - em->min[2];
            e_index = fluid_get_index(ex, em->res[0], ey, em->res[1], ez);

            /* get domain index */
            dx = gx - sds->res_min[0];
            dy = gy - sds->res_min[1];
            dz = gz - sds->res_min[2];
            d_index = fluid_get_index(dx, sds->res[0], dy, sds->res[1], dz);
            /* make sure emission cell is inside the new domain boundary */
            if (dx < 0 || dy < 0 || dz < 0 || dx >= sds->res[0] || dy >= sds->res[1] ||
                dz >= sds->res[2]) {
              continue;
            }

            /* sync inflow grids with actual simulation grids, inflow computation needs information from actual simulation */
            if (density) {
              density_in[d_index] = density[d_index];
            }
            if (heat) {
              heat_in[d_index] = heat[d_index];
            }
            if (color_r) {
              color_r_in[d_index] = color_r[d_index];
              color_g_in[d_index] = color_g[d_index];
              color_b_in[d_index] = color_b[d_index];
            }
            if (fuel) {
              fuel_in[d_index] = fuel[d_index];
              react_in[d_index] = react[d_index];
            }

            if (sfs->behavior == FLUID_FLOW_BEHAVIOR_OUTFLOW) {  // outflow
              apply_outflow_fields(d_index,
                                   distance_map[e_index],
                                   density_in,
                                   heat_in,
                                   fuel_in,
                                   react_in,
                                   color_r_in,
                                   color_g_in,
                                   color_b_in,
                                   phiout_in);
            }
            else if (sfs->behavior == FLUID_FLOW_BEHAVIOR_GEOMETRY && smd2->time > 2) {
              apply_inflow_fields(sfs,
                                  0.0f,
                                  9999.0f,
                                  d_index,
                                  density_in,
                                  heat_in,
                                  fuel_in,
                                  react_in,
                                  color_r_in,
                                  color_g_in,
                                  color_b_in,
                                  phi_in,
                                  emission_in);
            }
            else if (sfs->behavior == FLUID_FLOW_BEHAVIOR_INFLOW ||
                     sfs->behavior == FLUID_FLOW_BEHAVIOR_GEOMETRY) {  // inflow
              /* only apply inflow if enabled */
              if (sfs->flags & FLUID_FLOW_USE_INFLOW) {
                apply_inflow_fields(sfs,
                                    emission_map[e_index],
                                    distance_map[e_index],
                                    d_index,
                                    density_in,
                                    heat_in,
                                    fuel_in,
                                    react_in,
                                    color_r_in,
                                    color_g_in,
                                    color_b_in,
                                    phi_in,
                                    emission_in);
                /* initial velocity */
                if (sfs->flags & FLUID_FLOW_INITVELOCITY) {
                  velx_initial[d_index] = velocity_map[e_index * 3];
                  vely_initial[d_index] = velocity_map[e_index * 3 + 1];
                  velz_initial[d_index] = velocity_map[e_index * 3 + 2];
                }
              }
            }
          }  // low res loop
        }
      }

      // free emission maps
      em_freeData(em);

    }  // end emission
  }

  BKE_collision_objects_free(flowobjs);
  if (emaps) {
    MEM_freeN(emaps);
  }
}

typedef struct UpdateEffectorsData {
  Scene *scene;
  SmokeDomainSettings *sds;
  ListBase *effectors;

  float *density;
  float *fuel;
  float *force_x;
  float *force_y;
  float *force_z;
  float *velocity_x;
  float *velocity_y;
  float *velocity_z;
  int *flags;
  float *phiObsIn;
} UpdateEffectorsData;

static void update_effectors_task_cb(void *__restrict userdata,
                                     const int x,
                                     const ParallelRangeTLS *__restrict UNUSED(tls))
{
  UpdateEffectorsData *data = userdata;
  SmokeDomainSettings *sds = data->sds;

  for (int y = 0; y < sds->res[1]; y++) {
    for (int z = 0; z < sds->res[2]; z++) {
      EffectedPoint epoint;
      float mag;
      float voxelCenter[3] = {0, 0, 0}, vel[3] = {0, 0, 0}, retvel[3] = {0, 0, 0};
      const unsigned int index = fluid_get_index(x, sds->res[0], y, sds->res[1], z);

      if ((data->fuel && MAX2(data->density[index], data->fuel[index]) < FLT_EPSILON) ||
          (data->density && data->density[index] < FLT_EPSILON) ||
          (data->phiObsIn && data->phiObsIn[index] < 0.0f) ||
          data->flags[index] & 2)  // mantaflow convention: 2 == FlagObstacle
      {
        continue;
      }

      /* get velocities from manta grid space and convert to blender units */
      vel[0] = data->velocity_x[index];
      vel[1] = data->velocity_y[index];
      vel[2] = data->velocity_z[index];
      mul_v3_fl(vel, sds->dx);

      /* convert vel to global space */
      mag = len_v3(vel);
      mul_mat3_m4_v3(sds->obmat, vel);
      normalize_v3(vel);
      mul_v3_fl(vel, mag);

      voxelCenter[0] = sds->p0[0] + sds->cell_size[0] * ((float)(x + sds->res_min[0]) + 0.5f);
      voxelCenter[1] = sds->p0[1] + sds->cell_size[1] * ((float)(y + sds->res_min[1]) + 0.5f);
      voxelCenter[2] = sds->p0[2] + sds->cell_size[2] * ((float)(z + sds->res_min[2]) + 0.5f);
      mul_m4_v3(sds->obmat, voxelCenter);

      /* do effectors */
      pd_point_from_loc(data->scene, voxelCenter, vel, index, &epoint);
      BKE_effectors_apply(data->effectors, NULL, sds->effector_weights, &epoint, retvel, NULL);

      /* convert retvel to local space */
      mag = len_v3(retvel);
      mul_mat3_m4_v3(sds->imat, retvel);
      normalize_v3(retvel);
      mul_v3_fl(retvel, mag);

      /* constrain forces to interval -1 to 1 */
      data->force_x[index] = min_ff(max_ff(-1.0f, retvel[0] * 0.2f), 1.0f);
      data->force_y[index] = min_ff(max_ff(-1.0f, retvel[1] * 0.2f), 1.0f);
      data->force_z[index] = min_ff(max_ff(-1.0f, retvel[2] * 0.2f), 1.0f);
    }
  }
}

static void update_effectors(
    Depsgraph *depsgraph, Scene *scene, Object *ob, SmokeDomainSettings *sds, float UNUSED(dt))
{
  ListBase *effectors;
  /* make sure smoke flow influence is 0.0f */
  sds->effector_weights->weight[PFIELD_SMOKEFLOW] = 0.0f;
  effectors = BKE_effectors_create(depsgraph, ob, NULL, sds->effector_weights);

  if (effectors) {
    // precalculate wind forces
    UpdateEffectorsData data;
    data.scene = scene;
    data.sds = sds;
    data.effectors = effectors;
    data.density = smoke_get_density(sds->fluid);
    data.fuel = smoke_get_fuel(sds->fluid);
    data.force_x = fluid_get_force_x(sds->fluid);
    data.force_y = fluid_get_force_y(sds->fluid);
    data.force_z = fluid_get_force_z(sds->fluid);
    data.velocity_x = fluid_get_velocity_x(sds->fluid);
    data.velocity_y = fluid_get_velocity_y(sds->fluid);
    data.velocity_z = fluid_get_velocity_z(sds->fluid);
    data.flags = smoke_get_obstacle(sds->fluid);
    data.phiObsIn = fluid_get_phiobs_in(sds->fluid);

    ParallelRangeSettings settings;
    BLI_parallel_range_settings_defaults(&settings);
    settings.scheduling_mode = TASK_SCHEDULING_DYNAMIC;
    BLI_task_parallel_range(0, sds->res[0], &data, update_effectors_task_cb, &settings);
  }

  BKE_effectors_free(effectors);
}

static Mesh *createLiquidGeometry(SmokeDomainSettings *sds, Mesh *orgmesh, Object *ob)
{
  Mesh *me;
  MVert *mverts;
  MPoly *mpolys;
  MLoop *mloops;
  short *normals, *no_s;
  float no[3];
  float min[3];
  float max[3];
  float size[3];
  float cell_size_scaled[3];

  /* assign material + flags to new dm
   * if there's no faces in original dm, keep materials and flags unchanged */
  MPoly *mpoly;
  MPoly mp_example = {0};
  mpoly = orgmesh->mpoly;
  if (mpoly) {
    mp_example = *mpoly;
  }
  /* else leave NULL'd */

  const short mp_mat_nr = mp_example.mat_nr;
  const char mp_flag = mp_example.flag;

  int i;
  int num_verts, num_normals, num_faces;

  if (!sds->fluid) {
    return NULL;
  }

  num_verts = liquid_get_num_verts(sds->fluid);
  num_normals = liquid_get_num_normals(sds->fluid);
  num_faces = liquid_get_num_triangles(sds->fluid);

  //printf("num_verts: %d, num_normals: %d, num_faces: %d\n", num_verts, num_normals, num_faces);

  if (!num_verts || !num_faces) {
    return NULL;
  }

  me = BKE_mesh_new_nomain(num_verts, 0, 0, num_faces * 3, num_faces);
  mverts = me->mvert;
  mpolys = me->mpoly;
  mloops = me->mloop;
  if (!me) {
    return NULL;
  }

  // Get size (dimension) but considering scaling scaling
  copy_v3_v3(cell_size_scaled, sds->cell_size);
  mul_v3_v3(cell_size_scaled, ob->scale);
  madd_v3fl_v3fl_v3fl_v3i(min, sds->p0, cell_size_scaled, sds->res_min);
  madd_v3fl_v3fl_v3fl_v3i(max, sds->p0, cell_size_scaled, sds->res_max);
  sub_v3_v3v3(size, max, min);

  // Biggest dimension will be used for upscaling
  float max_size = MAX3(size[0], size[1], size[2]);

  // Vertices
  for (i = 0; i < num_verts; i++, mverts++) {
    // read raw data. is normalized cube around domain origin
    mverts->co[0] = liquid_get_vertex_x_at(sds->fluid, i);
    mverts->co[1] = liquid_get_vertex_y_at(sds->fluid, i);
    mverts->co[2] = liquid_get_vertex_z_at(sds->fluid, i);

    // if reading raw data directly from manta, normalize now
    if ((sds->cache_flag & FLUID_DOMAIN_BAKED_MESH) == 0) {
      // normalize to unit cube around 0
      mverts->co[0] -= ((float)sds->res[0] * sds->mesh_scale) * 0.5f;
      mverts->co[1] -= ((float)sds->res[1] * sds->mesh_scale) * 0.5f;
      mverts->co[2] -= ((float)sds->res[2] * sds->mesh_scale) * 0.5f;
      mverts->co[0] *= sds->dx / sds->mesh_scale;
      mverts->co[1] *= sds->dx / sds->mesh_scale;
      mverts->co[2] *= sds->dx / sds->mesh_scale;
    }

    mverts->co[0] *= max_size / fabsf(ob->scale[0]);
    mverts->co[1] *= max_size / fabsf(ob->scale[1]);
    mverts->co[2] *= max_size / fabsf(ob->scale[2]);

    //printf("mverts->co[0]: %f, mverts->co[1]: %f, mverts->co[2]: %f\n", mverts->co[0], mverts->co[1], mverts->co[2]);
  }

  // Normals
  normals = MEM_callocN(sizeof(short) * num_normals * 3, "Fluidmesh_tmp_normals");

  for (i = 0, no_s = normals; i < num_normals; no_s += 3, i++) {
    no[0] = liquid_get_normal_x_at(sds->fluid, i);
    no[1] = liquid_get_normal_y_at(sds->fluid, i);
    no[2] = liquid_get_normal_z_at(sds->fluid, i);

    normal_float_to_short_v3(no_s, no);

    //printf("no_s[0]: %d, no_s[1]: %d, no_s[2]: %d\n", no_s[0], no_s[1], no_s[2]);
  }

  // Triangles
  for (i = 0; i < num_faces; i++, mpolys++, mloops += 3) {
    /* initialize from existing face */
    mpolys->mat_nr = mp_mat_nr;
    mpolys->flag = mp_flag;

    mpolys->loopstart = i * 3;
    mpolys->totloop = 3;

    mloops[0].v = liquid_get_triangle_x_at(sds->fluid, i);
    mloops[1].v = liquid_get_triangle_y_at(sds->fluid, i);
    mloops[2].v = liquid_get_triangle_z_at(sds->fluid, i);

    //printf("mloops[0].v: %d, mloops[1].v: %d, mloops[2].v: %d\n", mloops[0].v, mloops[1].v, mloops[2].v);
  }

  BKE_mesh_ensure_normals(me);
  BKE_mesh_calc_edges(me, false, false);
  BKE_mesh_apply_vert_normals(me, (short(*)[3])normals);

  MEM_freeN(normals);

  /* return early if no mesh vert velocities required */
  if ((sds->flags & FLUID_DOMAIN_USE_SPEED_VECTORS) == 0) {
    return me;
  }

  if (sds->mesh_velocities) {
    MEM_freeN(sds->mesh_velocities);
  }

  sds->mesh_velocities = MEM_calloc_arrayN(
      num_verts, sizeof(SmokeVertexVelocity), "Fluidmesh_vertvelocities");
  sds->totvert = num_verts;

  SmokeVertexVelocity *velarray = NULL;
  velarray = sds->mesh_velocities;

  float time_mult = 25.f * DT_DEFAULT;

  for (i = 0; i < num_verts; i++, mverts++) {
    velarray[i].vel[0] = liquid_get_vertvel_x_at(sds->fluid, i) * (sds->dx / time_mult);
    velarray[i].vel[1] = liquid_get_vertvel_y_at(sds->fluid, i) * (sds->dx / time_mult);
    velarray[i].vel[2] = liquid_get_vertvel_z_at(sds->fluid, i) * (sds->dx / time_mult);

    //printf("velarray[%d].vel[0]: %f, velarray[%d].vel[1]: %f, velarray[%d].vel[2]: %f\n", i, velarray[i].vel[0], i, velarray[i].vel[1], i, velarray[i].vel[2]);
  }

  return me;
}

static Mesh *createSmokeGeometry(SmokeDomainSettings *sds, Mesh *orgmesh, Object *ob)
{
  Mesh *result;
  MVert *mverts;
  MPoly *mpolys;
  MLoop *mloops;
  float min[3];
  float max[3];
  float *co;
  MPoly *mp;
  MLoop *ml;

  int num_verts = 8;
  int num_faces = 6;
  int i;
  float ob_loc[3] = {0};
  float ob_cache_loc[3] = {0};

  /* just copy existing mesh if there is no content or if the adaptive domain is not being used */
  if (sds->total_cells <= 1 || (sds-> flags & FLUID_DOMAIN_USE_ADAPTIVE_DOMAIN) == 0) {
    return BKE_mesh_copy_for_eval(orgmesh, false);
  }

  result = BKE_mesh_new_nomain(num_verts, 0, 0, num_faces * 4, num_faces);
  mverts = result->mvert;
  mpolys = result->mpoly;
  mloops = result->mloop;

  if (num_verts) {
    /* volume bounds */
    madd_v3fl_v3fl_v3fl_v3i(min, sds->p0, sds->cell_size, sds->res_min);
    madd_v3fl_v3fl_v3fl_v3i(max, sds->p0, sds->cell_size, sds->res_max);

    /* set vertices */
    /* top slab */
    co = mverts[0].co;
    co[0] = min[0];
    co[1] = min[1];
    co[2] = max[2];
    co = mverts[1].co;
    co[0] = max[0];
    co[1] = min[1];
    co[2] = max[2];
    co = mverts[2].co;
    co[0] = max[0];
    co[1] = max[1];
    co[2] = max[2];
    co = mverts[3].co;
    co[0] = min[0];
    co[1] = max[1];
    co[2] = max[2];
    /* bottom slab */
    co = mverts[4].co;
    co[0] = min[0];
    co[1] = min[1];
    co[2] = min[2];
    co = mverts[5].co;
    co[0] = max[0];
    co[1] = min[1];
    co[2] = min[2];
    co = mverts[6].co;
    co[0] = max[0];
    co[1] = max[1];
    co[2] = min[2];
    co = mverts[7].co;
    co[0] = min[0];
    co[1] = max[1];
    co[2] = min[2];

    /* create faces */
    /* top */
    mp = &mpolys[0];
    ml = &mloops[0 * 4];
    mp->loopstart = 0 * 4;
    mp->totloop = 4;
    ml[0].v = 0;
    ml[1].v = 1;
    ml[2].v = 2;
    ml[3].v = 3;
    /* right */
    mp = &mpolys[1];
    ml = &mloops[1 * 4];
    mp->loopstart = 1 * 4;
    mp->totloop = 4;
    ml[0].v = 2;
    ml[1].v = 1;
    ml[2].v = 5;
    ml[3].v = 6;
    /* bottom */
    mp = &mpolys[2];
    ml = &mloops[2 * 4];
    mp->loopstart = 2 * 4;
    mp->totloop = 4;
    ml[0].v = 7;
    ml[1].v = 6;
    ml[2].v = 5;
    ml[3].v = 4;
    /* left */
    mp = &mpolys[3];
    ml = &mloops[3 * 4];
    mp->loopstart = 3 * 4;
    mp->totloop = 4;
    ml[0].v = 0;
    ml[1].v = 3;
    ml[2].v = 7;
    ml[3].v = 4;
    /* front */
    mp = &mpolys[4];
    ml = &mloops[4 * 4];
    mp->loopstart = 4 * 4;
    mp->totloop = 4;
    ml[0].v = 3;
    ml[1].v = 2;
    ml[2].v = 6;
    ml[3].v = 7;
    /* back */
    mp = &mpolys[5];
    ml = &mloops[5 * 4];
    mp->loopstart = 5 * 4;
    mp->totloop = 4;
    ml[0].v = 1;
    ml[1].v = 0;
    ml[2].v = 4;
    ml[3].v = 5;

    /* calculate required shift to match domain's global position
     * it was originally simulated at (if object moves without smoke step) */
    invert_m4_m4(ob->imat, ob->obmat);
    mul_m4_v3(ob->obmat, ob_loc);
    mul_m4_v3(sds->obmat, ob_cache_loc);
    sub_v3_v3v3(sds->obj_shift_f, ob_cache_loc, ob_loc);
    /* convert shift to local space and apply to vertices */
    mul_mat3_m4_v3(ob->imat, sds->obj_shift_f);
    /* apply */
    for (i = 0; i < num_verts; i++) {
      add_v3_v3(mverts[i].co, sds->obj_shift_f);
    }
  }

  BKE_mesh_calc_edges(result, false, false);
  result->runtime.cd_dirty_vert |= CD_MASK_NORMAL;
  return result;
}

static void smoke_step(Depsgraph *depsgraph,
                       Scene *scene,
                       Object *ob,
                       Mesh *me,
                       SmokeModifierData *smd,
                       int frame,
                       bool is_first_frame)
{
  SmokeDomainSettings *sds = smd->domain;
  float fps = scene->r.frs_sec / scene->r.frs_sec_base;
  float dt;
  float sdt;  // dt after adapted timestep
  float time_per_frame;

  /* update object state */
  invert_m4_m4(sds->imat, ob->obmat);
  copy_m4_m4(sds->obmat, ob->obmat);
  smoke_set_domain_from_mesh(sds, ob, me, (sds->flags & FLUID_DOMAIN_USE_ADAPTIVE_DOMAIN) != 0);

  /* adapt timestep for different framerates, dt = 0.1 is at 25fps */
  dt = DT_DEFAULT * (25.0f / fps);

  time_per_frame = 0;

  BLI_mutex_lock(&object_update_lock);

  // loop as long as time_per_frame (sum of sudivdt) does not exceed dt (actual framelength)
  while (time_per_frame < dt) {
    fluid_adapt_timestep(sds->fluid);
    sdt = fluid_get_timestep(sds->fluid);
    time_per_frame += sdt;

    // Calculate inflow geometry
    update_flowsfluids(depsgraph, scene, ob, sds, time_per_frame, sdt, frame, is_first_frame);

    // Calculate obstacle geometry
    update_obstacles(depsgraph, scene, ob, sds, time_per_frame, dt, frame, sdt);

    if (sds->total_cells > 1) {
      update_effectors(
          depsgraph,
          scene,
          ob,
          sds,
          sdt);  // DG TODO? problem --> uses forces instead of velocity, need to check how they need to be changed with variable dt
      fluid_bake_data(sds->fluid, smd, frame);
    }
  }
  if (sds->type == FLUID_DOMAIN_TYPE_GAS) {
    smoke_calc_transparency(sds, DEG_get_evaluated_view_layer(depsgraph));
  }
  BLI_mutex_unlock(&object_update_lock);
}

static void smoke_guiding(
    Depsgraph *depsgraph, Scene *scene, Object *ob, SmokeModifierData *smd, int frame)
{
  SmokeDomainSettings *sds = smd->domain;
  float fps = scene->r.frs_sec / scene->r.frs_sec_base;
  float dt = DT_DEFAULT * (25.0f / fps);

  BLI_mutex_lock(&object_update_lock);

  update_obstacles(depsgraph, scene, ob, sds, dt, dt, frame, dt);
  fluid_bake_guiding(sds->fluid, smd, frame);

  BLI_mutex_unlock(&object_update_lock);
}

static void smokeModifier_process(
    SmokeModifierData *smd, Depsgraph *depsgraph, Scene *scene, Object *ob, Mesh *me)
{
  const int scene_framenr = (int)DEG_get_ctime(depsgraph);

  if ((smd->type & MOD_SMOKE_TYPE_FLOW)) {
    if (scene_framenr >= smd->time) {
      smokeModifier_init(smd, depsgraph, ob, scene, me);
    }

    if (smd->flow) {
      if (smd->flow->mesh) {
        BKE_id_free(NULL, smd->flow->mesh);
      }
      smd->flow->mesh = BKE_mesh_copy_for_eval(me, false);
    }

    if (scene_framenr > smd->time) {
      smd->time = scene_framenr;
    }
    else if (scene_framenr < smd->time) {
      smd->time = scene_framenr;
      smokeModifier_reset_ex(smd, false);
    }
  }
  else if (smd->type & MOD_SMOKE_TYPE_EFFEC) {
    if (scene_framenr >= smd->time) {
      smokeModifier_init(smd, depsgraph, ob, scene, me);
    }

    if (smd->effec) {
      if (smd->effec->mesh) {
        BKE_id_free(NULL, smd->effec->mesh);
      }
      smd->effec->mesh = BKE_mesh_copy_for_eval(me, false);
    }

    if (scene_framenr > smd->time) {
      smd->time = scene_framenr;
    }
    else if (scene_framenr < smd->time) {
      smd->time = scene_framenr;
      smokeModifier_reset_ex(smd, false);
    }
  }
  else if (smd->type & MOD_SMOKE_TYPE_DOMAIN) {
    SmokeDomainSettings *sds = smd->domain;
    int startframe, endframe;
    Object *guiding_parent = NULL;
    Object **objs = NULL;
    unsigned int numobj = 0;
    SmokeModifierData *smd_parent = NULL;
    bool is_first_frame;
    startframe = sds->cache_frame_start;
    endframe = sds->cache_frame_end;

    is_first_frame = (scene_framenr == startframe);

    bool is_baking = (sds->cache_flag & (FLUID_DOMAIN_BAKING_DATA | FLUID_DOMAIN_BAKING_NOISE |
                                         FLUID_DOMAIN_BAKING_MESH | FLUID_DOMAIN_BAKING_PARTICLES |
                                         FLUID_DOMAIN_BAKING_GUIDING));
    bool is_baked = (sds->cache_flag & (FLUID_DOMAIN_BAKED_DATA | FLUID_DOMAIN_BAKED_NOISE |
                                        FLUID_DOMAIN_BAKED_MESH | FLUID_DOMAIN_BAKED_PARTICLES |
                                        FLUID_DOMAIN_BAKED_GUIDING));

    /* Reset fluid if no fluid present (obviously)
     * or if timeline gets reset to startframe when no (!) baking is running
     * or if no baking is running and also there is no baked data present */
    if (!sds->fluid || (scene_framenr == startframe && !is_baking) || (!is_baking && !is_baked)) {
      smokeModifier_reset_ex(smd, false);
    }

    if (smokeModifier_init(smd, depsgraph, ob, scene, me) == 0) {
      return;
    }

    /* Guiding parent res pointer needs initialization */
    guiding_parent = sds->guiding_parent;
    if (guiding_parent) {
      smd_parent = (SmokeModifierData *)modifiers_findByType(guiding_parent, eModifierType_Smoke);
      if (smd_parent->domain) {
        sds->guide_res = smd_parent->domain->res;
      }
    }

    /* Cache does not keep track of active fields yet. So refresh them here */
    objs = BKE_collision_objects_create(
        depsgraph, ob, sds->fluid_group, &numobj, eModifierType_Smoke);
    update_flowsflags(sds, objs, numobj);
    if (objs) {
      MEM_freeN(objs);
    }

    objs = BKE_collision_objects_create(
        depsgraph, ob, sds->coll_group, &numobj, eModifierType_Smoke);
    update_obstacleflags(sds, objs, numobj);
    if (objs) {
      MEM_freeN(objs);
    }

    /* Read cache. For liquids update data directly (i.e. not via python) */
    if (!is_baking) {
      if (sds->cache_flag & FLUID_DOMAIN_BAKED_DATA)
      {
        if (sds->type == FLUID_DOMAIN_TYPE_GAS) {
          if (fluid_read_config(sds->fluid, smd, scene_framenr)) {
            /* Adaptive domain might have changed resolution */
            if (fluid_needs_realloc(sds->fluid, smd)) {
              BKE_smoke_reallocate_fluid(sds, sds->res, 1);
            }
            fluid_read_data(sds->fluid, smd, scene_framenr);
          }
        }
        if (sds->type == FLUID_DOMAIN_TYPE_LIQUID) {
          fluid_update_liquid_structures(sds->fluid, smd, scene_framenr);
        }
      }
      if (sds->cache_flag & FLUID_DOMAIN_BAKED_NOISE)
      {
        if (fluid_read_config(sds->fluid, smd, scene_framenr)) {
          if (fluid_needs_realloc(sds->fluid, smd)) {
            BKE_smoke_reallocate_fluid(sds, sds->res, 1);
          }
          fluid_read_noise(sds->fluid, smd, scene_framenr);
        }
      }
      if (sds->cache_flag & FLUID_DOMAIN_BAKED_MESH)
      {
        //if (sds->type == FLUID_DOMAIN_TYPE_GAS)
        // TODO (sebbas): smoke as mesh
        if (sds->type == FLUID_DOMAIN_TYPE_LIQUID) {
          fluid_update_mesh_structures(sds->fluid, smd, scene_framenr);
        }
      }
      if (sds->cache_flag & FLUID_DOMAIN_BAKED_PARTICLES)
      {
        //if (sds->type == FLUID_DOMAIN_TYPE_GAS)
        // TODO (sebbas): fire particles
        if (sds->type == FLUID_DOMAIN_TYPE_LIQUID) {
          fluid_update_particle_structures(sds->fluid, smd, scene_framenr);
        }
      }
    }

    /* Simulate step and write cache. Optionally also read py objects once from previous frame (bake started from resume operator) */
    if (is_baking) {
      /* Ensure fresh variables at every animation step */
      fluid_update_variables(sds->fluid, smd);

      /* Export mantaflow python script on first frame (once only) and for any bake type */
      if ((sds->flags & FLUID_DOMAIN_EXPORT_MANTA_SCRIPT) &&
          scene_framenr == sds->cache_frame_start)
      {
        if (smd->domain && smd->domain->type == FLUID_DOMAIN_TYPE_GAS) {
          smoke_manta_export(smd->domain->fluid, smd);
        }
        if (smd->domain && smd->domain->type == FLUID_DOMAIN_TYPE_LIQUID) {
          liquid_manta_export(smd->domain->fluid, smd);
        }
      }

      if (sds->cache_flag & FLUID_DOMAIN_BAKING_DATA)
      {
        if (sds->flags & FLUID_DOMAIN_USE_GUIDING) {
          /* Load guiding vel from flow object (only if baked) or else from domain object */
          if (sds->guiding_source == FLUID_DOMAIN_GUIDING_SRC_EFFECTOR &&
              sds->cache_flag & FLUID_DOMAIN_BAKED_GUIDING) {
            fluid_read_guiding(sds->fluid, smd, scene_framenr, false);
          }
          else if (sds->guiding_source == FLUID_DOMAIN_GUIDING_SRC_DOMAIN && smd_parent) {
            fluid_read_guiding(sds->fluid, smd_parent, scene_framenr, true);
          }
        }

        /* Refresh all objects if we start baking from a resumed frame */
        if (sds->cache_frame_start != scene_framenr &&
            sds->cache_frame_pause_data == scene_framenr) {
          fluid_read_data(sds->fluid, smd, scene_framenr - 1);
        }

        /* Base step needs separated bake and write calls - reason being that transparency calculation is after fluid step */
        smoke_step(depsgraph, scene, ob, me, smd, scene_framenr, is_first_frame);
        fluid_write_config(sds->fluid, smd, scene_framenr);
        fluid_write_data(sds->fluid, smd, scene_framenr);
      }
      if (sds->cache_flag & FLUID_DOMAIN_BAKING_NOISE)
      {
        /* Refresh all objects if we start baking from a resumed frame */
        if (sds->cache_frame_start != scene_framenr &&
            sds->cache_frame_pause_data == scene_framenr) {
          fluid_read_noise(sds->fluid, smd, scene_framenr - 1);
        }
        if (sds->flags & FLUID_DOMAIN_USE_ADAPTIVE_DOMAIN) {
          int o_res[3], o_min[3], o_max[3], o_shift[3];
          copy_v3_v3_int(o_res, sds->res);
          copy_v3_v3_int(o_min, sds->res_min);
          copy_v3_v3_int(o_max, sds->res_max);
          copy_v3_v3_int(o_shift, sds->shift);
          fluid_read_config(sds->fluid, smd, scene_framenr);
          if (fluid_needs_realloc(sds->fluid, smd)) {
            adaptiveDomainCopy(sds, o_res, sds->res, o_min, sds->res_min, o_max, o_shift, sds->shift, 1);
          }
        }
        fluid_bake_noise(sds->fluid, smd, scene_framenr);
      }
      if (sds->cache_flag & FLUID_DOMAIN_BAKING_MESH)
      {
        /* Note: Mesh bake does not need object refresh from cache */
        fluid_bake_mesh(sds->fluid, smd, scene_framenr);
      }
      if (sds->cache_flag & FLUID_DOMAIN_BAKING_PARTICLES)
      {
        /* Refresh all objects if we start baking from a resumed frame */
        if (sds->cache_frame_start != scene_framenr &&
            sds->cache_frame_pause_data == scene_framenr) {
          fluid_read_particles(sds->fluid, smd, scene_framenr - 1);
        }
        fluid_bake_particles(sds->fluid, smd, scene_framenr);
      }
      if (sds->cache_flag & FLUID_DOMAIN_BAKING_GUIDING)
      {
        smoke_guiding(depsgraph, scene, ob, smd, scene_framenr);
      }
    }
    smd->time = scene_framenr;
  }
}

struct Mesh *smokeModifier_do(
    SmokeModifierData *smd, Depsgraph *depsgraph, Scene *scene, Object *ob, Mesh *me)
{
  /* lock so preview render does not read smoke data while it gets modified */
  if ((smd->type & MOD_SMOKE_TYPE_DOMAIN) && smd->domain) {
    BLI_rw_mutex_lock(smd->domain->fluid_mutex, THREAD_LOCK_WRITE);
  }

  smokeModifier_process(smd, depsgraph, scene, ob, me);

  if ((smd->type & MOD_SMOKE_TYPE_DOMAIN) && smd->domain) {
    BLI_rw_mutex_unlock(smd->domain->fluid_mutex);
  }

  /* return generated geometry for adaptive domain */
  Mesh *result = NULL;
  if (smd->type & MOD_SMOKE_TYPE_DOMAIN && smd->domain)
  {
    if (smd->domain->type == FLUID_DOMAIN_TYPE_LIQUID) {
      result = createLiquidGeometry(smd->domain, me, ob);
    }
    if (smd->domain->type == FLUID_DOMAIN_TYPE_GAS) {
      result = createSmokeGeometry(smd->domain, me, ob);
    }
  }
  if (!result) {
    result = BKE_mesh_copy_for_eval(me, false);
  }
  /* XXX This is really not a nice hack, but until root of the problem is understood,
   * this should be an acceptable workaround I think.
   * See T58492 for details on the issue. */
  result->texflag |= ME_AUTOSPACE;
  return result;
}

static float calc_voxel_transp(
    float *result, float *input, int res[3], int *pixel, float *tRay, float correct)
{
  const size_t index = fluid_get_index(pixel[0], res[0], pixel[1], res[1], pixel[2]);

  // T_ray *= T_vox
  *tRay *= expf(input[index] * correct);

  if (result[index] < 0.0f) {
    result[index] = *tRay;
  }

  return *tRay;
}

static void bresenham_linie_3D(int x1,
                               int y1,
                               int z1,
                               int x2,
                               int y2,
                               int z2,
                               float *tRay,
                               bresenham_callback cb,
                               float *result,
                               float *input,
                               int res[3],
                               float correct)
{
  int dx, dy, dz, i, l, m, n, x_inc, y_inc, z_inc, err_1, err_2, dx2, dy2, dz2;
  int pixel[3];

  pixel[0] = x1;
  pixel[1] = y1;
  pixel[2] = z1;

  dx = x2 - x1;
  dy = y2 - y1;
  dz = z2 - z1;

  x_inc = (dx < 0) ? -1 : 1;
  l = abs(dx);
  y_inc = (dy < 0) ? -1 : 1;
  m = abs(dy);
  z_inc = (dz < 0) ? -1 : 1;
  n = abs(dz);
  dx2 = l << 1;
  dy2 = m << 1;
  dz2 = n << 1;

  if ((l >= m) && (l >= n)) {
    err_1 = dy2 - l;
    err_2 = dz2 - l;
    for (i = 0; i < l; i++) {
      if (cb(result, input, res, pixel, tRay, correct) <= FLT_EPSILON) {
        break;
      }
      if (err_1 > 0) {
        pixel[1] += y_inc;
        err_1 -= dx2;
      }
      if (err_2 > 0) {
        pixel[2] += z_inc;
        err_2 -= dx2;
      }
      err_1 += dy2;
      err_2 += dz2;
      pixel[0] += x_inc;
    }
  }
  else if ((m >= l) && (m >= n)) {
    err_1 = dx2 - m;
    err_2 = dz2 - m;
    for (i = 0; i < m; i++) {
      if (cb(result, input, res, pixel, tRay, correct) <= FLT_EPSILON) {
        break;
      }
      if (err_1 > 0) {
        pixel[0] += x_inc;
        err_1 -= dy2;
      }
      if (err_2 > 0) {
        pixel[2] += z_inc;
        err_2 -= dy2;
      }
      err_1 += dx2;
      err_2 += dz2;
      pixel[1] += y_inc;
    }
  }
  else {
    err_1 = dy2 - n;
    err_2 = dx2 - n;
    for (i = 0; i < n; i++) {
      if (cb(result, input, res, pixel, tRay, correct) <= FLT_EPSILON) {
        break;
      }
      if (err_1 > 0) {
        pixel[1] += y_inc;
        err_1 -= dz2;
      }
      if (err_2 > 0) {
        pixel[0] += x_inc;
        err_2 -= dz2;
      }
      err_1 += dy2;
      err_2 += dx2;
      pixel[2] += z_inc;
    }
  }
  cb(result, input, res, pixel, tRay, correct);
}

static void smoke_calc_transparency(SmokeDomainSettings *sds, ViewLayer *view_layer)
{
  float bv[6] = {0};
  float light[3];
  int a, z, slabsize = sds->res[0] * sds->res[1], size = sds->res[0] * sds->res[1] * sds->res[2];
  float *density = smoke_get_density(sds->fluid);
  float *shadow = smoke_get_shadow(sds->fluid);
  float correct = -7.0f * sds->dx;

  if (!get_light(view_layer, light)) {
    return;
  }

  /* convert light pos to sim cell space */
  mul_m4_v3(sds->imat, light);
  light[0] = (light[0] - sds->p0[0]) / sds->cell_size[0] - 0.5f - (float)sds->res_min[0];
  light[1] = (light[1] - sds->p0[1]) / sds->cell_size[1] - 0.5f - (float)sds->res_min[1];
  light[2] = (light[2] - sds->p0[2]) / sds->cell_size[2] - 0.5f - (float)sds->res_min[2];

  for (a = 0; a < size; a++) {
    shadow[a] = -1.0f;
  }

  /* calculate domain bounds in sim cell space */
  // 0,2,4 = 0.0f
  bv[1] = (float)sds->res[0];  // x
  bv[3] = (float)sds->res[1];  // y
  bv[5] = (float)sds->res[2];  // z

  for (z = 0; z < sds->res[2]; z++) {
    size_t index = z * slabsize;
    int x, y;

    for (y = 0; y < sds->res[1]; y++) {
      for (x = 0; x < sds->res[0]; x++, index++) {
        float voxelCenter[3];
        float pos[3];
        int cell[3];
        float tRay = 1.0;

        if (shadow[index] >= 0.0f) {
          continue;
        }
        voxelCenter[0] = (float)x;
        voxelCenter[1] = (float)y;
        voxelCenter[2] = (float)z;

        // get starting cell (light pos)
        if (BLI_bvhtree_bb_raycast(bv, light, voxelCenter, pos) > FLT_EPSILON) {
          // we're ouside -> use point on side of domain
          cell[0] = (int)floor(pos[0]);
          cell[1] = (int)floor(pos[1]);
          cell[2] = (int)floor(pos[2]);
        }
        else {
          // we're inside -> use light itself
          cell[0] = (int)floor(light[0]);
          cell[1] = (int)floor(light[1]);
          cell[2] = (int)floor(light[2]);
        }
        /* clamp within grid bounds */
        CLAMP(cell[0], 0, sds->res[0] - 1);
        CLAMP(cell[1], 0, sds->res[1] - 1);
        CLAMP(cell[2], 0, sds->res[2] - 1);

        bresenham_linie_3D(cell[0],
                           cell[1],
                           cell[2],
                           x,
                           y,
                           z,
                           &tRay,
                           calc_voxel_transp,
                           shadow,
                           density,
                           sds->res,
                           correct);

        // convention -> from a RGBA float array, use G value for tRay
        shadow[index] = tRay;
      }
    }
  }
}

/* get smoke velocity and density at given coordinates
 * returns fluid density or -1.0f if outside domain. */
float BKE_smoke_get_velocity_at(struct Object *ob, float position[3], float velocity[3])
{
  SmokeModifierData *smd = (SmokeModifierData *)modifiers_findByType(ob, eModifierType_Smoke);
  zero_v3(velocity);

  if (smd && (smd->type & MOD_SMOKE_TYPE_DOMAIN) && smd->domain && smd->domain->fluid) {
    SmokeDomainSettings *sds = smd->domain;
    float time_mult = 25.f * DT_DEFAULT;
    float vel_mag;
    float *velX = fluid_get_velocity_x(sds->fluid);
    float *velY = fluid_get_velocity_y(sds->fluid);
    float *velZ = fluid_get_velocity_z(sds->fluid);
    float density = 0.0f, fuel = 0.0f;
    float pos[3];
    copy_v3_v3(pos, position);
    smoke_pos_to_cell(sds, pos);

    /* check if point is outside domain max bounds */
    if (pos[0] < sds->res_min[0] || pos[1] < sds->res_min[1] || pos[2] < sds->res_min[2]) {
      return -1.0f;
    }
    if (pos[0] > sds->res_max[0] || pos[1] > sds->res_max[1] || pos[2] > sds->res_max[2]) {
      return -1.0f;
    }

    /* map pos between 0.0 - 1.0 */
    pos[0] = (pos[0] - sds->res_min[0]) / ((float)sds->res[0]);
    pos[1] = (pos[1] - sds->res_min[1]) / ((float)sds->res[1]);
    pos[2] = (pos[2] - sds->res_min[2]) / ((float)sds->res[2]);

    /* check if point is outside active area */
    if (smd->domain->flags & FLUID_DOMAIN_USE_ADAPTIVE_DOMAIN) {
      if (pos[0] < 0.0f || pos[1] < 0.0f || pos[2] < 0.0f) {
        return 0.0f;
      }
      if (pos[0] > 1.0f || pos[1] > 1.0f || pos[2] > 1.0f) {
        return 0.0f;
      }
    }

    /* get interpolated velocity */
    velocity[0] = BLI_voxel_sample_trilinear(velX, sds->res, pos) * sds->global_size[0] *
                  time_mult;
    velocity[1] = BLI_voxel_sample_trilinear(velY, sds->res, pos) * sds->global_size[1] *
                  time_mult;
    velocity[2] = BLI_voxel_sample_trilinear(velZ, sds->res, pos) * sds->global_size[2] *
                  time_mult;

    /* convert velocity direction to global space */
    vel_mag = len_v3(velocity);
    mul_mat3_m4_v3(sds->obmat, velocity);
    normalize_v3(velocity);
    mul_v3_fl(velocity, vel_mag);

    /* use max value of fuel or smoke density */
    density = BLI_voxel_sample_trilinear(smoke_get_density(sds->fluid), sds->res, pos);
    if (smoke_has_fuel(sds->fluid)) {
      fuel = BLI_voxel_sample_trilinear(smoke_get_fuel(sds->fluid), sds->res, pos);
    }
    return MAX2(density, fuel);
  }
  return -1.0f;
}

int BKE_smoke_get_data_flags(SmokeDomainSettings *sds)
{
  int flags = 0;

  if (sds->fluid) {
    if (smoke_has_heat(sds->fluid)) {
      flags |= FLUID_DOMAIN_ACTIVE_HEAT;
    }
    if (smoke_has_fuel(sds->fluid)) {
      flags |= FLUID_DOMAIN_ACTIVE_FIRE;
    }
    if (smoke_has_colors(sds->fluid)) {
      flags |= FLUID_DOMAIN_ACTIVE_COLORS;
    }
  }

  return flags;
}

#endif /* WITH_MANTA */
