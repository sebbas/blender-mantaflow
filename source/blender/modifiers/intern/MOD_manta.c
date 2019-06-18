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
 * along with this program; if not, write to the Free Software  Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * The Original Code is Copyright (C) 2005 by the Blender Foundation.
 * All rights reserved.
 */

/** \file
 * \ingroup modifiers
 */

#include <stddef.h>

#include "MEM_guardedalloc.h"

#include "BLI_utildefines.h"

#include "DNA_collection_types.h"
#include "DNA_object_types.h"
#include "DNA_scene_types.h"
#include "DNA_manta_types.h"
#include "DNA_object_force_types.h"
#include "DNA_mesh_types.h"

#include "BKE_cdderivedmesh.h"
#include "BKE_layer.h"
#include "BKE_library_query.h"
#include "BKE_modifier.h"
#include "BKE_manta.h"

#include "DEG_depsgraph.h"
#include "DEG_depsgraph_build.h"
#include "DEG_depsgraph_physics.h"
#include "DEG_depsgraph_query.h"

#include "MOD_modifiertypes.h"

static void initData(ModifierData *md)
{
  MantaModifierData *mmd = (MantaModifierData *)md;

  mmd->domain = NULL;
  mmd->flow = NULL;
  mmd->effec = NULL;
  mmd->type = 0;
  mmd->time = -1;
}

static void copyData(const ModifierData *md, ModifierData *target, const int flag)
{
  const MantaModifierData *mmd = (const MantaModifierData *)md;
  MantaModifierData *tmmd = (MantaModifierData *)target;

  mantaModifier_free(tmmd);
  mantaModifier_copy(mmd, tmmd, flag);
}

static void freeData(ModifierData *md)
{
  MantaModifierData *mmd = (MantaModifierData *)md;

  mantaModifier_free(mmd);
}

static void requiredDataMask(Object *UNUSED(ob),
                             ModifierData *md,
                             CustomData_MeshMasks *r_cddata_masks)
{
  MantaModifierData *mmd = (MantaModifierData *)md;

  if (mmd && (mmd->type & MOD_MANTA_TYPE_FLOW) && mmd->flow) {
    if (mmd->flow->source == FLUID_FLOW_SOURCE_MESH) {
      /* vertex groups */
      if (mmd->flow->vgroup_density) {
        r_cddata_masks->vmask |= CD_MASK_MDEFORMVERT;
      }
      /* uv layer */
      if (mmd->flow->texture_type == FLUID_FLOW_TEXTURE_MAP_UV) {
        r_cddata_masks->fmask |= CD_MASK_MTFACE;
      }
    }
  }
}

static Mesh *applyModifier(ModifierData *md, const ModifierEvalContext *ctx, Mesh *me)
{
  MantaModifierData *mmd = (MantaModifierData *)md;
  Mesh *result = NULL;

  if (ctx->flag & MOD_APPLY_ORCO) {
    return me;
  }

  Scene *scene = DEG_get_evaluated_scene(ctx->depsgraph);

  result = mantaModifier_do(mmd, ctx->depsgraph, scene, ctx->object, me);
  return result ? result : me;
}

static bool dependsOnTime(ModifierData *UNUSED(md))
{
  return true;
}

static bool is_flow_cb(Object *UNUSED(ob), ModifierData *md)
{
  MantaModifierData *mmd = (MantaModifierData *)md;
  return (mmd->type & MOD_MANTA_TYPE_FLOW) && mmd->flow;
}

static bool is_coll_cb(Object *UNUSED(ob), ModifierData *md)
{
  MantaModifierData *mmd = (MantaModifierData *)md;
  return (mmd->type & MOD_MANTA_TYPE_EFFEC) && mmd->effec;
}

static void updateDepsgraph(ModifierData *md, const ModifierUpdateDepsgraphContext *ctx)
{
  MantaModifierData *mmd = (MantaModifierData *)md;

  if (mmd && (mmd->type & MOD_MANTA_TYPE_DOMAIN) && mmd->domain) {
    DEG_add_collision_relations(ctx->node,
                                ctx->object,
                                mmd->domain->fluid_group,
                                eModifierType_Manta,
                                is_flow_cb,
                                "Manta Flow");
    DEG_add_collision_relations(ctx->node,
                                ctx->object,
                                mmd->domain->coll_group,
                                eModifierType_Manta,
                                is_coll_cb,
                                "Manta Coll");
    DEG_add_forcefield_relations(ctx->node,
                                 ctx->object,
                                 mmd->domain->effector_weights,
                                 true,
                                 PFIELD_SMOKEFLOW,
                                 "Manta Force Field");

    if (mmd->domain->guiding_parent != NULL) {
      DEG_add_object_relation(
          ctx->node, mmd->domain->guiding_parent, DEG_OB_COMP_TRANSFORM, "Fluid Guiding Object");
      DEG_add_object_relation(
          ctx->node, mmd->domain->guiding_parent, DEG_OB_COMP_GEOMETRY, "Fluid Guiding Object");
    }
  }
}

static void foreachIDLink(ModifierData *md, Object *ob, IDWalkFunc walk, void *userData)
{
  MantaModifierData *mmd = (MantaModifierData *)md;

  if (mmd->type == MOD_MANTA_TYPE_DOMAIN && mmd->domain) {
    walk(userData, ob, (ID **)&mmd->domain->coll_group, IDWALK_CB_NOP);
    walk(userData, ob, (ID **)&mmd->domain->fluid_group, IDWALK_CB_NOP);
    walk(userData, ob, (ID **)&mmd->domain->eff_group, IDWALK_CB_NOP);

    if (mmd->domain->guiding_parent) {
      walk(userData, ob, (ID **)&mmd->domain->guiding_parent, IDWALK_CB_NOP);
    }

    if (mmd->domain->effector_weights) {
      walk(userData, ob, (ID **)&mmd->domain->effector_weights->group, IDWALK_CB_NOP);
    }
  }

  if (mmd->type == MOD_MANTA_TYPE_FLOW && mmd->flow) {
    walk(userData, ob, (ID **)&mmd->flow->noise_texture, IDWALK_CB_USER);
  }
}

ModifierTypeInfo modifierType_Manta = {
    /* name */ "Manta",
    /* structName */ "MantaModifierData",
    /* structSize */ sizeof(MantaModifierData),
    /* type */ eModifierTypeType_Constructive,
    /* flags */ eModifierTypeFlag_AcceptsMesh | eModifierTypeFlag_UsesPointCache |
        eModifierTypeFlag_Single,

    /* copyData */ copyData,

    /* deformVerts */ NULL,
    /* deformMatrices */ NULL,
    /* deformVertsEM */ NULL,
    /* deformMatricesEM */ NULL,
    /* applyModifier */ applyModifier,

    /* initData */ initData,
    /* requiredDataMask */ requiredDataMask,
    /* freeData */ freeData,
    /* isDisabled */ NULL,
    /* updateDepsgraph */ updateDepsgraph,
    /* dependsOnTime */ dependsOnTime,
    /* dependsOnNormals */ NULL,
    /* foreachObjectLink */ NULL,
    /* foreachIDLink */ foreachIDLink,
    /* foreachTexLink */ NULL,
    /* freeRuntimeData */ NULL,
};
