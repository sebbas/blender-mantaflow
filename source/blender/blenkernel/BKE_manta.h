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

#ifndef __BKE_MANTA_H__
#define __BKE_MANTA_H__

/** \file
 * \ingroup bke
 */

struct Scene;
struct MantaDomainSettings;
struct MantaModifierData;

typedef float (*bresenham_callback)(
    float *result, float *input, int res[3], int *pixel, float *tRay, float correct);

struct Mesh *mantaModifier_do(struct MantaModifierData *mmd,
                              struct Depsgraph *depsgraph,
                              struct Scene *scene,
                              struct Object *ob,
                              struct Mesh *me);

void mantaModifier_free(struct MantaModifierData *mmd);
void mantaModifier_reset(struct MantaModifierData *mmd);
void mantaModifier_createType(struct MantaModifierData *mmd);
void mantaModifier_copy(const struct MantaModifierData *mmd,
                        struct MantaModifierData *tmmd,
                        const int flag);

void BKE_manta_reallocate_fluid(struct MantaDomainSettings *mds, int res[3], int free_old);
void BKE_manta_reallocate_copy_fluid(struct MantaDomainSettings *mds,
                                     int o_res[3],
                                     int n_res[3],
                                     int o_min[3],
                                     int n_min[3],
                                     int o_max[3],
                                     int o_shift[3],
                                     int n_shift[3]);
void BKE_manta_cache_free(struct MantaDomainSettings *mds, struct Object *ob, int cache_map);

float BKE_manta_get_velocity_at(struct Object *ob, float position[3], float velocity[3]);
int BKE_manta_get_data_flags(struct MantaDomainSettings *mds);

#endif /* __BKE_MANTA_H__ */
