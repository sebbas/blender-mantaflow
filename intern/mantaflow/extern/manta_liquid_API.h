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

/** \file mantaflow/extern/manta_liquid_API.h
 *  \ingroup mantaflow
 */

#ifndef MANTA_LIQUID_API_H_
#define MANTA_LIQUID_API_H_

#ifdef __cplusplus
extern "C" {
#endif

struct LIQUID;
struct LIQUID *liquid_init();
void liquid_step(struct LIQUID *liquid, int currentFrame);

#ifdef __cplusplus
}
#endif

#endif /* MANTA_LIQUID_API_H_ */