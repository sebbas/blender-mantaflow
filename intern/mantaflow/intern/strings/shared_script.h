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

/** \file mantaflow/intern/strings/shared_script.h
 *  \ingroup mantaflow
 */

#include <string>

const std::string manta_import = "\
from manta import *\n\
import os, shutil, math, sys, gc\n";

const std::string solver_low = "\n\
# solver low params\n\
mantaMsg('Solver low')\n\
dim    = $SOLVER_DIM$\n\
res    = $RES$\n\
gs     = vec3($RESX$,$RESY$,$RESZ$)\n\
if dim == 2: gs.z = 1\n\
s      = Solver(name='main', gridSize=gs, dim=dim)\n";

const std::string solver_high = "\n\
# solver high params\n\
mantaMsg('Solver high')\n\
upres  = $UPRES$\n\
xl_gs  = vec3($HRESX$, $HRESY$, $HRESZ$)\n\
if dim == 2: xl_gs.z = 1\n\
xl     = Solver(name='larger', gridSize=xl_gs)\n";

const std::string adaptive_time_stepping_low = "\n\
# adaptive time stepping\n\
mantaMsg('Adaptive time stepping low')\n\
dt_default    = 0.1\n\
dt_factor     = $DT_FACTOR$\n\
fps           = $FPS$\n\
dt0           = dt_default * (25.0 / fps) * dt_factor\n\
s.frameLength = dt0\n\
s.timestepMin = dt0 / 10\n\
s.timestepMax = dt0\n\
s.cfl         = 4.0\n\
s.timestep    = dt0\n";

const std::string adaptive_time_stepping_high = "\n\
mantaMsg('Adaptive time stepping high')\n\
xl.frameLength = s.frameLength\n\
xl.timestepMin = s.timestepMin / 10\n\
xl.timestepMax = s.timestepMax\n\
xl.cfl         = s.cfl\n";

const std::string fluid_variables = "\n\
doOpen          = $DO_OPEN$\n\
boundConditions = '$BOUNDCONDITIONS$'\n\
boundaryWidth   = 1\n";

const std::string prep_domain_low = "\n\
# prepare domain low\n\
mantaMsg('Fluid domain low')\n\
flags.initDomain(boundaryWidth=boundaryWidth)\n\
flags.fillGrid()\n\
if doOpen:\n\
    setOpenBound(flags=flags, bWidth=boundaryWidth, openBound=boundConditions, type=FlagOutflow|FlagEmpty)\n";

const std::string prep_domain_high = "\n\
# prepare domain high\n\
mantaMsg('Fluid domain high')\n\
xl_flags.initDomain(boundaryWidth=boundaryWidth)\n\
xl_flags.fillGrid()\n\
if doOpen:\n\
    setOpenBound(flags=xl_flags, bWidth=boundaryWidth, openBound=boundConditions, type=FlagOutflow|FlagEmpty)\n";

