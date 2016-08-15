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

//////////////////////////////////////////////////////////////////////
// LIBRARIES
//////////////////////////////////////////////////////////////////////

const std::string manta_import = "\
from manta import *\n\
import os, shutil, math, sys, gc\n";

//////////////////////////////////////////////////////////////////////
// SOLVERS
//////////////////////////////////////////////////////////////////////

const std::string fluid_solver_low = "\n\
mantaMsg('Solver low')\n\
dim     = $SOLVER_DIM$\n\
res     = $RES$\n\
gravity = vec3($GRAVITY_X$, $GRAVITY_Y$, $GRAVITY_Z$)\n\
gs      = vec3($RESX$, $RESY$, $RESZ$)\n\
\n\
if dim == 2:\n\
    gs.z    = 1\n\
    gravity = vec3($GRAVITY_X$,$GRAVITY_Z$,0)\n\
\n\
doOpen          = $DO_OPEN$\n\
boundConditions = '$BOUNDCONDITIONS$'\n\
boundaryWidth   = 0\n\
\n\
s = Solver(name='main', gridSize=gs, dim=dim)\n";

const std::string fluid_solver_high = "\n\
mantaMsg('Solver high')\n\
upres  = $UPRES$\n\
xl_gs  = vec3($HRESX$, $HRESY$, $HRESZ$)\n\
\n\
if dim == 2:\n\
    xl_gs.z = 1\n\
\n\
xl = Solver(name='larger', gridSize=xl_gs)\n";

//////////////////////////////////////////////////////////////////////
// ADAPTIVE TIME STEPPING
//////////////////////////////////////////////////////////////////////

const std::string fluid_adaptive_time_stepping_low = "\n\
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

const std::string fluid_adaptive_time_stepping_high = "\n\
mantaMsg('Adaptive time stepping high')\n\
xl.frameLength = s.frameLength\n\
xl.timestepMin = s.timestepMin / 10\n\
xl.timestepMax = s.timestepMax\n\
xl.cfl         = s.cfl\n";

//////////////////////////////////////////////////////////////////////
// DESTRUCTION
//////////////////////////////////////////////////////////////////////

const std::string fluid_delete_variables_low = "\n\
if 'dim'             in globals() : del dim\n\
if 'res'             in globals() : del gravity\n\
if 'gravity'         in globals() : del gravity\n\
if 'gs'              in globals() : del gs\n\
if 'gravity'         in globals() : del gravity\n\
if 'gs'              in globals() : del gs\n\
if 'doOpen'          in globals() : del doOpen\n\
if 'boundConditions' in globals() : del boundConditions\n\
if 'boundaryWidth'   in globals() : del boundaryWidth\n\
if 'dt_default'      in globals() : del dt_default\n\
if 'dt_factor'       in globals() : del dt_factor\n\
if 'fps'             in globals() : del fps\n\
if 'dt0'             in globals() : del dt0\n";

const std::string fluid_delete_variables_high = "\n\
if 'upres'           in globals() : del upres\n\
if 'xl_gs'           in globals() : del xl_gs\n";

const std::string fluid_delete_solver_low = "\n\
mantaMsg('Deleting solver low')\n\
if 's' in globals() : del s\n";

const std::string fluid_delete_solver_high = "\n\
mantaMsg('Deleting solver high')\n\
if 'xl' in globals() : del xl\n";

const std::string gc_collect = "\n\
gc.collect()\n";
