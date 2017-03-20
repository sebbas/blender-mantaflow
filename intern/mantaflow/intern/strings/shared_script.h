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
s$ID$ = Solver(name='solver_s$ID$', gridSize=gs_s$ID$, dim=dim_s$ID$)\n";

const std::string fluid_solver_high = "\n\
mantaMsg('Solver high')\n\
xl$ID$ = Solver(name='solver_xl$ID$', gridSize=gs_xl$ID$)\n";

//////////////////////////////////////////////////////////////////////
// VARIABLES
//////////////////////////////////////////////////////////////////////

const std::string fluid_variables_low = "\n\
dim_s$ID$     = $SOLVER_DIM$\n\
res_s$ID$     = $RES$\n\
gravity_s$ID$ = vec3($GRAVITY_X$, $GRAVITY_Y$, $GRAVITY_Z$)\n\
gs_s$ID$      = vec3($RESX$, $RESY$, $RESZ$)\n\
\n\
if dim_s$ID$ == 2:\n\
    gs_s$ID$.z    = 1\n\
    gravity_s$ID$ = vec3($GRAVITY_X$,$GRAVITY_Z$,0)\n\
\n\
doOpen_s$ID$          = $DO_OPEN$\n\
boundConditions_s$ID$ = '$BOUNDCONDITIONS$'\n\
boundaryWidth_s$ID$   = 1\n\
\n\
using_highres_s$ID$   = $USING_HIGHRES$\n\
using_adaptTime_s$ID$ = True # adaptive time stepping disabled for now\n";

const std::string fluid_variables_high= "\n\
upres_xl$ID$  = $UPRES$\n\
gs_xl$ID$     = vec3($HRESX$, $HRESY$, $HRESZ$)\n\
\n\
if dim_s$ID$ == 2:\n\
    gs_xl$ID$.z = 1\n";

//////////////////////////////////////////////////////////////////////
// ADAPTIVE TIME STEPPING
//////////////////////////////////////////////////////////////////////

const std::string fluid_adaptive_time_stepping_low = "\n\
mantaMsg('Adaptive time stepping low')\n\
dt_default_s$ID$  = 0.1\n\
dt_factor_s$ID$   = $DT_FACTOR$\n\
fps_s$ID$         = $FPS$\n\
dt0_s$ID$         = dt_default_s$ID$ * (25.0 / fps_s$ID$) * dt_factor_s$ID$\n\
s$ID$.frameLength = dt0_s$ID$\n\
s$ID$.timestepMin = dt0_s$ID$ / 10\n\
s$ID$.timestepMax = dt0_s$ID$\n\
s$ID$.cfl         = 4.0\n\
s$ID$.timestep    = dt0_s$ID$\n";

const std::string fluid_adaptive_time_stepping_high = "\n\
mantaMsg('Adaptive time stepping high')\n\
xl$ID$.frameLength = s$ID$.frameLength\n\
xl$ID$.timestepMin = s$ID$.timestepMin\n\
xl$ID$.timestepMax = s$ID$.timestepMax\n\
xl$ID$.cfl         = s$ID$.cfl\n";

//////////////////////////////////////////////////////////////////////
// DESTRUCTION
//////////////////////////////////////////////////////////////////////

const std::string fluid_delete_variables_low = "\n\
mantaMsg('Deleting fluid variables low')\n\
if 'dim_s$ID$'             in globals() : del dim_s$ID$\n\
if 'res_s$ID$'             in globals() : del res_s$ID$\n\
if 'gs_s$ID$'              in globals() : del gs_s$ID$\n\
if 'gravity_s$ID$'         in globals() : del gravity_s$ID$\n\
if 'doOpen_s$ID$'          in globals() : del doOpen_s$ID$\n\
if 'boundConditions_s$ID$' in globals() : del boundConditions_s$ID$\n\
if 'boundaryWidth_s$ID$'   in globals() : del boundaryWidth_s$ID$\n\
if 'dt_default_s$ID$'      in globals() : del dt_default_s$ID$\n\
if 'dt_factor_s$ID$'       in globals() : del dt_factor_s$ID$\n\
if 'fps_s$ID$'             in globals() : del fps_s$ID$\n\
if 'dt0_s$ID$'             in globals() : del dt0_s$ID$\n";

const std::string fluid_delete_variables_high = "\n\
mantaMsg('Deleting fluid variables high')\n\
if 'upres_xl$ID$'          in globals() : del upres_xl$ID$\n\
if 'gs_xl$ID$'             in globals() : del gs_xl$ID$\n";

const std::string fluid_delete_solver_low = "\n\
mantaMsg('Deleting solver low')\n\
if 's$ID$' in globals() : del s$ID$\n";

const std::string fluid_delete_solver_high = "\n\
mantaMsg('Deleting solver high')\n\
if 'xl$ID$' in globals() : del xl$ID$\n";

const std::string gc_collect = "\n\
gc.collect()\n";

//////////////////////////////////////////////////////////////////////
// STANDALONE MODE
//////////////////////////////////////////////////////////////////////

const std::string fluid_standalone = "\n\
if (GUI):\n\
    gui=Gui()\n\
    gui.show()\n\
    gui.pause()\n\
\n\
start_frame = $CURRENT_FRAME$\n\
end_frame = 1000\n\
\n\
# All low and high res steps\n\
while start_frame <= end_frame:\n\
    manta_step_$ID$(start_frame)\n\
    start_frame += 1\n";
