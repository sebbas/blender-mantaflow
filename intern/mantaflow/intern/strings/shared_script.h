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
// DEBUG
//////////////////////////////////////////////////////////////////////

const std::string manta_debuglevel = "\n\
def set_manta_debuglevel(level):\n\
    setDebugLevel(level=level)\n # level 0 = mute all output from manta\n";

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
mantaMsg('Fluid variables low')\n\
dim_s$ID$     = $SOLVER_DIM$\n\
res_s$ID$     = $RES$\n\
gravity_s$ID$ = vec3($GRAVITY_X$, $GRAVITY_Y$, $GRAVITY_Z$)\n\
gs_s$ID$      = vec3($RESX$, $RESY$, $RESZ$)\n\
\n\
if dim_s$ID$ == 2:\n\
    gs_s$ID$.z    = 1\n\
    gravity_s$ID$ = vec3($GRAVITY_X$,$GRAVITY_Z$,0)\n\
\n\
doOpen_s$ID$              = $DO_OPEN$\n\
boundConditions_s$ID$     = '$BOUNDCONDITIONS$'\n\
boundaryWidth_s$ID$       = 1\n\
\n\
using_highres_s$ID$   = $USING_HIGHRES$\n\
using_adaptTime_s$ID$ = $USING_ADAPTIVETIME$\n\
using_obstacle_s$ID$  = $USING_OBSTACLE$\n\
using_guiding_s$ID$   = $USING_GUIDING$\n\
using_invel_s$ID$     = $USING_INVEL$\n\
using_sndparts_s$ID$  = $USING_SNDPARTS$\n\
\n\
# fluid guiding params\n\
alpha_s$ID$ = $GUIDING_ALPHA$\n\
beta_s$ID$  = $GUIDING_BETA$\n\
tau_s$ID$   = 1.0\n\
sigma_s$ID$ = 0.99/tau_s$ID$\n\
theta_s$ID$ = 1.0\n\
\n\
# fluid diffusion / viscosity\n\
domainSize_s$ID$ = $FLUID_DOMAIN_SIZE$ # longest domain side in cm\n\
if domainSize_s$ID$ == 0: domainSize_s$ID$ = 100 # TODO (sebbas): just for versioning, remove with proper 2.8 versioning\n\
viscosity_s$ID$ = $FLUID_VISCOSITY$ / (domainSize_s$ID$*domainSize_s$ID$) # kinematic viscosity in m^2/s\n";

const std::string fluid_variables_high= "\n\
mantaMsg('Fluid variables high')\n\
upres_xl$ID$  = $UPRES$\n\
gs_xl$ID$     = vec3($HRESX$, $HRESY$, $HRESZ$)\n\
\n\
if dim_s$ID$ == 2:\n\
    gs_xl$ID$.z = 1\n";

const std::string fluid_with_obstacle = "\n\
using_obstacle_s$ID$ = True\n";

const std::string fluid_with_guiding = "\n\
using_guiding_s$ID$ = True\n";

const std::string fluid_with_invel = "\n\
using_invel_s$ID$ = True\n";

const std::string fluid_with_sndparts = "\n\
using_sndparts_s$ID$ = True\n";

//////////////////////////////////////////////////////////////////////
// ADAPTIVE TIME STEPPING
//////////////////////////////////////////////////////////////////////

const std::string fluid_adaptive_time_stepping_low = "\n\
mantaMsg('Fluid adaptive time stepping low')\n\
dt_default_s$ID$  = 0.1 # dt is 0.1 at 25fps\n\
dt_factor_s$ID$   = $DT_FACTOR$\n\
fps_s$ID$         = $FPS$\n\
dt0_s$ID$         = dt_default_s$ID$ * (25.0 / fps_s$ID$) * dt_factor_s$ID$\n\
s$ID$.frameLength = dt0_s$ID$ \n\
s$ID$.timestepMin = s$ID$.frameLength / 10.\n\
s$ID$.timestepMax = s$ID$.frameLength\n\
s$ID$.cfl         = $CFL$\n\
s$ID$.timestep    = s$ID$.frameLength\n";

const std::string fluid_adaptive_time_stepping_high = "\n\
mantaMsg('Adaptive time stepping high')\n\
xl$ID$.frameLength = s$ID$.frameLength\n\
xl$ID$.timestepMin = s$ID$.timestepMin\n\
xl$ID$.timestepMax = s$ID$.timestepMax\n\
xl$ID$.cfl         = s$ID$.cfl\n";

const std::string fluid_adapt_time_step_low = "\n\
def fluid_adapt_time_step_low():\n\
    maxVel_s$ID$ = vel_s$ID$.getMax() if vel_s$ID$ else 0\n\
    if using_adaptTime_s$ID$:\n\
        mantaMsg('Adapt timestep')\n\
        s$ID$.adaptTimestep(maxVel_s$ID$)\n";

const std::string fluid_adapt_time_step_high = "\n\
def fluid_adapt_time_step_high():\n\
    xl$ID$.timestep = s$ID$.timestep\n";

//////////////////////////////////////////////////////////////////////
// GRIDS
//////////////////////////////////////////////////////////////////////

const std::string fluid_alloc_low = "\n\
mantaMsg('Fluid alloc low')\n\
flags_s$ID$       = s$ID$.create(FlagGrid)\n\
vel_s$ID$         = s$ID$.create(MACGrid)\n\
x_vel_s$ID$       = s$ID$.create(RealGrid)\n\
y_vel_s$ID$       = s$ID$.create(RealGrid)\n\
z_vel_s$ID$       = s$ID$.create(RealGrid)\n\
pressure_s$ID$    = s$ID$.create(RealGrid)\n\
phiObs_s$ID$      = s$ID$.create(LevelsetGrid)\n\
phiOutIn_s$ID$    = s$ID$.create(LevelsetGrid)\n\
forces_s$ID$      = s$ID$.create(Vec3Grid)\n\
x_force_s$ID$     = s$ID$.create(RealGrid)\n\
y_force_s$ID$     = s$ID$.create(RealGrid)\n\
z_force_s$ID$     = s$ID$.create(RealGrid)\n";

const std::string fluid_alloc_high = "\n\
mantaMsg('Fluid alloc high')\n\
flags_xl$ID$     = xl$ID$.create(FlagGrid)\n";

const std::string fluid_alloc_obstacle_low = "\n\
mantaMsg('Allocating obstacle low')\n\
numObs_s$ID$     = s$ID$.create(IntGrid)\n\
phiObsIn_s$ID$   = s$ID$.create(LevelsetGrid)\n\
obvel_s$ID$      = s$ID$.create(MACGrid)\n\
obvelC_s$ID$     = s$ID$.create(Vec3Grid)\n\
x_obvel_s$ID$    = s$ID$.create(RealGrid)\n\
y_obvel_s$ID$    = s$ID$.create(RealGrid)\n\
z_obvel_s$ID$    = s$ID$.create(RealGrid)\n";

const std::string fluid_alloc_guiding_low = "\n\
mantaMsg('Allocating guiding low')\n\
numGuides_s$ID$   = s$ID$.create(IntGrid)\n\
phiGuideIn_s$ID$  = s$ID$.create(LevelsetGrid)\n\
guidevel_s$ID$    = s$ID$.create(MACGrid)\n\
guidevelC_s$ID$   = s$ID$.create(Vec3Grid)\n\
x_guidevel_s$ID$  = s$ID$.create(RealGrid)\n\
y_guidevel_s$ID$  = s$ID$.create(RealGrid)\n\
z_guidevel_s$ID$  = s$ID$.create(RealGrid)\n\
weightGuide_s$ID$ = s$ID$.create(RealGrid)\n";

const std::string fluid_alloc_invel_low = "\n\
mantaMsg('Allocating initial velocity low')\n\
invel_s$ID$   = s$ID$.create(VecGrid)\n\
x_invel_s$ID$ = s$ID$.create(RealGrid)\n\
y_invel_s$ID$ = s$ID$.create(RealGrid)\n\
z_invel_s$ID$ = s$ID$.create(RealGrid)\n";

const std::string fluid_alloc_sndparts_low = "\n\
mantaMsg('Allocating snd parts low')\n\
ppSnd_s$ID$     = s$ID$.create(BasicParticleSystem)\n\
pVelSnd_pp$ID$  = ppSnd_s$ID$.create(PdataVec3)\n\
pLifeSnd_pp$ID$ = ppSnd_s$ID$.create(PdataReal)\n\
pForceSnd_pp$ID$ = ppSnd_s$ID$.create(PdataVec3)\n";

//////////////////////////////////////////////////////////////////////
// DESTRUCTION
//////////////////////////////////////////////////////////////////////

const std::string fluid_delete_all = "\n\
mantaMsg('Deleting fluid')\n\
# Delete childs from pp object first\n\
for var in list(globals()):\n\
    if var.endswith('_pp$ID$'):\n\
        del globals()[var]\n\
# Now delete childs from s and xl object\n\
for var in list(globals()):\n\
    if var.endswith('_s$ID$') or var.endswith('_xl$ID$'):\n\
        del globals()[var]\n\
\n\
# Extra cleanup for multigrid and fluid guiding\n\
mantaMsg('Release multigrid')\n\
if 's$ID$' in globals(): releaseMG(s$ID$)\n\
if 'xl$ID$' in globals(): releaseMG(xl$ID$)\n\
mantaMsg('Release fluid guiding')\n\
releaseBlurPrecomp()\n\
\n\
# Release unreferenced memory (if there is some left, can in fact happen)\n\
gc.collect()\n\
\n\
# Now it is safe to delete solver objects (always need to be deleted last)\n\
mantaMsg('Delete solver low')\n\
if 's$ID$' in globals(): del s$ID$\n\
mantaMsg('Delete solver high')\n\
if 'xl$ID$' in globals(): del xl$ID$\n\
\n\
# Release unreferenced memory (if there is some left)\n\
gc.collect()\n";

const std::string fluid_delete_solver_low = "\n\
mantaMsg('Deleting solver low')\n\
if 's$ID$' in globals() : del s$ID$\n";

const std::string fluid_delete_solver_high = "\n\
mantaMsg('Deleting solver high')\n\
if 'xl$ID$' in globals() : del xl$ID$\n";

const std::string fluid_multigrid_cleanup_low = "\n\
mantaMsg('Cleanup multigrid low')\n\
if 's$ID$' in globals() : releaseMG(s$ID$)\n";

const std::string fluid_multigrid_cleanup_high = "\n\
mantaMsg('Cleanup multigrid high')\n\
if 'xl$ID$' in globals() : releaseMG(xl$ID$)\n";

const std::string fluid_guiding_cleanup_low = "\n\
mantaMsg('Cleanup guiding low')\n\
releaseBlurPrecomp()\n";

const std::string gc_collect = "\n\
gc.collect()\n";

//////////////////////////////////////////////////////////////////////
// IMPORT / EXPORT
//////////////////////////////////////////////////////////////////////

const std::string fluid_obstacle_import_low = "\n\
def load_fluid_obstacle_data_low_$ID$(path):\n\
    numObs_s$ID$.load(path + '_numObs.uni')\n\
    phiObsIn_s$ID$.load(path + '_phiObsIn.uni')\n\
    obvel_s$ID$.load(path + '_obvel.uni')\n\
    x_obvel_s$ID$.load(path + '_x_obvel.uni')\n\
    y_obvel_s$ID$.load(path + '_y_obvel.uni')\n\
    z_obvel_s$ID$.load(path + '_z_obvel.uni')\n";

const std::string fluid_guiding_import_low = "\n\
def load_fluid_guiding_data_low_$ID$(path):\n\
    numGuides_s$ID$.load(path + '_numGuides.uni')\n\
    phiGuideIn_s$ID$.load(path + '_phiGuideIn.uni')\n\
    guidevel_s$ID$.load(path + '_guidevel.uni')\n\
    x_guidevel_s$ID$.load(path + '_x_guidevel.uni')\n\
    y_guidevel_s$ID$.load(path + '_y_guidevel.uni')\n\
    z_guidevel_s$ID$.load(path + '_z_guidevel.uni')\n\
    weightGuide_s$ID$.load(path + '_weightGuide.uni')\n";

const std::string fluid_invel_import_low = "\n\
def load_fluid_invel_data_low_$ID$(path):\n\
    invel_s$ID$.load(path + '_invel.uni')\n\
    x_invel_s$ID$.load(path + '_x_invel.uni')\n\
    y_invel_s$ID$.load(path + '_y_invel.uni')\n\
    z_invel_s$ID$.load(path + '_z_invel.uni')\n";

const std::string fluid_sndparts_import_low = "\n\
def load_fluid_sndparts_data_low_$ID$(path):\n\
    ppSnd_s$ID$.load(path + '_ppSnd.uni')\n\
    pVelSnd_pp$ID$.load(path + '_pVelSnd.uni')\n\
    pLifeSnd_pp$ID$.load(path + '_pLifeSnd.uni')\n";

const std::string fluid_obstacle_export_low = "\n\
def save_fluid_obstacle_data_low_$ID$(path):\n\
    numObs_s$ID$.save(path + '_numObs.uni')\n\
    phiObsIn_s$ID$.save(path + '_phiObsIn.uni')\n\
    obvel_s$ID$.save(path + '_obvel.uni')\n\
    x_obvel_s$ID$.save(path + '_x_obvel.uni')\n\
    y_obvel_s$ID$.save(path + '_y_obvel.uni')\n\
    z_obvel_s$ID$.save(path + '_z_obvel.uni')\n";

const std::string fluid_guiding_export_low = "\n\
def save_fluid_guiding_data_low_$ID$(path):\n\
    numGuides_s$ID$.save(path + '_numGuides.uni')\n\
    phiGuideIn_s$ID$.save(path + '_phiGuideIn.uni')\n\
    guidevel_s$ID$.save(path + '_guidevel.uni')\n\
    x_guidevel_s$ID$.save(path + '_x_guidevel.uni')\n\
    y_guidevel_s$ID$.save(path + '_y_guidevel.uni')\n\
    z_guidevel_s$ID$.save(path + '_z_guidevel.uni')\n\
    weightGuide_s$ID$.save(path + '_weightGuide.uni')\n";

const std::string fluid_invel_export_low = "\n\
def save_fluid_invel_data_low_$ID$(path):\n\
    invel_s$ID$.save(path + '_invel.uni')\n\
    x_invel_s$ID$.save(path + '_x_invel.uni')\n\
    y_invel_s$ID$.save(path + '_y_invel.uni')\n\
    z_invel_s$ID$.save(path + '_z_invel.uni')\n";

const std::string fluid_sndparts_export_low = "\n\
def save_fluid_sndparts_data_low_$ID$(path):\n\
    ppSnd_s$ID$.save(path + '_ppSnd.uni')\n\
    pVelSnd_pp$ID$.save(path + '_pVelSnd.uni')\n\
    pLifeSnd_pp$ID$.save(path + '_pLifeSnd.uni')\n";

//////////////////////////////////////////////////////////////////////
// STANDALONE MODE
//////////////////////////////////////////////////////////////////////

const std::string fluid_standalone_load = "\n\
if using_obstacle_s$ID$:\n\
    load_fluid_obstacle_data_low_$ID$(path_prefix_$ID$)\n\
if using_guiding_s$ID$:\n\
    load_fluid_guiding_data_low_$ID$(path_prefix_$ID$)\n\
if using_invel_s$ID$:\n\
    load_fluid_invel_data_low_$ID$(path_prefix_$ID$)\n\
if using_sndparts_s$ID$:\n\
    load_fluid_sndparts_data_low_$ID$(path_prefix_$ID$)\n";

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
