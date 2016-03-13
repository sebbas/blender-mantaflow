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

/** \file mantaflow/intern/strings/smoke.h
 *  \ingroup mantaflow
 */

#include <string>
using namespace std;

//////////////////////////////////////////////////////////////////////
// GENERAL SETUP
//////////////////////////////////////////////////////////////////////

const string manta_import = "\
from manta import *\n\
import os, shutil, math, sys\n";

const string flags = "\n\
using_colors = $USING_COLORS$\n\
using_heat = $USING_HEAT$\n\
using_fire = $USING_FIRE$\n\
low_flags_updated = False\n\
using_wavelets = $USE_WAVELETS$\n";

const string uv_setup = "\n\
# create the array of uv grids\n\
uv = []\n\
for i in range(uvs):\n\
  uvGrid = s.create(VecGrid)\n\
  uv.append(uvGrid)\n\
  resetUvGrid(uv[i])\n";

//////////////////////////////////////////////////////////////////////
// LOW RESOLUTION SETUP
//////////////////////////////////////////////////////////////////////

const string solver_setup_low = "\n\
# solver low params\n\
dim = $SOLVER_DIM$\n\
doOpen = $DO_OPEN$\n\
boundConditions = '$BOUNDCONDITIONS$'\n\
res = $RES$\n\
gs = vec3($RESX$,$RESY$,$RESZ$)\n\
if dim == 2:\n\
  gs.z = 1\n\
s = FluidSolver(name='main', gridSize=gs, dim=dim)\n\
dt_default = 0.1\n\
dt_factor = $DT_FACTOR$\n\
fps = $FPS$\n\
dt0 = dt_default * (25.0 / fps) * dt_factor\n\
s.frameLength = dt0\n\
s.timestepMin = dt0 / 10\n\
s.timestepMax = dt0\n\
s.cfl = 4.0\n\
s.timestep = dt0\n\
timings = Timings()\n\
vorticity = $VORTICITY$\n\
boundaryWidth = 1\n\
uvs = 2\n";

const string alloc_base_grids_low = "\n\
# prepare grids low\n\
flags = s.create(FlagGrid)\n\
vel = s.create(MACGrid)\n\
x_vel = s.create(RealGrid)\n\
y_vel = s.create(RealGrid)\n\
z_vel = s.create(RealGrid)\n\
density = s.create(LevelsetGrid)\n\
pressure = s.create(RealGrid)\n\
energy = s.create(RealGrid)\n\
forces = s.create(MACGrid)\n\
inflow_grid = s.create(LevelsetGrid)\n\
fuel_inflow = s.create(LevelsetGrid)\n";

const string prep_domain_low = "\n\
# prepare domain low\n\
flags.initDomain(boundaryWidth=boundaryWidth)\n\
flags.fillGrid()\n\
if doOpen:\n\
  setOpenBound(flags=flags, bWidth=boundaryWidth, openBound=boundConditions, type=FlagOutflow|FlagEmpty)\n";

//////////////////////////////////////////////////////////////////////
// HIGH RESOLUTION SETUP
//////////////////////////////////////////////////////////////////////

const string solver_setup_high = "\n\
# solver high params\n\
upres = $UPRES$\n\
xl_gs = vec3($HRESX$, $HRESY$, $HRESZ$)\n\
if dim == 2:\n\
  xl_gs.z = 1\n\
xl = Solver(name = 'larger', gridSize = xl_gs)\n\
xl.frameLength = s.frameLength\n\
xl.timestepMin = s.timestepMin / 10\n\
xl.timestepMax = s.timestepMax\n\
xl.cfl = s.cfl\n\
wltStrength = $WLT_STR$\n\
octaves = 0\n\
if upres == 1:\n\
  octaves = int(math.log(upres+1)/ math.log(2.0) + 0.5)\n\
elif upres > 1:\n\
  octaves = int(math.log(upres)/ math.log(2.0) + 0.5)\n";

const string alloc_base_grids_high = "\n\
# prepare grids high\n\
xl_flags = xl.create(FlagGrid)\n\
xl_vel = xl.create(MACGrid)\n\
xl_x_vel = s.create(RealGrid)\n\
xl_y_vel = s.create(RealGrid)\n\
xl_z_vel = s.create(RealGrid)\n\
xl_density = xl.create(RealGrid)\n\
xl_weight  = xl.create(RealGrid)\n";

const string prep_domain_high = "\n\
# prepare domain high\n\
xl_flags.initDomain(boundaryWidth=boundaryWidth)\n\
xl_flags.fillGrid()\n\
if doOpen:\n\
  setOpenBound(flags=xl_flags, bWidth=boundaryWidth, openBound=boundConditions, type=FlagOutflow|FlagEmpty)\n";

const string wavelet_turbulence_noise = "\n\
# wavelet turbulence noise field\n\
xl_wltnoise = s.create(NoiseField, loadFromFile=True)\n\
xl_wltnoise.posScale = vec3(int(1.0*gs.x)) / $NOISE_POSSCALE$\n\
xl_wltnoise.timeAnim = $NOISE_TIMEANIM$\n\
if(upres>0):\n\
  xl_wltnoise.posScale = xl_wltnoise.posScale * (1./upres)\n\
  xl_wltnoise.timeAnim = xl_wltnoise.timeAnim * upres\n";

//////////////////////////////////////////////////////////////////////
// ADDITIONAL GRIDS
//////////////////////////////////////////////////////////////////////

const string alloc_colors_low = "\n\
mantaMsg('Allocating colors low')\n\
color_r = s.create(RealGrid)\n\
color_g = s.create(RealGrid)\n\
color_b = s.create(RealGrid)\n";

const string alloc_colors_high = "\
mantaMsg('Allocating colors high')\n\
xl_color_r = xl.create(RealGrid)\n\
xl_color_g = xl.create(RealGrid)\n\
xl_color_b = xl.create(RealGrid)\n";

const string init_colors_low = "\n\
mantaMsg('Initializing colors low')\n\
color_r.copyFrom(density) \n\
color_r.multConst(manta_color_r) \n\
color_g.copyFrom(density) \n\
color_g.multConst(manta_color_g) \n\
color_b.copyFrom(density) \n\
color_b.multConst(manta_color_b)\n";

const string init_colors_high = "\n\
mantaMsg('Initializing colors high')\n\
xl_color_r.copyFrom(xl_density) \n\
xl_color_r.multConst(manta_color_r) \n\
xl_color_g.copyFrom(xl_density) \n\
xl_color_g.multConst(manta_color_g) \n\
xl_color_b.copyFrom(xl_density) \n\
xl_color_b.multConst(manta_color_b)\n";

const string alloc_heat_low = "\n\
mantaMsg('Allocating heat low')\n\
heat = s.create(RealGrid)\n";

const string alloc_fire_low = "\n\
mantaMsg('Allocating fire low')\n\
flame = s.create(RealGrid)\n\
fuel = s.create(RealGrid)\n\
react = s.create(RealGrid)\n";

const string alloc_fire_high = "\n\
mantaMsg('Allocating fire high')\n\
xl_flame = xl.create(RealGrid)\n\
xl_fuel = xl.create(RealGrid)\n\
xl_react = xl.create(RealGrid)\n";

const string with_heat = "\n\
using_heat = True\n";

const string with_colors = "\n\
using_colors = True\n";

const string with_fire = "\n\
using_fire = True\n";

//////////////////////////////////////////////////////////////////////
// STANDALONE MODE
//////////////////////////////////////////////////////////////////////

const string standalone = "\n\
if (GUI):\n\
  gui=Gui()\n\
  gui.show()\n\
  gui.pause()\n\
\n\
# import *.uni files\n\
import_grids_low()\n\
if using_wavelets:\n\
  import_grids_high()\n\
\n\
for step in range(1000):\n\
  apply_inflow()\n\
  \n\
  mantaMsg('Low step '+ str(step))\n\
  if using_fire:\n\
    process_burn_low()\n\
  step_low()\n\
  if using_fire:\n\
    update_flame_low()\n\
  \n\
  if using_wavelets:\n\
    mantaMsg('High step '+ str(step))\n\
    if using_fire:\n\
      process_burn_high()\n\
    step_high()\n\
    if using_fire:\n\
      update_flame_high()\n";

//////////////////////////////////////////////////////////////////////
// DESTRUCTION
//////////////////////////////////////////////////////////////////////

const string del_colors_low = "\n\
mantaMsg('Deleting colors low')\n\
if 'color_r' in globals() : del color_r\n\
if 'color_g' in globals() : del color_g\n\
if 'color_b' in globals() : del color_b\n";

const string del_colors_high = "\n\
mantaMsg('Deleting colors high')\n\
if 'xl_color_r' in globals() : del xl_color_r\n\
if 'xl_color_g' in globals() : del xl_color_g\n\
if 'xl_color_b' in globals() : del xl_color_b\n";

const string del_fire_low = "\n\
mantaMsg('Deleting fire low')\n\
if 'flame' in globals() : del flame\n\
if 'fuel' in globals() : del fuel\n\
if 'react' in globals() : del react\n";

const string del_fire_high = "\n\
mantaMsg('Deleting fire high')\n\
if 'xl_flame' in globals() : del xl_flame\n\
if 'xl_fuel' in globals() : del xl_fuel\n\
if 'xl_react' in globals() : del xl_react\n";

const string del_heat_low = "\n\
mantaMsg('Deleting heat low')\n\
if 'heat' in globals() : del heat\n";

const string del_base_grids_low = "\n\
mantaMsg('Deleting base grids low')\n\
if 'flags' in globals() : del flags\n\
if 'uvs' in globals() : del uvs\n\
if 'vel' in globals() : del vel\n\
if 'x_vel' in globals() : del x_vel\n\
if 'y_vel' in globals() : del y_vel\n\
if 'z_vel' in globals() : del z_vel\n\
if 'density' in globals() : del density\n\
if 'pressure' in globals() : del pressure\n\
if 'energy' in globals() : del energy\n\
if 'tempFlag' in globals() : del tempFlag\n\
if 'forces' in globals() : del forces\n\
if 'inflow_grid' in globals() : del inflow_grid\n\
if 'fuel_inflow' in globals() : del fuel_inflow\n";

const string del_base_grids_high = "\n\
mantaMsg('Deleting base grids high')\n\
if 'xl_flags' in globals() : del xl_flags\n\
if 'xl_vel' in globals() : del xl_vel\n\
if 'xl_x_vel' in globals() : del xl_x_vel\n\
if 'xl_y_vel' in globals() : del xl_y_vel\n\
if 'xl_z_vel' in globals() : del xl_z_vel\n\
if 'xl_density' in globals() : del xl_density\n\
if 'xl_wltnoise' in globals() : del xl_wltnoise\n";

const string del_vars_low = "\n\
mantaMsg('Deleting variables low')\n\
if 'res' in globals() : del res\n\
if 'dim' in globals() : del dim\n\
if 'gs' in globals() : del gs\n\
if 'doOpen' in globals() : del doOpen\n\
if 'boundConditions' in globals() : del boundConditions\n\
if 'dt_default' in globals() : del dt_default\n\
if 'dt_factor' in globals() : del dt_factor\n\
if 'fps' in globals() : del fps\n\
if 'dt0' in globals() : del dt0\n\
if 'vorticity' in globals() : del vorticity\n\
if 'boundaryWidth' in globals() : del boundaryWidth\n\
if 's' in globals() : del s\n\
if 'timings' in globals() : del timings\n\
if 'using_colors' in globals() : del using_colors\n\
if 'using_heat' in globals() : del using_heat\n\
if 'using_fire' in globals() : del using_fire\n\
if 'last_frame' in globals() : del last_frame\n\
if 'maxvel' in globals() : del maxvel\n\
if 'gravity' in globals() : del gravity\n\
if 'uv' in globals() : del uv\n";

const string del_vars_high = "\n\
mantaMsg('Deleting variables high')\n\
if 'upres' in globals() : del upres\n\
if 'xl_gs' in globals() : del xl_gs\n\
if 'xl' in globals() : del xl\n\
if 'wltStrength' in globals() : del wltStrength\n\
if 'octaves' in globals() : del octaves\n";

//////////////////////////////////////////////////////////////////////
// MANTA STEP
//////////////////////////////////////////////////////////////////////

const string manta_step = "\n\
def manta_step():\n\
  last_frame = s.frame\n\
  while s.frame == last_frame:\n\
    mantaMsg('Adapt timestep')\n\
    maxvel = vel.getMaxValue()\n\
    s.adaptTimestep(maxvel)\n\
    \n\
    mantaMsg('Low step / s.frame: ' + str(s.frame))\n\
    if using_fire:\n\
      process_burn_low()\n\
    step_low()\n\
    if using_fire:\n\
      update_flame_low()\n\
    \n\
    if using_wavelets:\n\
      xl.timestep = s.timestep\n\
      mantaMsg('High step / s.frame: ' + str(s.frame))\n\
      if using_fire:\n\
        process_burn_high()\n\
      step_high()\n\
      if using_fire:\n\
        update_flame_high()\n\
    s.step()\n";

//////////////////////////////////////////////////////////////////////
// STEP FUNCTIONS LOW
//////////////////////////////////////////////////////////////////////

const string smoke_step_low = "\n\
def step_low():\n\
  mantaMsg('Step low')\n\
  copyRealToMac(sourceX=x_vel, sourceY=y_vel, sourceZ=z_vel, target=vel)\n\
  if dim == 2:\n\
    density.add(inflow_grid)\n\
  \n\
  mantaMsg('Advecting density')\n\
  advectSemiLagrange(flags=flags, vel=vel, grid=density, order=$ADVECT_ORDER$)\n\
  \n\
  if using_heat:\n\
    mantaMsg('Advecting heat')\n\
    advectSemiLagrange(flags=flags, vel=vel, grid=heat, order=$ADVECT_ORDER$)\n\
  \n\
  if using_fire:\n\
    mantaMsg('Advecting fire')\n\
    advectSemiLagrange(flags=flags, vel=vel, grid=fuel, order=$ADVECT_ORDER$)\n\
    advectSemiLagrange(flags=flags, vel=vel, grid=react, order=$ADVECT_ORDER$)\n\
  \n\
  if using_colors:\n\
    mantaMsg('Advecting colors')\n\
    advectSemiLagrange(flags=flags, vel=vel, grid=color_r, order=$ADVECT_ORDER$)\n\
    advectSemiLagrange(flags=flags, vel=vel, grid=color_g, order=$ADVECT_ORDER$)\n\
    advectSemiLagrange(flags=flags, vel=vel, grid=color_b, order=$ADVECT_ORDER$)\n\
  \n\
  mantaMsg('Advecting velocity')\n\
  advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=$ADVECT_ORDER$, openBounds=doOpen, boundaryWidth=boundaryWidth)\n\
  \n\
  for i in range(uvs):\n\
    mantaMsg('Advecting UV and updating UVWeight')\n\
    advectSemiLagrange(flags=flags, vel=vel, grid=uv[i], order=$ADVECT_ORDER$)\n\
    updateUvWeight(resetTime=16.5 , index=i, numUvs=uvs, uv=uv[i])\n\
  \n\
  if doOpen:\n\
    resetOutflow(flags=flags, real=density)\n\
  mantaMsg('Vorticity')\n\
  if vorticity > 0.01:\n\
    vorticityConfinement(vel=vel, flags=flags, strength=$VORTICITY$)\n\
  \n\
  if using_heat:\n\
    mantaMsg('Adding heat buoyancy')\n\
    gravity=vec3(0,0,-1) if dim==3 else vec3(0,-0.0981,0)\n\
    addBuoyancy(flags=flags, density=density, vel=vel, gravity=gravity, coefficient=$ALPHA$*(-1))\n\
    addBuoyancy(flags=flags, density=heat, vel=vel, gravity=gravity, coefficient=$BETA$*(-1))\n\
  else:\n\
    mantaMsg('Adding buoyancy')\n\
    gravity=vec3(0,0,-0.01 * $ALPHA$) if dim==3 else vec3(0,-0.01* $ALPHA$,0)\n\
    addBuoyancy(density=density, vel=vel, gravity=gravity, flags=flags)\n\
  \n\
  mantaMsg('Walls')\n\
  setWallBcs(flags=flags, vel=vel)\n\
  \n\
  mantaMsg('Pressure')\n\
  solvePressure(flags=flags, vel=vel, pressure=pressure)\n\
  \n\
  mantaMsg('Energy')\n\
  computeEnergy(flags=flags, vel=vel, energy=energy)\n\
  \n\
  # TODO: mantaMsg('Forcefield')\n\
  # TODO: addForceField(flags=flags, vel=vel, force=forces)\n\
  # TODO: forces.clear()\n\
  copyMacToReal(source=vel, targetX=x_vel, targetY=y_vel, targetZ=z_vel)\n\
  # --> removing solver step for now from here \n\
  # s.step()\n\
\n\
def process_burn_low():\n\
  mantaMsg('Process burn low')\n\
  if (using_colors):\n\
    processBurn(fuel=fuel, density=density, react=react, red=color_r, green=color_g, blue=color_b, heat=heat, burningRate=$BURNING_RATE$, flameSmoke=$FLAME_SMOKE$, ignitionTemp=$IGNITION_TEMP$, maxTemp=$MAX_TEMP$, dt=dt0, flameSmokeColor=vec3($FLAME_SMOKE_COLOR_X$,$FLAME_SMOKE_COLOR_Y$,$FLAME_SMOKE_COLOR_Z$))\n\
  else:\n\
    processBurn(fuel=fuel, density=density, react=react, heat=heat, burningRate=$BURNING_RATE$, flameSmoke=$FLAME_SMOKE$, ignitionTemp=$IGNITION_TEMP$, maxTemp=$MAX_TEMP$, dt=dt0, flameSmokeColor=vec3($FLAME_SMOKE_COLOR_X$,$FLAME_SMOKE_COLOR_Y$,$FLAME_SMOKE_COLOR_Z$))\n\
\n\
def update_flame_low():\n\
  mantaMsg('Update flame low')\n\
  updateFlame(react=react, flame=flame)\n";

//////////////////////////////////////////////////////////////////////
// STEP FUNCTIONS HIGH
//////////////////////////////////////////////////////////////////////

const string smoke_step_high = "\n\
def step_high():\n\
  mantaMsg('Step high')\n\
  interpolateMACGrid(source=vel, target=xl_vel)\n\
  sStr = 1.0 * wltStrength\n\
  sPos = 2.0\n\
  \n\
  mantaMsg('Octaves')\n\
  for o in range(octaves):\n\
    for i in range(uvs):\n\
      uvWeight = getUvWeight(uv[i])\n\
      applyNoiseVec3(flags=xl_flags, target=xl_vel, noise=xl_wltnoise, scale=sStr * uvWeight, scaleSpatial=sPos , weight=energy, uv=uv[i])\n\
    sStr *= 0.06 # magic kolmogorov factor \n\
    sPos *= 2.0 \n\
  \n\
  for substep in range(upres):\n\
    if using_colors: \n\
      # mantaMsg ('Advecting colors high')\n\
      advectSemiLagrange(flags=xl_flags, vel=xl_vel, grid=xl_color_r, order=$ADVECT_ORDER$, openBounds=doOpen)\n\
      advectSemiLagrange(flags=xl_flags, vel=xl_vel, grid=xl_color_g, order=$ADVECT_ORDER$, openBounds=doOpen)\n\
      advectSemiLagrange(flags=xl_flags, vel=xl_vel, grid=xl_color_b, order=$ADVECT_ORDER$, openBounds=doOpen)\n\
    \n\
    if using_fire: \n\
      # mantaMsg ('Advecting fire high')\n\
      advectSemiLagrange(flags=xl_flags, vel=xl_vel, grid=xl_fuel, order=$ADVECT_ORDER$, openBounds=doOpen)\n\
      advectSemiLagrange(flags=xl_flags, vel=xl_vel, grid=xl_react, order=$ADVECT_ORDER$, openBounds=doOpen)\n\
    \n\
    mantaMsg('Advecting density high')\n\
    advectSemiLagrange(flags=xl_flags, vel=xl_vel, grid=xl_density, order=$ADVECT_ORDER$, openBounds=doOpen)\n\
  \n\
  # --> removing solver step for now from here\n\
  # xl.step()\n\
\n\
def process_burn_high():\n\
  mantaMsg('Process burn high')\n\
  if (using_colors):\n\
    processBurn(fuel=xl_fuel, density=xl_density, react=xl_react, red=xl_color_r, green=xl_color_g, blue=xl_color_b, burningRate=$BURNING_RATE$, flameSmoke=$FLAME_SMOKE$, ignitionTemp=$IGNITION_TEMP$, maxTemp=$MAX_TEMP$, dt=dt0, flameSmokeColor=vec3($FLAME_SMOKE_COLOR_X$,$FLAME_SMOKE_COLOR_Y$,$FLAME_SMOKE_COLOR_Z$))\n\
  else:\n\
    processBurn(fuel=xl_fuel, density=xl_density, react=xl_react, burningRate=$BURNING_RATE$, flameSmoke=$FLAME_SMOKE$, ignitionTemp=$IGNITION_TEMP$, maxTemp=$MAX_TEMP$, dt=dt0, flameSmokeColor=vec3($FLAME_SMOKE_COLOR_X$,$FLAME_SMOKE_COLOR_Y$,$FLAME_SMOKE_COLOR_Z$))\n\
\n\
def update_flame_high():\n\
  mantaMsg('Update flame high')\n\
  updateFlame(react=xl_react, flame=xl_flame)\n";

//////////////////////////////////////////////////////////////////////
// STEP FUNCTIONS LIQUID
//////////////////////////////////////////////////////////////////////

const string liquid_step_low = "\n\
def sim_step_low(t):\n\
#update flags from density on first step\n\
  setWallBcs(flags=flags, vel=vel)\n\
  density.multConst(-1.)\n\
  mantaMsg(using_colors)\n\
  global low_flags_updated\n\
  if not low_flags_updated:\n\
    mantaMsg('Updating Flags from Levelset on startup!')\n\
    flags.updateFromLevelset(density)\n\
  low_flags_updated = True \n\
  setWallBcs(flags=flags, vel=vel)\n\
  density.reinitMarching(flags=flags, velTransport=vel)\n\
  advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2)\n\
  flags.updateFromLevelset(density)\n\
  \n\
  advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2)\n\
  addGravity(flags=flags, vel=vel, gravity=vec3(0,0,-0.981))\n\
  \n\
  # mantaMsg current maximal velocity\n\
  maxvel = vel.getMaxValue()\n\
  mantaMsg('Current max velocity %f ' % maxvel)\n\
  \n\
  # pressure solve\n\
  setWallBcs(flags=flags, vel=vel)\n\
  solvePressure(flags=flags, vel=vel, pressure=pressure, cgMaxIterFac=0.5, useResNorm=True) \n\
  setWallBcs(flags=flags, vel=vel)\n\
  s.step()\n\
  density.multConst(-1.)\n";

//////////////////////////////////////////////////////////////////////
// EXPORT GRIDS
//////////////////////////////////////////////////////////////////////

const string smoke_export_low = "\n\
import os\n\
mantaMsg('Exporting grids low')\n\
density.save(os.path.join('$MANTA_EXPORT_PATH$','density.uni'))\n\
flags.save(os.path.join('$MANTA_EXPORT_PATH$','flags.uni'))\n\
vel.save(os.path.join('$MANTA_EXPORT_PATH$','vel.uni'))\n\
forces.save(os.path.join('$MANTA_EXPORT_PATH$','forces.uni'))\n\
inflow_grid.save(os.path.join('$MANTA_EXPORT_PATH$','inflow_low.uni'))\n\
fuel_inflow.save(os.path.join('$MANTA_EXPORT_PATH$','fuel_inflow.uni'))\n\
if using_colors:\n\
  color_r.save(os.path.join('$MANTA_EXPORT_PATH$','color_r.uni'))\n\
  color_g.save(os.path.join('$MANTA_EXPORT_PATH$','color_g.uni'))\n\
  color_b.save(os.path.join('$MANTA_EXPORT_PATH$','color_b.uni'))\n\
if using_heat:\n\
  heat.save(os.path.join('$MANTA_EXPORT_PATH$','heat.uni'))\n\
if using_fire:\n\
  flame.save(os.path.join('$MANTA_EXPORT_PATH$','flame.uni'))\n\
  fuel.save(os.path.join('$MANTA_EXPORT_PATH$','fuel.uni'))\n\
  react.save(os.path.join('$MANTA_EXPORT_PATH$','react.uni'))\n";

const string smoke_export_high = "\n\
mantaMsg('Exporting grids high')\n\
xl_density.save(os.path.join('$MANTA_EXPORT_PATH$','xl_density.uni'))\n\
xl_flags.save(os.path.join('$MANTA_EXPORT_PATH$','xl_flags.uni'))\n\
if using_colors:\n\
  xl_color_r.save(os.path.join('$MANTA_EXPORT_PATH$','xl_color_r.uni'))\n\
  xl_color_g.save(os.path.join('$MANTA_EXPORT_PATH$','xl_color_g.uni'))\n\
  xl_color_b.save(os.path.join('$MANTA_EXPORT_PATH$','xl_color_b.uni'))\n\
if using_fire:\n\
  xl_flame.save(os.path.join('$MANTA_EXPORT_PATH$','xl_flame.uni'))\n\
  xl_fuel.save(os.path.join('$MANTA_EXPORT_PATH$','xl_fuel.uni'))\n\
  xl_react.save(os.path.join('$MANTA_EXPORT_PATH$','xl_react.uni'))\n";

//////////////////////////////////////////////////////////////////////
// IMPORT GRIDS
//////////////////////////////////////////////////////////////////////

const string smoke_import_low = "\n\
def import_grids_low():\n\
  mantaMsg('Importing grids low')\n\
  density.load('$MANTA_EXPORT_PATH$density.uni')\n\
  flags.load('$MANTA_EXPORT_PATH$flags.uni')\n\
  vel.save(os.path.join('$MANTA_EXPORT_PATH$','vel.uni'))\n\
  forces.load('$MANTA_EXPORT_PATH$forces.uni')\n\
  inflow_grid.load('$MANTA_EXPORT_PATH$inflow_low.uni')\n\
  fuel_inflow.load('$MANTA_EXPORT_PATH$fuel_inflow.uni')\n\
  if using_colors:\n\
    color_r.load('$MANTA_EXPORT_PATH$color_r.uni')\n\
    color_g.load('$MANTA_EXPORT_PATH$color_g.uni')\n\
    color_b.load('$MANTA_EXPORT_PATH$color_b.uni')\n\
  if using_heat:\n\
    heat.load('$MANTA_EXPORT_PATH$heat.uni')\n\
  if using_fire:\n\
    flame.load('$MANTA_EXPORT_PATH$flame.uni')\n\
    fuel.load('$MANTA_EXPORT_PATH$fuel.uni')\n\
    react.load('$MANTA_EXPORT_PATH$react.uni')\n";

const string smoke_import_high = "\n\
def import_grids_high():\n\
  mantaMsg('Importing grids high')\n\
  xl_density.load('$MANTA_EXPORT_PATH$xl_density.uni')\n\
  xl_flags.load('$MANTA_EXPORT_PATH$xl_flags.uni')\n\
  if using_colors:\n\
    xl_color_r.load('$MANTA_EXPORT_PATH$xl_color_r.uni')\n\
    xl_color_g.load('$MANTA_EXPORT_PATH$xl_color_g.uni')\n\
    xl_color_b.load('$MANTA_EXPORT_PATH$xl_color_b.uni')\n\
  if using_fire:\n\
    xl_flame.load('$MANTA_EXPORT_PATH$xl_flame.uni')\n\
    xl_fuel.load('$MANTA_EXPORT_PATH$xl_fuel.uni')\n\
    xl_react.load('$MANTA_EXPORT_PATH$xl_react.uni')\n";

//////////////////////////////////////////////////////////////////////
// INFLOW
//////////////////////////////////////////////////////////////////////

const string smoke_inflow_low = "\n\
def apply_inflow():\n\
  mantaMsg('Applying inflow')\n\
  #inflow_grid.multConst(0.1)\n\
  #fuel_inflow.multConst(0.1)\n\
  density.add(inflow_grid)\n\
  if using_fire:\n\
    fuel.add(fuel_inflow)\n";

const string smoke_inflow_high = "\n\
 # TODO\n";
