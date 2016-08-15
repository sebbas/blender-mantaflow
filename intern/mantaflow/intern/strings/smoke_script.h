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

//////////////////////////////////////////////////////////////////////
// BOUNDS
//////////////////////////////////////////////////////////////////////

const std::string smoke_bounds_low = "\n\
# prepare domain low\n\
mantaMsg('Smoke domain low')\n\
flags.initDomain(boundaryWidth=boundaryWidth)\n\
flags.fillGrid()\n\
if doOpen:\n\
    setOpenBound(flags=flags, bWidth=boundaryWidth, openBound=boundConditions, type=FlagOutflow|FlagEmpty)\n";

const std::string smoke_bounds_high = "\n\
# prepare domain high\n\
mantaMsg('Smoke domain high')\n\
xl_flags.initDomain(boundaryWidth=boundaryWidth)\n\
xl_flags.fillGrid()\n\
if doOpen:\n\
    setOpenBound(flags=xl_flags, bWidth=boundaryWidth, openBound=boundConditions, type=FlagOutflow|FlagEmpty)\n";

//////////////////////////////////////////////////////////////////////
// VARIABLES
//////////////////////////////////////////////////////////////////////

const std::string smoke_variables_low = "\n\
mantaMsg('Smoke variables low')\n\
using_colors    = $USING_COLORS$\n\
using_heat      = $USING_HEAT$\n\
using_fire      = $USING_FIRE$\n\
using_wavelets  = $USE_WAVELETS$\n\
vorticity       = $VORTICITY$\n";

const std::string smoke_variables_high = "\n\
mantaMsg('Smoke variables high')\n\
wltStrength = $WLT_STR$\n\
octaves     = 0\n\
uvs         = 2\n\
\n\
if upres == 1:\n\
    octaves = int(math.log(upres+1)/ math.log(2.0) + 0.5)\n\
elif upres > 1:\n\
    octaves = int(math.log(upres)/ math.log(2.0) + 0.5)\n";

const std::string smoke_with_heat = "\n\
using_heat = True\n";

const std::string smoke_with_colors = "\n\
using_colors = True\n";

const std::string smoke_with_fire = "\n\
using_fire = True\n";

//////////////////////////////////////////////////////////////////////
// GRIDS
//////////////////////////////////////////////////////////////////////

const std::string smoke_alloc_low = "\n\
# prepare grids low\n\
mantaMsg('Smoke alloc grids low')\n\
flags       = s.create(FlagGrid)\n\
vel         = s.create(MACGrid)\n\
x_vel       = s.create(RealGrid)\n\
y_vel       = s.create(RealGrid)\n\
z_vel       = s.create(RealGrid)\n\
obvel       = s.create(MACGrid)\n\
x_obvel     = s.create(RealGrid)\n\
y_obvel     = s.create(RealGrid)\n\
z_obvel     = s.create(RealGrid)\n\
density     = s.create(LevelsetGrid)\n\
pressure    = s.create(RealGrid)\n\
forces      = s.create(MACGrid)\n\
x_force     = s.create(RealGrid)\n\
y_force     = s.create(RealGrid)\n\
z_force     = s.create(RealGrid)\n\
inflow_grid = s.create(LevelsetGrid)\n\
fuel_inflow = s.create(LevelsetGrid)\n";

const std::string smoke_alloc_high = "\n\
# prepare grids high\n\
mantaMsg('Smoke alloc grids high')\n\
xl_flags   = xl.create(FlagGrid)\n\
xl_vel     = xl.create(MACGrid)\n\
xl_density = xl.create(RealGrid)\n\
energy     = s.create(RealGrid)\n\
tempFlag   = s.create(FlagGrid)\n\
texture_u  = s.create(RealGrid)\n\
texture_v  = s.create(RealGrid)\n\
texture_w  = s.create(RealGrid)\n\
texture_u2 = s.create(RealGrid)\n\
texture_v2 = s.create(RealGrid)\n\
texture_w2 = s.create(RealGrid)\n";

//////////////////////////////////////////////////////////////////////
// ADDITIONAL GRIDS
//////////////////////////////////////////////////////////////////////

const std::string smoke_set_color_codes = "\n\
mantaMsg('Setting color codes')\n\
manta_color_r = $COLOR_R$\n\
manta_color_g = $COLOR_G$\n\
manta_color_b = $COLOR_B$\n";

const std::string smoke_alloc_colors_low = "\n\
mantaMsg('Allocating colors low')\n\
color_r = s.create(RealGrid)\n\
color_g = s.create(RealGrid)\n\
color_b = s.create(RealGrid)\n";

const std::string smoke_alloc_colors_high = "\
mantaMsg('Allocating colors high')\n\
xl_color_r = xl.create(RealGrid)\n\
xl_color_g = xl.create(RealGrid)\n\
xl_color_b = xl.create(RealGrid)\n";

const std::string smoke_init_colors_low = "\n\
mantaMsg('Initializing colors low')\n\
color_r.copyFrom(density) \n\
color_r.multConst(manta_color_r) \n\
color_g.copyFrom(density) \n\
color_g.multConst(manta_color_g) \n\
color_b.copyFrom(density) \n\
color_b.multConst(manta_color_b)\n";

const std::string smoke_init_colors_high = "\n\
mantaMsg('Initializing colors high')\n\
xl_color_r.copyFrom(xl_density) \n\
xl_color_r.multConst(manta_color_r) \n\
xl_color_g.copyFrom(xl_density) \n\
xl_color_g.multConst(manta_color_g) \n\
xl_color_b.copyFrom(xl_density) \n\
xl_color_b.multConst(manta_color_b)\n";

const std::string smoke_alloc_heat_low = "\n\
mantaMsg('Allocating heat low')\n\
heat = s.create(RealGrid)\n";

const std::string smoke_alloc_fire_low = "\n\
mantaMsg('Allocating fire low')\n\
flame = s.create(RealGrid)\n\
fuel  = s.create(RealGrid)\n\
react = s.create(RealGrid)\n";

const std::string smoke_alloc_fire_high = "\n\
mantaMsg('Allocating fire high')\n\
xl_flame = xl.create(RealGrid)\n\
xl_fuel  = xl.create(RealGrid)\n\
xl_react = xl.create(RealGrid)\n";

//////////////////////////////////////////////////////////////////////
// STEP FUNCTIONS
//////////////////////////////////////////////////////////////////////

const std::string smoke_adaptive_step = "\n\
def manta_step(start_frame):\n\
    s.frame = start_frame\n\
    s.timeTotal = s.frame * dt0\n\
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

const std::string smoke_step_low = "\n\
def step_low():\n\
    mantaMsg('Step low')\n\
    copyRealToVec3(sourceX=x_vel, sourceY=y_vel, sourceZ=z_vel, target=vel)\n\
    copyRealToVec3(sourceX=x_obvel, sourceY=y_obvel, sourceZ=z_obvel, target=obvel)\n\
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
    if doOpen:\n\
        resetOutflow(flags=flags, real=density)\n\
    \n\
    mantaMsg('Vorticity')\n\
    if vorticity > 0.01:\n\
        vorticityConfinement(vel=vel, flags=flags, strength=$VORTICITY$)\n\
    \n\
    if using_heat:\n\
        mantaMsg('Adding heat buoyancy')\n\
        addBuoyancy(flags=flags, density=density, vel=vel, gravity=gravity, coefficient=$ALPHA$)\n\
        addBuoyancy(flags=flags, density=heat, vel=vel, gravity=gravity, coefficient=$BETA$)\n\
    else:\n\
        mantaMsg('Adding buoyancy')\n\
        addBuoyancy(density=density, vel=vel, gravity=gravity, flags=flags)\n\
    \n\
    copyRealToVec3(sourceX=x_force, sourceY=y_force, sourceZ=z_force, target=forces)\n\
    mantaMsg('Adding forces')\n\
    addForceField(flags=flags, vel=vel, force=forces)\n\
    forces.clear()\n\
    \n\
    mantaMsg('Walls')\n\
    setWallBcs(flags=flags, vel=vel)\n\
    \n\
    mantaMsg('Pressure')\n\
    solvePressure(flags=flags, vel=vel, pressure=pressure)\n\
    \n\
    copyVec3ToReal(source=vel, targetX=x_vel, targetY=y_vel, targetZ=z_vel)\n\
    copyVec3ToReal(source=obvel, targetX=x_obvel, targetY=y_obvel, targetZ=z_obvel)\n\
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

const std::string smoke_step_high = "\n\
def step_high():\n\
    mantaMsg('Step high')\n\
    copyRealToVec3(sourceX=texture_u, sourceY=texture_v, sourceZ=texture_w, target=uv[0])\n\
    copyRealToVec3(sourceX=texture_u2, sourceY=texture_v2, sourceZ=texture_w2, target=uv[1])\n\
    \n\
    interpolateMACGrid(source=vel, target=xl_vel)\n\
    for i in range(uvs):\n\
        mantaMsg('Advecting UV')\n\
        advectSemiLagrange(flags=flags, vel=vel, grid=uv[i], order=$ADVECT_ORDER$)\n\
        mantaMsg('Updating UVWeight')\n\
        updateUvWeight(resetTime=10.0 , index=i, numUvs=uvs, uv=uv[i])\n\
    \n\
    mantaMsg('Energy')\n\
    computeEnergy(flags=flags, vel=vel, energy=energy)\n\
    \n\
    tempFlag.copyFrom(flags)\n\
    extrapolateSimpleFlags(flags=flags, val=tempFlag, distance=2, flagFrom=FlagObstacle, flagTo=FlagFluid)\n\
    extrapolateSimpleFlags(flags=tempFlag, val=energy, distance=6, flagFrom=FlagFluid, flagTo=FlagObstacle)\n\
    computeWaveletCoeffs(energy)\n\
    \n\
    sStr = 1.0 * wltStrength\n\
    sPos = 2.0\n\
    \n\
    mantaMsg('Applying noise vec')\n\
    for o in range(octaves):\n\
        for i in range(uvs):\n\
            uvWeight = getUvWeight(uv[i])\n\
            applyNoiseVec3(flags=xl_flags, target=xl_vel, noise=xl_wltnoise, scale=sStr * uvWeight, scaleSpatial=sPos , weight=energy, uv=uv[i])\n\
        sStr *= 0.06 # magic kolmogorov factor \n\
        sPos *= 2.0 \n\
    \n\
    for substep in range(upres):\n\
        if using_colors: \n\
            mantaMsg('Advecting colors high')\n\
            advectSemiLagrange(flags=xl_flags, vel=xl_vel, grid=xl_color_r, order=$ADVECT_ORDER$, openBounds=doOpen)\n\
            advectSemiLagrange(flags=xl_flags, vel=xl_vel, grid=xl_color_g, order=$ADVECT_ORDER$, openBounds=doOpen)\n\
            advectSemiLagrange(flags=xl_flags, vel=xl_vel, grid=xl_color_b, order=$ADVECT_ORDER$, openBounds=doOpen)\n\
        \n\
        if using_fire: \n\
            mantaMsg('Advecting fire high')\n\
            advectSemiLagrange(flags=xl_flags, vel=xl_vel, grid=xl_fuel, order=$ADVECT_ORDER$, openBounds=doOpen)\n\
            advectSemiLagrange(flags=xl_flags, vel=xl_vel, grid=xl_react, order=$ADVECT_ORDER$, openBounds=doOpen)\n\
        \n\
        mantaMsg('Advecting density high')\n\
        advectSemiLagrange(flags=xl_flags, vel=xl_vel, grid=xl_density, order=$ADVECT_ORDER$, openBounds=doOpen)\n\
    \n\
    copyVec3ToReal(source=uv[0], targetX=texture_u, targetY=texture_v, targetZ=texture_w)\n\
    copyVec3ToReal(source=uv[1], targetX=texture_u2, targetY=texture_v2, targetZ=texture_w2)\n\
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
// IMPORT / EXPORT
//////////////////////////////////////////////////////////////////////

const std::string smoke_import_low = "\n\
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

const std::string smoke_import_high = "\n\
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

const std::string smoke_export_low = "\n\
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

const std::string smoke_export_high = "\n\
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
// DESTRUCTION
//////////////////////////////////////////////////////////////////////

const std::string smoke_delete_colors_low = "\n\
mantaMsg('Deleting colors low')\n\
if 'color_r' in globals() : del color_r\n\
if 'color_g' in globals() : del color_g\n\
if 'color_b' in globals() : del color_b\n";

const std::string smoke_delete_colors_high = "\n\
mantaMsg('Deleting colors high')\n\
if 'xl_color_r' in globals() : del xl_color_r\n\
if 'xl_color_g' in globals() : del xl_color_g\n\
if 'xl_color_b' in globals() : del xl_color_b\n";

const std::string smoke_delete_fire_low = "\n\
mantaMsg('Deleting fire low')\n\
if 'flame' in globals() : del flame\n\
if 'fuel'  in globals() : del fuel\n\
if 'react' in globals() : del react\n";

const std::string smoke_delete_fire_high = "\n\
mantaMsg('Deleting fire high')\n\
if 'xl_flame' in globals() : del xl_flame\n\
if 'xl_fuel'  in globals() : del xl_fuel\n\
if 'xl_react' in globals() : del xl_react\n";

const std::string smoke_delete_heat_low = "\n\
mantaMsg('Deleting heat low')\n\
if 'heat' in globals() : del heat\n";

const std::string smoke_delete_grids_low = "\n\
mantaMsg('Deleting base grids low')\n\
if 'flags'       in globals() : del flags\n\
if 'vel'         in globals() : del vel\n\
if 'x_vel'       in globals() : del x_vel\n\
if 'y_vel'       in globals() : del y_vel\n\
if 'z_vel'       in globals() : del z_vel\n\
if 'obvel'       in globals() : del obvel\n\
if 'x_obvel'     in globals() : del x_obvel\n\
if 'y_obvel'     in globals() : del y_obvel\n\
if 'z_obvel'     in globals() : del z_obvel\n\
if 'density'     in globals() : del density\n\
if 'pressure'    in globals() : del pressure\n\
if 'forces'      in globals() : del forces\n\
if 'x_force'     in globals() : del x_force\n\
if 'y_force'     in globals() : del y_force\n\
if 'z_force'     in globals() : del z_force\n\
if 'inflow_grid' in globals() : del inflow_grid\n\
if 'fuel_inflow' in globals() : del fuel_inflow\n";

const std::string smoke_delete_grids_high = "\n\
mantaMsg('Deleting base grids high')\n\
if 'xl_flags'    in globals() : del xl_flags\n\
if 'xl_vel'      in globals() : del xl_vel\n\
if 'xl_density'  in globals() : del xl_density\n\
if 'energy'      in globals() : del energy\n\
if 'tempFlag'    in globals() : del tempFlag\n\
if 'uvGrid'      in globals() : del uvGrid\n\
if 'texture_u'   in globals() : del texture_u\n\
if 'texture_v'   in globals() : del texture_v\n\
if 'texture_w'   in globals() : del texture_w\n\
if 'texture_u2'  in globals() : del texture_u2\n\
if 'texture_v2'  in globals() : del texture_v2\n\
if 'texture_w2'  in globals() : del texture_w2\n\
if 'xl_wltnoise' in globals() : del xl_wltnoise\n";

const std::string smoke_delete_variables_low = "\n\
mantaMsg('Deleting variables low')\n\
if 'res'             in globals() : del res\n\
if 'dim'             in globals() : del dim\n\
if 'gs'              in globals() : del gs\n\
if 'doOpen'          in globals() : del doOpen\n\
if 'boundConditions' in globals() : del boundConditions\n\
if 'dt_default'      in globals() : del dt_default\n\
if 'dt_factor'       in globals() : del dt_factor\n\
if 'fps'             in globals() : del fps\n\
if 'dt0'             in globals() : del dt0\n\
if 'vorticity'       in globals() : del vorticity\n\
if 'boundaryWidth'   in globals() : del boundaryWidth\n\
if 'using_colors'    in globals() : del using_colors\n\
if 'using_heat'      in globals() : del using_heat\n\
if 'using_fire'      in globals() : del using_fire\n\
if 'last_frame'      in globals() : del last_frame\n\
if 'maxvel'          in globals() : del maxvel\n";

const std::string smoke_delete_variables_high = "\n\
mantaMsg('Deleting variables high')\n\
if 'upres'       in globals() : del upres\n\
if 'xl_gs'       in globals() : del xl_gs\n\
if 'wltStrength' in globals() : del wltStrength\n\
if 'uvs'         in globals() : del uvs\n\
if 'uv'          in globals() : del uv\n\
if 'octaves'     in globals() : del octaves\n";

//////////////////////////////////////////////////////////////////////
// OTHER SETUPS
//////////////////////////////////////////////////////////////////////

const std::string smoke_uv_setup = "\n\
# create the array of uv grids\n\
uv = []\n\
mantaMsg('Initializing UV Grids')\n\
for i in range(uvs):\n\
    uvGrid = s.create(VecGrid)\n\
    uv.append(uvGrid)\n\
    resetUvGrid(uv[i])\n\
\n\
# Need to initialize helper grids for uvw as well\n\
copyVec3ToReal(source=uv[0], targetX=texture_u, targetY=texture_v, targetZ=texture_w)\n\
copyVec3ToReal(source=uv[1], targetX=texture_u2, targetY=texture_v2, targetZ=texture_w2)\n";

const std::string smoke_wavelet_turbulence_noise = "\n\
# wavelet turbulence noise field\n\
mantaMsg('Smoke wavelet noise')\n\
xl_wltnoise = NoiseField(parent=xl, loadFromFile=True)\n\
xl_wltnoise.posScale = vec3(int(1.0*gs.x)) / $NOISE_POSSCALE$\n\
xl_wltnoise.timeAnim = $NOISE_TIMEANIM$\n";

const std::string smoke_inflow_low = "\n\
def apply_inflow():\n\
    mantaMsg('Applying inflow')\n\
    #inflow_grid.multConst(0.1)\n\
    #fuel_inflow.multConst(0.1)\n\
    density.add(inflow_grid)\n\
    if using_fire:\n\
        fuel.add(fuel_inflow)\n";

const std::string smoke_inflow_high = "\n\
# TODO\n";

//////////////////////////////////////////////////////////////////////
// STANDALONE MODE
//////////////////////////////////////////////////////////////////////

const std::string smoke_standalone = "\n\
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
start_frame = $CURRENT_FRAME$\n\
# All low and high res steps\n\
manta_step(start_frame)\n";
