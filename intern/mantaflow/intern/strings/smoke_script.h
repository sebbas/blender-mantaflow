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
flags_s$ID$.initDomain(boundaryWidth=boundaryWidth_s$ID$)\n\
flags_s$ID$.fillGrid()\n\
if doOpen_s$ID$:\n\
    setOpenBound(flags=flags_s$ID$, bWidth=boundaryWidth_s$ID$, openBound=boundConditions_s$ID$, type=FlagOutflow|FlagEmpty)\n";

const std::string smoke_bounds_high = "\n\
# prepare domain high\n\
mantaMsg('Smoke domain high')\n\
flags_xl$ID$.initDomain(boundaryWidth=boundaryWidth_s$ID$)\n\
flags_xl$ID$.fillGrid()\n\
if doOpen_s$ID$:\n\
    setOpenBound(flags=flags_xl$ID$, bWidth=boundaryWidth_s$ID$, openBound=boundConditions_s$ID$, type=FlagOutflow|FlagEmpty)\n";

//////////////////////////////////////////////////////////////////////
// VARIABLES
//////////////////////////////////////////////////////////////////////

const std::string smoke_variables_low = "\n\
mantaMsg('Smoke variables low')\n\
preconditioner_s$ID$  = PcMGStatic\n\
using_colors_s$ID$    = $USING_COLORS$\n\
using_heat_s$ID$      = $USING_HEAT$\n\
using_fire_s$ID$      = $USING_FIRE$\n\
vorticity_s$ID$       = $VORTICITY$\n";

const std::string smoke_variables_high = "\n\
mantaMsg('Smoke variables high')\n\
wltStrength_s$ID$ = $WLT_STR$\n\
octaves_s$ID$     = 0\n\
uvs_s$ID$         = 2\n\
\n\
if upres_xl$ID$ == 1:\n\
    octaves_s$ID$ = int(math.log(upres_xl$ID$+1)/ math.log(2.0) + 0.5)\n\
elif upres_xl$ID$ > 1:\n\
    octaves_s$ID$ = int(math.log(upres_xl$ID$)/ math.log(2.0) + 0.5)\n";

const std::string smoke_with_heat = "\n\
using_heat_s$ID$ = True\n";

const std::string smoke_with_colors = "\n\
using_colors_s$ID$ = True\n";

const std::string smoke_with_fire = "\n\
using_fire_s$ID$ = True\n";

//////////////////////////////////////////////////////////////////////
// GRIDS
//////////////////////////////////////////////////////////////////////

const std::string smoke_alloc_low = "\n\
# prepare grids low\n\
mantaMsg('Smoke alloc low')\n\
flags_s$ID$       = s$ID$.create(FlagGrid)\n\
vel_s$ID$         = s$ID$.create(MACGrid)\n\
x_vel_s$ID$       = s$ID$.create(RealGrid)\n\
y_vel_s$ID$       = s$ID$.create(RealGrid)\n\
z_vel_s$ID$       = s$ID$.create(RealGrid)\n\
invel_s$ID$       = s$ID$.create(VecGrid)\n\
x_invel_s$ID$     = s$ID$.create(RealGrid)\n\
y_invel_s$ID$     = s$ID$.create(RealGrid)\n\
z_invel_s$ID$     = s$ID$.create(RealGrid)\n\
density_s$ID$     = s$ID$.create(RealGrid)\n\
pressure_s$ID$    = s$ID$.create(RealGrid)\n\
phiObs_s$ID$      = s$ID$.create(LevelsetGrid)\n\
phiOutIn_s$ID$    = s$ID$.create(LevelsetGrid)\n\
forces_s$ID$      = s$ID$.create(Vec3Grid)\n\
x_force_s$ID$     = s$ID$.create(RealGrid)\n\
y_force_s$ID$     = s$ID$.create(RealGrid)\n\
z_force_s$ID$     = s$ID$.create(RealGrid)\n\
inflow_s$ID$      = s$ID$.create(RealGrid)\n";

const std::string smoke_alloc_high = "\n\
# prepare grids high\n\
mantaMsg('Smoke alloc high')\n\
flags_xl$ID$     = xl$ID$.create(FlagGrid)\n\
vel_xl$ID$       = xl$ID$.create(MACGrid)\n\
density_xl$ID$   = xl$ID$.create(RealGrid)\n\
phiObs_xl$ID$    = xl$ID$.create(LevelsetGrid)\n\
energy_s$ID$     = s$ID$.create(RealGrid)\n\
tempFlag_s$ID$   = s$ID$.create(FlagGrid)\n\
texture_u_s$ID$  = s$ID$.create(RealGrid)\n\
texture_v_s$ID$  = s$ID$.create(RealGrid)\n\
texture_w_s$ID$  = s$ID$.create(RealGrid)\n\
texture_u2_s$ID$ = s$ID$.create(RealGrid)\n\
texture_v2_s$ID$ = s$ID$.create(RealGrid)\n\
texture_w2_s$ID$ = s$ID$.create(RealGrid)\n";

//////////////////////////////////////////////////////////////////////
// ADDITIONAL GRIDS
//////////////////////////////////////////////////////////////////////

const std::string smoke_alloc_colors_low = "\n\
mantaMsg('Allocating colors low')\n\
color_r_s$ID$ = s$ID$.create(RealGrid)\n\
color_g_s$ID$ = s$ID$.create(RealGrid)\n\
color_b_s$ID$ = s$ID$.create(RealGrid)\n";

const std::string smoke_alloc_colors_high = "\
mantaMsg('Allocating colors high')\n\
color_r_xl$ID$ = xl$ID$.create(RealGrid)\n\
color_g_xl$ID$ = xl$ID$.create(RealGrid)\n\
color_b_xl$ID$ = xl$ID$.create(RealGrid)\n";

const std::string smoke_init_colors_low = "\n\
mantaMsg('Initializing colors low')\n\
color_r_s$ID$.copyFrom(density_s$ID$) \n\
color_r_s$ID$.multConst($COLOR_R$) \n\
color_g_s$ID$.copyFrom(density_s$ID$) \n\
color_g_s$ID$.multConst($COLOR_G$) \n\
color_b_s$ID$.copyFrom(density_s$ID$) \n\
color_b_s$ID$.multConst($COLOR_B$)\n";

const std::string smoke_init_colors_high = "\n\
mantaMsg('Initializing colors high')\n\
color_r_xl$ID$.copyFrom(density_xl$ID$) \n\
color_r_xl$ID$.multConst($COLOR_R$) \n\
color_g_xl$ID$.copyFrom(density_xl$ID$) \n\
color_g_xl$ID$.multConst($COLOR_G$) \n\
color_b_xl$ID$.copyFrom(density_xl$ID$) \n\
color_b_xl$ID$.multConst($COLOR_B$)\n";

const std::string smoke_alloc_heat_low = "\n\
mantaMsg('Allocating heat low')\n\
heat_s$ID$ = s$ID$.create(RealGrid)\n";

const std::string smoke_alloc_fire_low = "\n\
mantaMsg('Allocating fire low')\n\
flame_s$ID$ = s$ID$.create(RealGrid)\n\
fuel_s$ID$  = s$ID$.create(RealGrid)\n\
react_s$ID$ = s$ID$.create(RealGrid)\n";

const std::string smoke_alloc_fire_high = "\n\
mantaMsg('Allocating fire high')\n\
flame_xl$ID$ = xl$ID$.create(RealGrid)\n\
fuel_xl$ID$  = xl$ID$.create(RealGrid)\n\
react_xl$ID$ = xl$ID$.create(RealGrid)\n";

//////////////////////////////////////////////////////////////////////
// PRE / POST STEP
//////////////////////////////////////////////////////////////////////

const std::string smoke_pre_step_low = "\n\
def smoke_pre_step_low_$ID$():\n\
    mantaMsg('Smoke pre step low')\n\
    # translate obvels (world space) to grid space\n\
    if using_obstacle_s$ID$:\n\
        x_obvel_s$ID$.multConst(gs_s$ID$.x)\n\
        y_obvel_s$ID$.multConst(gs_s$ID$.y)\n\
        z_obvel_s$ID$.multConst(gs_s$ID$.z)\n\
    \n\
    # translate guiding vels (world space) to grid space\n\
    if using_guiding_s$ID$:\n\
        x_guidevel_s$ID$.multConst(gs_s$ID$.x)\n\
        y_guidevel_s$ID$.multConst(gs_s$ID$.y)\n\
        z_guidevel_s$ID$.multConst(gs_s$ID$.z)\n\
        copyRealToVec3(sourceX=x_guidevel_s$ID$, sourceY=y_guidevel_s$ID$, sourceZ=z_guidevel_s$ID$, target=guidevelC_s$ID$)\n\
    \n\
    # translate invels (world space) to grid space\n\
    x_invel_s$ID$.multConst(gs_s$ID$.x)\n\
    y_invel_s$ID$.multConst(gs_s$ID$.y)\n\
    z_invel_s$ID$.multConst(gs_s$ID$.z)\n\
    \n\
    copyRealToVec3(sourceX=x_vel_s$ID$, sourceY=y_vel_s$ID$, sourceZ=z_vel_s$ID$, target=vel_s$ID$)\n\
    copyRealToVec3(sourceX=x_invel_s$ID$, sourceY=y_invel_s$ID$, sourceZ=z_invel_s$ID$, target=invel_s$ID$)\n\
    copyRealToVec3(sourceX=x_force_s$ID$, sourceY=y_force_s$ID$, sourceZ=z_force_s$ID$, target=forces_s$ID$)\n\
    \n\
    if using_obstacle_s$ID$:\n\
        copyRealToVec3(sourceX=x_obvel_s$ID$, sourceY=y_obvel_s$ID$, sourceZ=z_obvel_s$ID$, target=obvelC_s$ID$)\n\
    \n\
    # If obstacle has velocity, i.e. is moving switch to dynamic preconditioner\n\
    if using_obstacle_s$ID$ and obvelC_s$ID$.getMaxValue() > 0:\n\
        mantaMsg('Using dynamic preconditioner')\n\
        preconditioner_s$ID$ = PcMGDynamic\n\
    else:\n\
        mantaMsg('Using static preconditioner')\n\
        preconditioner_s$ID$ = PcMGStatic\n";

const std::string smoke_pre_step_high = "\n\
def smoke_pre_step_high_$ID$():\n\
    copyRealToVec3(sourceX=texture_u_s$ID$, sourceY=texture_v_s$ID$, sourceZ=texture_w_s$ID$, target=uv_s$ID$[0])\n\
    copyRealToVec3(sourceX=texture_u2_s$ID$, sourceY=texture_v2_s$ID$, sourceZ=texture_w2_s$ID$, target=uv_s$ID$[1])\n";

const std::string smoke_post_step_low = "\n\
def smoke_post_step_low_$ID$():\n\
    forces_s$ID$.clear()\n\
    invel_s$ID$.clear()\n\
    phiOutIn_s$ID$.setConst(9999)\n\
    \n\
    if using_obstacle_s$ID$:\n\
        phiObsIn_s$ID$.setConst(9999)\n\
        obvelC_s$ID$.clear()\n\
    \n\
    if using_guiding_s$ID$:\n\
        phiGuideIn_s$ID$.setConst(9999)\n\
        guidevelC_s$ID$.clear()\n\
    copyVec3ToReal(source=vel_s$ID$, targetX=x_vel_s$ID$, targetY=y_vel_s$ID$, targetZ=z_vel_s$ID$)\n";

const std::string smoke_post_step_high = "\n\
def smoke_post_step_high_$ID$():\n\
    copyVec3ToReal(source=uv_s$ID$[0], targetX=texture_u_s$ID$, targetY=texture_v_s$ID$, targetZ=texture_w_s$ID$)\n\
    copyVec3ToReal(source=uv_s$ID$[1], targetX=texture_u2_s$ID$, targetY=texture_v2_s$ID$, targetZ=texture_w2_s$ID$)\n";

//////////////////////////////////////////////////////////////////////
// STEP FUNCTIONS
//////////////////////////////////////////////////////////////////////

const std::string smoke_adaptive_step = "\n\
def manta_step_$ID$(framenr):\n\
    s$ID$.frame = framenr\n\
    s$ID$.timeTotal = s$ID$.frame * dt0_s$ID$\n\
    last_frame_s$ID$ = s$ID$.frame\n\
    \n\
    smoke_pre_step_low_$ID$()\n\
    if using_highres_s$ID$:\n\
        smoke_pre_step_high_$ID$()\n\
    \n\
    while s$ID$.frame == last_frame_s$ID$:\n\
        \n\
        if using_obstacle_s$ID$: # TODO (sebbas): allow outflow objects when no obstacle set\n\
            phiObs_s$ID$.join(phiObsIn_s$ID$)\n\
        setObstacleFlags(flags=flags_s$ID$, phiObs=phiObs_s$ID$, phiOut=phiOutIn_s$ID$)\n\
        flags_s$ID$.fillGrid()\n\
        \n\
        if using_adaptTime_s$ID$:\n\
            mantaMsg('Adapt timestep')\n\
            maxvel_s$ID$ = vel_s$ID$.getMaxValue()\n\
            s$ID$.adaptTimestep(maxvel_s$ID$)\n\
        \n\
        mantaMsg('Low step / s$ID$.frame: ' + str(s$ID$.frame))\n\
        if using_fire_s$ID$:\n\
            process_burn_low_$ID$()\n\
        step_low_$ID$()\n\
        if using_fire_s$ID$:\n\
            update_flame_low_$ID$()\n\
        \n\
        if using_highres_s$ID$:\n\
            xl$ID$.timestep = s$ID$.timestep\n\
            mantaMsg('High step / s$ID$.frame: ' + str(s$ID$.frame))\n\
            if using_fire_s$ID$:\n\
                process_burn_high_$ID$()\n\
            step_high_$ID$()\n\
            if using_fire_s$ID$:\n\
                update_flame_high_$ID$()\n\
        s$ID$.step()\n\
    \n\
    smoke_post_step_low_$ID$()\n\
    if using_highres_s$ID$:\n\
        smoke_post_step_high_$ID$()\n";

const std::string smoke_step_low = "\n\
def step_low_$ID$():\n\
    mantaMsg('Smoke step low')\n\
    \n\
    mantaMsg('Advecting density')\n\
    advectSemiLagrange(flags=flags_s$ID$, vel=vel_s$ID$, grid=density_s$ID$, order=$ADVECT_ORDER$)\n\
    \n\
    if using_heat_s$ID$:\n\
        mantaMsg('Advecting heat')\n\
        advectSemiLagrange(flags=flags_s$ID$, vel=vel_s$ID$, grid=heat_s$ID$, order=$ADVECT_ORDER$)\n\
    \n\
    if using_fire_s$ID$:\n\
        mantaMsg('Advecting fire')\n\
        advectSemiLagrange(flags=flags_s$ID$, vel=vel_s$ID$, grid=fuel_s$ID$, order=$ADVECT_ORDER$)\n\
        advectSemiLagrange(flags=flags_s$ID$, vel=vel_s$ID$, grid=react_s$ID$, order=$ADVECT_ORDER$)\n\
    \n\
    if using_colors_s$ID$:\n\
        mantaMsg('Advecting colors')\n\
        advectSemiLagrange(flags=flags_s$ID$, vel=vel_s$ID$, grid=color_r_s$ID$, order=$ADVECT_ORDER$)\n\
        advectSemiLagrange(flags=flags_s$ID$, vel=vel_s$ID$, grid=color_g_s$ID$, order=$ADVECT_ORDER$)\n\
        advectSemiLagrange(flags=flags_s$ID$, vel=vel_s$ID$, grid=color_b_s$ID$, order=$ADVECT_ORDER$)\n\
    \n\
    mantaMsg('Advecting velocity')\n\
    advectSemiLagrange(flags=flags_s$ID$, vel=vel_s$ID$, grid=vel_s$ID$, order=$ADVECT_ORDER$, openBounds=doOpen_s$ID$, boundaryWidth=boundaryWidth_s$ID$)\n\
    \n\
    # Create interpolated version of original phi grid for later use in (optional) high-res step\n\
    if using_obstacle_s$ID$ and using_highres_s$ID$:\n\
        interpolateGrid(target=phiObs_xl$ID$, source=phiObs_s$ID$)\n\
    if doOpen_s$ID$:\n\
        resetOutflow(flags=flags_s$ID$, real=density_s$ID$)\n\
    \n\
    mantaMsg('Vorticity')\n\
    vorticityConfinement(vel=vel_s$ID$, flags=flags_s$ID$, strength=$VORTICITY$)\n\
    \n\
    if using_heat_s$ID$:\n\
        mantaMsg('Adding heat buoyancy')\n\
        addBuoyancy(flags=flags_s$ID$, density=density_s$ID$, vel=vel_s$ID$, gravity=gravity_s$ID$, coefficient=$ALPHA$)\n\
        addBuoyancy(flags=flags_s$ID$, density=heat_s$ID$, vel=vel_s$ID$, gravity=gravity_s$ID$, coefficient=$BETA$)\n\
    else:\n\
        mantaMsg('Adding buoyancy')\n\
        addBuoyancy(density=density_s$ID$, vel=vel_s$ID$, gravity=gravity_s$ID$, flags=flags_s$ID$)\n\
    \n\
    mantaMsg('Adding forces')\n\
    addForceField(flags=flags_s$ID$, vel=vel_s$ID$, force=forces_s$ID$)\n\
    \n\
    # add initial velocity\n\
    setInitialVelocity(flags=flags_s$ID$, vel=vel_s$ID$, invel=invel_s$ID$)\n\
    \n\
    if using_obstacle_s$ID$:\n\
        mantaMsg('Extrapolating object velocity')\n\
        # ensure velocities inside of obs object, slightly add obvels outside of obs object\n\
        extrapolateVec3Simple(vel=obvelC_s$ID$, phi=phiObsIn_s$ID$, distance=int(res_s$ID$/2), inside=True)\n\
        extrapolateVec3Simple(vel=obvelC_s$ID$, phi=phiObsIn_s$ID$, distance=1, inside=False)\n\
        resampleVec3ToMac(source=obvelC_s$ID$, target=obvel_s$ID$)\n\
    \n\
    if using_guiding_s$ID$:\n\
        mantaMsg('Extrapolating guiding velocity')\n\
        # ensure velocities inside of guiding object, slightly add guiding vels outside of object too\n\
        extrapolateVec3Simple(vel=guidevelC_s$ID$, phi=phiGuideIn_s$ID$, distance=int(res_s$ID$/2), inside=True)\n\
        extrapolateVec3Simple(vel=guidevelC_s$ID$, phi=phiGuideIn_s$ID$, distance=1, inside=False)\n\
        resampleVec3ToMac(source=guidevelC_s$ID$, target=guidevel_s$ID$)\n\
    \n\
    mantaMsg('Walls')\n\
    if using_obstacle_s$ID$:\n\
        setWallBcs(flags=flags_s$ID$, vel=vel_s$ID$, obvel=obvel_s$ID$)\n\
    else:\n\
        setWallBcs(flags=flags_s$ID$, vel=vel_s$ID$)\n\
    \n\
    if using_guiding_s$ID$:\n\
        mantaMsg('Guiding and pressure')\n\
        weightGuide_s$ID$.multConst(0)\n\
        weightGuide_s$ID$.addConst(alpha_s$ID$)\n\
        PD_fluid_guiding(vel=vel_s$ID$, velT=guidevel_s$ID$, flags=flags_s$ID$, weight=weightGuide_s$ID$, blurRadius=beta_s$ID$, pressure=pressure_s$ID$, tau=tau_s$ID$, sigma=sigma_s$ID$, theta=theta_s$ID$, preconditioner=preconditioner_s$ID$, zeroPressureFixing=not doOpen_s$ID$)\n\
    else:\n\
        mantaMsg('Pressure')\n\
        solvePressure(flags=flags_s$ID$, vel=vel_s$ID$, pressure=pressure_s$ID$, preconditioner=preconditioner_s$ID$, zeroPressureFixing=not doOpen_s$ID$) # closed domains require pressure fixing\n\
\n\
def process_burn_low_$ID$():\n\
    mantaMsg('Process burn low')\n\
    if (using_colors_s$ID$):\n\
        processBurn(fuel=fuel_s$ID$, density=density_s$ID$, react=react_s$ID$, red=color_r_s$ID$, green=color_g_s$ID$, blue=color_b_s$ID$, heat=heat_s$ID$, burningRate=$BURNING_RATE$, flameSmoke=$FLAME_SMOKE$, ignitionTemp=$IGNITION_TEMP$, maxTemp=$MAX_TEMP$, dt=dt0_s$ID$, flameSmokeColor=vec3($FLAME_SMOKE_COLOR_X$,$FLAME_SMOKE_COLOR_Y$,$FLAME_SMOKE_COLOR_Z$))\n\
    else:\n\
        processBurn(fuel=fuel_s$ID$, density=density_s$ID$, react=react_s$ID$, heat=heat_s$ID$, burningRate=$BURNING_RATE$, flameSmoke=$FLAME_SMOKE$, ignitionTemp=$IGNITION_TEMP$, maxTemp=$MAX_TEMP$, dt=dt0_s$ID$, flameSmokeColor=vec3($FLAME_SMOKE_COLOR_X$,$FLAME_SMOKE_COLOR_Y$,$FLAME_SMOKE_COLOR_Z$))\n\
\n\
def update_flame_low_$ID$():\n\
    mantaMsg('Update flame low')\n\
    updateFlame(react=react_s$ID$, flame=flame_s$ID$)\n";

const std::string smoke_step_high = "\n\
def step_high_$ID$():\n\
    mantaMsg('Smoke step high')\n\
    if using_obstacle_s$ID$: # TODO (sebbas): allow outflow objects when no obstacle set\n\
        setObstacleFlags(flags=flags_xl$ID$, phiObs=phiObs_xl$ID$)\n\
    flags_xl$ID$.fillGrid()\n\
    \n\
    interpolateMACGrid(source=vel_s$ID$, target=vel_xl$ID$)\n\
    for i in range(uvs_s$ID$):\n\
        mantaMsg('Advecting UV')\n\
        advectSemiLagrange(flags=flags_s$ID$, vel=vel_s$ID$, grid=uv_s$ID$[i], order=$ADVECT_ORDER$)\n\
        mantaMsg('Updating UVWeight')\n\
        updateUvWeight(resetTime=10.0 , index=i, numUvs=uvs_s$ID$, uv=uv_s$ID$[i])\n\
    \n\
    mantaMsg('Energy')\n\
    computeEnergy(flags=flags_s$ID$, vel=vel_s$ID$, energy=energy_s$ID$)\n\
    \n\
    tempFlag_s$ID$.copyFrom(flags_s$ID$)\n\
    extrapolateSimpleFlags(flags=flags_s$ID$, val=tempFlag_s$ID$, distance=2, flagFrom=FlagObstacle, flagTo=FlagFluid)\n\
    extrapolateSimpleFlags(flags=tempFlag_s$ID$, val=energy_s$ID$, distance=6, flagFrom=FlagFluid, flagTo=FlagObstacle)\n\
    computeWaveletCoeffs(energy_s$ID$)\n\
    \n\
    sStr_s$ID$ = 1.0 * wltStrength_s$ID$\n\
    sPos_s$ID$ = 2.0\n\
    \n\
    mantaMsg('Applying noise vec')\n\
    for o in range(octaves_s$ID$):\n\
        for i in range(uvs_s$ID$):\n\
            uvWeight_s$ID$ = getUvWeight(uv_s$ID$[i])\n\
            applyNoiseVec3(flags=flags_xl$ID$, target=vel_xl$ID$, noise=wltnoise_xl$ID$, scale=sStr_s$ID$ * uvWeight_s$ID$, scaleSpatial=sPos_s$ID$ , weight=energy_s$ID$, uv=uv_s$ID$[i])\n\
        sStr_s$ID$ *= 0.06 # magic kolmogorov factor \n\
        sPos_s$ID$ *= 2.0 \n\
    \n\
    for substep in range(upres_xl$ID$):\n\
        if using_colors_s$ID$: \n\
            mantaMsg('Advecting colors high')\n\
            advectSemiLagrange(flags=flags_xl$ID$, vel=vel_xl$ID$, grid=color_r_xl$ID$, order=$ADVECT_ORDER$, openBounds=doOpen_s$ID$)\n\
            advectSemiLagrange(flags=flags_xl$ID$, vel=vel_xl$ID$, grid=color_g_xl$ID$, order=$ADVECT_ORDER$, openBounds=doOpen_s$ID$)\n\
            advectSemiLagrange(flags=flags_xl$ID$, vel=vel_xl$ID$, grid=color_b_xl$ID$, order=$ADVECT_ORDER$, openBounds=doOpen_s$ID$)\n\
        \n\
        if using_fire_s$ID$: \n\
            mantaMsg('Advecting fire high')\n\
            advectSemiLagrange(flags=flags_xl$ID$, vel=vel_xl$ID$, grid=fuel_xl$ID$, order=$ADVECT_ORDER$, openBounds=doOpen_s$ID$)\n\
            advectSemiLagrange(flags=flags_xl$ID$, vel=vel_xl$ID$, grid=react_xl$ID$, order=$ADVECT_ORDER$, openBounds=doOpen_s$ID$)\n\
        \n\
        mantaMsg('Advecting density high')\n\
        advectSemiLagrange(flags=flags_xl$ID$, vel=vel_xl$ID$, grid=density_xl$ID$, order=$ADVECT_ORDER$, openBounds=doOpen_s$ID$)\n\
\n\
def process_burn_high_$ID$():\n\
    mantaMsg('Process burn high')\n\
    if (using_colors_s$ID$):\n\
        processBurn(fuel=fuel_xl$ID$, density=density_xl$ID$, react=react_xl$ID$, red=color_r_xl$ID$, green=color_g_xl$ID$, blue=color_b_xl$ID$, burningRate=$BURNING_RATE$, flameSmoke=$FLAME_SMOKE$, ignitionTemp=$IGNITION_TEMP$, maxTemp=$MAX_TEMP$, dt=dt0_s$ID$, flameSmokeColor=vec3($FLAME_SMOKE_COLOR_X$,$FLAME_SMOKE_COLOR_Y$,$FLAME_SMOKE_COLOR_Z$))\n\
    else:\n\
        processBurn(fuel=fuel_xl$ID$, density=density_xl$ID$, react=react_xl$ID$, burningRate=$BURNING_RATE$, flameSmoke=$FLAME_SMOKE$, ignitionTemp=$IGNITION_TEMP$, maxTemp=$MAX_TEMP$, dt=dt0_s$ID$, flameSmokeColor=vec3($FLAME_SMOKE_COLOR_X$,$FLAME_SMOKE_COLOR_Y$,$FLAME_SMOKE_COLOR_Z$))\n\
\n\
def update_flame_high_$ID$():\n\
    mantaMsg('Update flame high')\n\
    updateFlame(react=react_xl$ID$, flame=flame_xl$ID$)\n";

//////////////////////////////////////////////////////////////////////
// IMPORT / EXPORT
//////////////////////////////////////////////////////////////////////

const std::string smoke_import_low = "\n\
def load_smoke_data_low_$ID$(path):\n\
    density_s$ID$.load(path + '_density.uni')\n\
    flags_s$ID$.load(path + '_flags.uni')\n\
    vel_s$ID$.load(path + '_vel.uni')\n\
    invel_s$ID$.load(path + '_invel.uni')\n\
    pressure_s$ID$.load(path + '_pressure.uni')\n\
    forces_s$ID$.load(path + '_forces.uni')\n\
    x_force_s$ID$.load(path + '_x_force.uni')\n\
    y_force_s$ID$.load(path + '_y_force.uni')\n\
    z_force_s$ID$.load(path + '_z_force.uni')\n\
    inflow_s$ID$.load(path + '_inflow.uni')\n\
    x_vel_s$ID$.load(path + '_x_vel.uni')\n\
    y_vel_s$ID$.load(path + '_y_vel.uni')\n\
    z_vel_s$ID$.load(path + '_z_vel.uni')\n\
    x_invel_s$ID$.load(path + '_x_invel.uni')\n\
    y_invel_s$ID$.load(path + '_y_invel.uni')\n\
    z_invel_s$ID$.load(path + '_z_invel.uni')\n\
    phiObs_s$ID$.load(path + '_phiObs.uni')\n\
    phiOutIn_s$ID$.load(path + '_phiOutIn.uni')\n\
    if using_colors_s$ID$:\n\
        color_r_s$ID$.load(path + '_color_r.uni')\n\
        color_g_s$ID$.load(path + '_color_g.uni')\n\
        color_b_s$ID$.load(path + '_color_b.uni')\n\
    if using_heat_s$ID$:\n\
        heat_s$ID$.load(path + '_heat.uni')\n\
    if using_fire_s$ID$:\n\
        flame_s$ID$.load(path + '_flame.uni')\n\
        fuel_s$ID$.load(path + '_fuel.uni')\n\
        react_s$ID$.load(path + '_react.uni')\n";

const std::string smoke_import_high = "\n\
def load_smoke_data_high_$ID$(path):\n\
    density_xl$ID$.load(path + '_density_xl.uni')\n\
    flags_xl$ID$.load(path + '_flags_xl.uni')\n\
    \n\
    texture_u_s$ID$.load(path + '_texture_u.uni')\n\
    texture_v_s$ID$.load(path + '_texture_v.uni')\n\
    texture_w_s$ID$.load(path + '_texture_w.uni')\n\
    texture_u2_s$ID$.load(path + '_texture_u2.uni')\n\
    texture_v2_s$ID$.load(path + '_texture_v2.uni')\n\
    texture_w2_s$ID$.load(path + '_texture_w2.uni')\n\
    \n\
    if using_colors_s$ID$:\n\
        color_r_xl$ID$.load(path + '_color_r_xl.uni')\n\
        color_g_xl$ID$.load(path + '_color_g_xl.uni')\n\
        color_b_xl$ID$.load(path + '_color_b_xl.uni')\n\
    if using_fire_s$ID$:\n\
        flame_xl$ID$.load(path + '_flame_xl.uni')\n\
        fuel_xl$ID$.load(path + '_fuel_xl.uni')\n\
        react_xl$ID$.load(path + '_react_xl.uni')\n";

const std::string smoke_export_low = "\n\
def save_smoke_data_low_$ID$(path):\n\
    density_s$ID$.save(path + '_density.uni')\n\
    flags_s$ID$.save(path + '_flags.uni')\n\
    vel_s$ID$.save(path + '_vel.uni')\n\
    invel_s$ID$.save(path + '_invel.uni')\n\
    pressure_s$ID$.save(path + '_pressure.uni')\n\
    forces_s$ID$.save(path + '_forces.uni')\n\
    x_force_s$ID$.save(path + '_x_force.uni')\n\
    y_force_s$ID$.save(path + '_y_force.uni')\n\
    z_force_s$ID$.save(path + '_z_force.uni')\n\
    inflow_s$ID$.save(path + '_inflow.uni')\n\
    x_vel_s$ID$.save(path + '_x_vel.uni')\n\
    y_vel_s$ID$.save(path + '_y_vel.uni')\n\
    z_vel_s$ID$.save(path + '_z_vel.uni')\n\
    x_invel_s$ID$.save(path + '_x_invel.uni')\n\
    y_invel_s$ID$.save(path + '_y_invel.uni')\n\
    z_invel_s$ID$.save(path + '_z_invel.uni')\n\
    phiObs_s$ID$.save(path + '_phiObs.uni')\n\
    phiOutIn_s$ID$.save(path + '_phiOutIn.uni')\n\
    if using_colors_s$ID$:\n\
        color_r_s$ID$.save(path + '_color_r.uni')\n\
        color_g_s$ID$.save(path + '_color_g.uni')\n\
        color_b_s$ID$.save(path + '_color_b.uni')\n\
    if using_heat_s$ID$:\n\
        heat_s$ID$.save(path + '_heat.uni')\n\
    if using_fire_s$ID$:\n\
        flame_s$ID$.save(path + '_flame.uni')\n\
        fuel_s$ID$.save(path + '_fuel.uni')\n\
        react_s$ID$.save(path + '_react.uni')\n";

const std::string smoke_export_high = "\n\
def save_smoke_data_high_$ID$(path):\n\
    density_xl$ID$.save(path + '_density_xl.uni')\n\
    flags_xl$ID$.save(path + '_flags_xl.uni')\n\
    \n\
    texture_u_s$ID$.save(path + '_texture_u.uni')\n\
    texture_v_s$ID$.save(path + '_texture_v.uni')\n\
    texture_w_s$ID$.save(path + '_texture_w.uni')\n\
    texture_u2_s$ID$.save(path + '_texture_u2.uni')\n\
    texture_v2_s$ID$.save(path + '_texture_v2.uni')\n\
    texture_w2_s$ID$.save(path + '_texture_w2.uni')\n\
    \n\
    if using_colors_s$ID$:\n\
        color_r_xl$ID$.save(path + '_color_r_xl.uni')\n\
        color_g_xl$ID$.save(path + '_color_g_xl.uni')\n\
        color_b_xl$ID$.save(path + '_color_b_xl.uni')\n\
    if using_fire_s$ID$:\n\
        flame_xl$ID$.save(path + '_flame_xl.uni')\n\
        fuel_xl$ID$.save(path + '_fuel_xl.uni')\n\
        react_xl$ID$.save(path + '_react_xl.uni')\n";

//////////////////////////////////////////////////////////////////////
// DESTRUCTION
//////////////////////////////////////////////////////////////////////

const std::string smoke_delete_colors_low = "\n\
mantaMsg('Deleting colors low')\n\
if 'color_r_s$ID$' in globals() : del color_r_s$ID$\n\
if 'color_g_s$ID$' in globals() : del color_g_s$ID$\n\
if 'color_b_s$ID$' in globals() : del color_b_s$ID$\n";

const std::string smoke_delete_colors_high = "\n\
mantaMsg('Deleting colors high')\n\
if 'color_r_xl$ID$' in globals() : del color_r_xl$ID$\n\
if 'color_g_xl$ID$' in globals() : del color_g_xl$ID$\n\
if 'color_b_xl$ID$' in globals() : del color_b_xl$ID$\n";

const std::string smoke_delete_fire_low = "\n\
mantaMsg('Deleting fire low')\n\
if 'flame_s$ID$' in globals() : del flame_s$ID$\n\
if 'fuel_s$ID$'  in globals() : del fuel_s$ID$\n\
if 'react_s$ID$' in globals() : del react_s$ID$\n";

const std::string smoke_delete_fire_high = "\n\
mantaMsg('Deleting fire high')\n\
if 'flame_xl$ID$' in globals() : del flame_xl$ID$\n\
if 'fuel_xl$ID$'  in globals() : del fuel_xl$ID$\n\
if 'react_xl$ID$' in globals() : del react_xl$ID$\n";

const std::string smoke_delete_heat_low = "\n\
mantaMsg('Deleting heat low')\n\
if 'heat_s$ID$' in globals() : del heat_s$ID$\n";

const std::string smoke_delete_grids_low = "\n\
mantaMsg('Deleting base grids low')\n\
if 'flags_s$ID$'       in globals() : del flags_s$ID$\n\
if 'vel_s$ID$'         in globals() : del vel_s$ID$\n\
if 'x_vel_s$ID$'       in globals() : del x_vel_s$ID$\n\
if 'y_vel_s$ID$'       in globals() : del y_vel_s$ID$\n\
if 'z_vel_s$ID$'       in globals() : del z_vel_s$ID$\n\
if 'invel_s$ID$'       in globals() : del invel_s$ID$\n\
if 'x_invel_s$ID$'     in globals() : del x_invel_s$ID$\n\
if 'y_invel_s$ID$'     in globals() : del y_invel_s$ID$\n\
if 'z_invel_s$ID$'     in globals() : del z_invel_s$ID$\n\
if 'density_s$ID$'     in globals() : del density_s$ID$\n\
if 'pressure_s$ID$'    in globals() : del pressure_s$ID$\n\
if 'phiOutIn_s$ID$'    in globals() : del phiOutIn_s$ID$\n\
if 'phiObs_s$ID$'      in globals() : del phiObs_s$ID$\n\
if 'forces_s$ID$'      in globals() : del forces_s$ID$\n\
if 'x_force_s$ID$'     in globals() : del x_force_s$ID$\n\
if 'y_force_s$ID$'     in globals() : del y_force_s$ID$\n\
if 'z_force_s$ID$'     in globals() : del z_force_s$ID$\n\
if 'inflow_s$ID$'      in globals() : del inflow_s$ID$\n";

const std::string smoke_delete_grids_high = "\n\
mantaMsg('Deleting base grids high')\n\
if 'flags_xl$ID$'      in globals() : del flags_xl$ID$\n\
if 'vel_xl$ID$'        in globals() : del vel_xl$ID$\n\
if 'density_xl$ID$'    in globals() : del density_xl$ID$\n\
if 'phiObs_xl$ID$'     in globals() : del phiObs_xl$ID$\n\
if 'energy_s$ID$'      in globals() : del energy_s$ID$\n\
if 'tempFlag_s$ID$'    in globals() : del tempFlag_s$ID$\n\
if 'uvGrid_s$ID$'      in globals() : del uvGrid_s$ID$\n\
if 'texture_u_s$ID$'   in globals() : del texture_u_s$ID$\n\
if 'texture_v_s$ID$'   in globals() : del texture_v_s$ID$\n\
if 'texture_w_s$ID$'   in globals() : del texture_w_s$ID$\n\
if 'texture_u2_s$ID$'  in globals() : del texture_u2_s$ID$\n\
if 'texture_v2_s$ID$'  in globals() : del texture_v2_s$ID$\n\
if 'texture_w2_s$ID$'  in globals() : del texture_w2_s$ID$\n\
if 'wltnoise_xl$ID$'   in globals() : del wltnoise_xl$ID$\n";

const std::string smoke_delete_variables_low = "\n\
mantaMsg('Deleting variables low')\n\
if 'res_s$ID$'             in globals() : del res_s$ID$\n\
if 'dim_s$ID$'             in globals() : del dim_s$ID$\n\
if 'gs_s$ID$'              in globals() : del gs_s$ID$\n\
if 'doOpen_s$ID$'          in globals() : del doOpen_s$ID$\n\
if 'boundConditions_s$ID$' in globals() : del boundConditions_s$ID$\n\
if 'dt_default_s$ID$'      in globals() : del dt_default_s$ID$\n\
if 'dt_factor_s$ID$'       in globals() : del dt_factor_s$ID$\n\
if 'fps_s$ID$'             in globals() : del fps_s$ID$\n\
if 'dt0_s$ID$'             in globals() : del dt0_s$ID$\n\
if 'vorticity_s$ID$'       in globals() : del vorticity_s$ID$\n\
if 'boundaryWidth_s$ID$'   in globals() : del boundaryWidth_s$ID$\n\
if 'using_colors_s$ID$'    in globals() : del using_colors_s$ID$\n\
if 'using_heat_s$ID$'      in globals() : del using_heat_s$ID$\n\
if 'using_fire_s$ID$'      in globals() : del using_fire_s$ID$\n\
if 'last_frame_s$ID$'      in globals() : del last_frame_s$ID$\n\
if 'maxvel_s$ID$'          in globals() : del maxvel_s$ID$\n";

const std::string smoke_delete_variables_high = "\n\
mantaMsg('Deleting variables high')\n\
if 'wltStrength_s$ID$'  in globals() : del wltStrength_s$ID$\n\
if 'uvs_s$ID$'          in globals() : del uvs_s$ID$\n\
if 'uv_s$ID$'           in globals() : del uv_s$ID$\n\
if 'octaves_s$ID$'      in globals() : del octaves_s$ID$\n";

//////////////////////////////////////////////////////////////////////
// OTHER SETUPS
//////////////////////////////////////////////////////////////////////

const std::string smoke_uv_setup = "\n\
# create the array of uv grids\n\
uv_s$ID$ = []\n\
mantaMsg('Initializing UV Grids')\n\
for i in range(uvs_s$ID$):\n\
    uvGrid_s$ID$ = s$ID$.create(VecGrid)\n\
    uv_s$ID$.append(uvGrid_s$ID$)\n\
    resetUvGrid(uv_s$ID$[i])\n\
\n\
# Need to initialize helper grids for uvw as well\n\
copyVec3ToReal(source=uv_s$ID$[0], targetX=texture_u_s$ID$, targetY=texture_v_s$ID$, targetZ=texture_w_s$ID$)\n\
copyVec3ToReal(source=uv_s$ID$[1], targetX=texture_u2_s$ID$, targetY=texture_v2_s$ID$, targetZ=texture_w2_s$ID$)\n";

const std::string smoke_wavelet_turbulence_noise = "\n\
# wavelet turbulence noise field\n\
mantaMsg('Smoke wavelet noise')\n\
wltnoise_xl$ID$ = xl$ID$.create(NoiseField, loadFromFile=True)\n\
wltnoise_xl$ID$.posScale = vec3(int(1.0*gs_s$ID$.x)) / $NOISE_POSSCALE$\n\
wltnoise_xl$ID$.timeAnim = $NOISE_TIMEANIM$\n";

const std::string smoke_inflow_low = "\n\
def apply_inflow_$ID$():\n\
    mantaMsg('Applying inflow')\n\
    if using_heat_s$ID$ and using_colors_s$ID$ and using_fire_s$ID$:\n\
        applyInflow(density=density_s$ID$, emission=inflow_s$ID$, heat=heat_s$ID$, fuel=fuel_s$ID$, react=react_s$ID$, red=color_r_s$ID$, green=color_g_s$ID$, blue=color_b_s$ID$)\n\
    elif using_heat_s$ID$ and using_colors_s$ID$ and not using_fire_s$ID$:\n\
        applyInflow(density=density_s$ID$, emission=inflow_s$ID$, heat=heat_s$ID$, red=color_r_s$ID$, green=color_g_s$ID$, blue=color_b_s$ID$)\n\
    elif using_heat_s$ID$ and not using_colors_s$ID$ and using_fire_s$ID$:\n\
        applyInflow(density=density_s$ID$, emission=inflow_s$ID$, heat=heat_s$ID$, fuel=fuel_s$ID$, react=react_s$ID$)\n\
    elif using_heat_s$ID$ and not using_colors_s$ID$ and not using_fire_s$ID$:\n\
        applyInflow(density=density_s$ID$, emission=inflow_s$ID$, heat=heat_s$ID$)\n\
    elif not using_heat_s$ID$ and using_colors_s$ID$ and using_fire_s$ID$:\n\
        applyInflow(density=density_s$ID$, emission=inflow_s$ID$, fuel=fuel_s$ID$, react=react_s$ID$, red=color_r_s$ID$, green=color_g_s$ID$, blue=color_b_s$ID$)\n\
    elif not using_heat_s$ID$ and using_colors_s$ID$ and not using_fire_s$ID$:\n\
        applyInflow(density=density_s$ID$, emission=inflow_s$ID$, red=color_r_s$ID$, green=color_g_s$ID$, blue=color_b_s$ID$)\n\
    elif not using_heat_s$ID$ and not using_colors_s$ID$ and using_fire_s$ID$:\n\
        applyInflow(density=density_s$ID$, emission=inflow_s$ID$, fuel=fuel_s$ID$, react=react_s$ID$)\n\
    else:\n\
        applyInflow(density=density_s$ID$, emission=inflow_s$ID$)\n";

const std::string smoke_inflow_high = "\n\
# TODO\n";

//////////////////////////////////////////////////////////////////////
// STANDALONE MODE
//////////////////////////////////////////////////////////////////////

const std::string smoke_standalone_load = "\n\
# import *.uni files\n\
path_prefix = '$MANTA_EXPORT_PATH$'\n\
load_smoke_data_low_$ID$(path_prefix)\n\
load_fluid_data_low_$ID$(path_prefix)\n\
if using_highres_s$ID$:\n\
    load_smoke_data_high_$ID$(path_prefix)\n";
