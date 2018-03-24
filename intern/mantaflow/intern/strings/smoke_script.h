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
vorticity_s$ID$       = $VORTICITY$\n\
absoluteFlow_s$ID$    = True\n";

const std::string smoke_variables_high = "\n\
mantaMsg('Smoke variables high')\n\
wltStrength_s$ID$ = $WLT_STR$\n\
octaves_s$ID$     = 0\n\
uvs_s$ID$         = 2\n\
uv_s$ID$          = [] # list for UV grids\n\
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
mantaMsg('Smoke alloc low')\n\
density_s$ID$    = s$ID$.create(RealGrid)\n\
emissionIn_s$ID$ = s$ID$.create(RealGrid)\n\
shadow_s$ID$     = s$ID$.create(RealGrid)\n\
heat_s$ID$    = 0 # allocated dynamically\n\
flame_s$ID$   = 0\n\
fuel_s$ID$    = 0\n\
react_s$ID$   = 0\n\
color_r_s$ID$ = 0\n\
color_g_s$ID$ = 0\n\
color_b_s$ID$ = 0\n\
\n\
densityIn_s$ID$ = s$ID$.create(RealGrid)\n\
heatIn_s$ID$    = s$ID$.create(RealGrid)\n\
fuelIn_s$ID$    = 0 # allocated dynamically\n\
colorIn_r_s$ID$ = 0\n\
colorIn_g_s$ID$ = 0\n\
colorIn_b_s$ID$ = 0\n";

const std::string smoke_alloc_high = "\n\
mantaMsg('Smoke alloc high')\n\
vel_xl$ID$       = xl$ID$.create(MACGrid)\n\
density_xl$ID$   = xl$ID$.create(RealGrid)\n\
phiOut_xl$ID$    = xl$ID$.create(LevelsetGrid)\n\
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
color_r_s$ID$   = s$ID$.create(RealGrid)\n\
color_g_s$ID$   = s$ID$.create(RealGrid)\n\
color_b_s$ID$   = s$ID$.create(RealGrid)\n\
colorIn_r_s$ID$ = s$ID$.create(RealGrid)\n\
colorIn_g_s$ID$ = s$ID$.create(RealGrid)\n\
colorIn_b_s$ID$ = s$ID$.create(RealGrid)\n";

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
heat_s$ID$   = s$ID$.create(RealGrid)\n\
heatIn_s$ID$ = s$ID$.create(RealGrid)\n";

const std::string smoke_alloc_fire_low = "\n\
mantaMsg('Allocating fire low')\n\
flame_s$ID$  = s$ID$.create(RealGrid)\n\
fuel_s$ID$   = s$ID$.create(RealGrid)\n\
fuelIn_s$ID$ = s$ID$.create(RealGrid)\n\
react_s$ID$  = s$ID$.create(RealGrid)\n";

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
        x_obvel_s$ID$.multConst(Real(gs_s$ID$.x))\n\
        y_obvel_s$ID$.multConst(Real(gs_s$ID$.y))\n\
        z_obvel_s$ID$.multConst(Real(gs_s$ID$.z))\n\
        copyRealToVec3(sourceX=x_obvel_s$ID$, sourceY=y_obvel_s$ID$, sourceZ=z_obvel_s$ID$, target=obvelC_s$ID$)\n\
    \n\
    # translate guiding vels (world space) to grid space\n\
    if using_guiding_s$ID$:\n\
        x_guidevel_s$ID$.multConst(Real(gs_s$ID$.x))\n\
        y_guidevel_s$ID$.multConst(Real(gs_s$ID$.y))\n\
        z_guidevel_s$ID$.multConst(Real(gs_s$ID$.z))\n\
        copyRealToVec3(sourceX=x_guidevel_s$ID$, sourceY=y_guidevel_s$ID$, sourceZ=z_guidevel_s$ID$, target=guidevelC_s$ID$)\n\
    \n\
    # translate invels (world space) to grid space\n\
    if using_invel_s$ID$:\n\
        x_invel_s$ID$.multConst(Real(gs_s$ID$.x))\n\
        y_invel_s$ID$.multConst(Real(gs_s$ID$.y))\n\
        z_invel_s$ID$.multConst(Real(gs_s$ID$.z))\n\
        copyRealToVec3(sourceX=x_invel_s$ID$, sourceY=y_invel_s$ID$, sourceZ=z_invel_s$ID$, target=invel_s$ID$)\n\
    \n\
    #x_vel_s$ID$.multConst(Real(gs_s$ID$.x))\n\
    #y_vel_s$ID$.multConst(Real(gs_s$ID$.y))\n\
    #z_vel_s$ID$.multConst(Real(gs_s$ID$.z))\n\
    #copyRealToVec3(sourceX=x_vel_s$ID$, sourceY=y_vel_s$ID$, sourceZ=z_vel_s$ID$, target=vel_s$ID$)\n\
    copyRealToVec3(sourceX=x_force_s$ID$, sourceY=y_force_s$ID$, sourceZ=z_force_s$ID$, target=forces_s$ID$)\n\
    \n\
    # If obstacle has velocity, i.e. is moving switch to dynamic preconditioner\n\
    if using_obstacle_s$ID$ and obvelC_s$ID$.getMax() > 0:\n\
        mantaMsg('Using dynamic preconditioner')\n\
        preconditioner_s$ID$ = PcMGDynamic\n\
    else:\n\
        mantaMsg('Using static preconditioner')\n\
        preconditioner_s$ID$ = PcMGStatic\n";

const std::string smoke_pre_step_high = "\n\
def smoke_pre_step_high_$ID$():\n\
    mantaMsg('Smoke pre step high')\n\
    global uv_s$ID$\n\
    if len(uv_s$ID$) != 0: # list of uvs already initialized?\n\
        copyRealToVec3(sourceX=texture_u_s$ID$, sourceY=texture_v_s$ID$, sourceZ=texture_w_s$ID$, target=uv_s$ID$[0])\n\
        copyRealToVec3(sourceX=texture_u2_s$ID$, sourceY=texture_v2_s$ID$, sourceZ=texture_w2_s$ID$, target=uv_s$ID$[1])\n\
    else:\n\
        mantaMsg('Initializing UV Grids')\n\
        for i in range(uvs_s$ID$):\n\
            uvGrid_s$ID$ = s$ID$.create(VecGrid)\n\
            uv_s$ID$.append(uvGrid_s$ID$)\n\
            resetUvGrid(uv_s$ID$[i])\n";

const std::string smoke_post_step_low = "\n\
def smoke_post_step_low_$ID$():\n\
    mantaMsg('Smoke post step low')\n\
    forces_s$ID$.clear()\n\
    if using_guiding_s$ID$:\n\
        weightGuide_s$ID$.clear()\n\
    if using_invel_s$ID$:\n\
        invel_s$ID$.clear()\n\
    \n\
    phiObs_s$ID$.setConst(9999)\n\
    phiOutIn_s$ID$.setConst(9999)\n\
    \n\
    #copyVec3ToReal(source=vel_s$ID$, targetX=x_vel_s$ID$, targetY=y_vel_s$ID$, targetZ=z_vel_s$ID$)\n\
    #x_vel_s$ID$.multConst( 1.0/Real(gs_s$ID$.x) )\n\
    #y_vel_s$ID$.multConst( 1.0/Real(gs_s$ID$.y) )\n\
    #z_vel_s$ID$.multConst( 1.0/Real(gs_s$ID$.z) )\n";

const std::string smoke_post_step_high = "\n\
def smoke_post_step_high_$ID$():\n\
    mantaMsg('Smoke post step high')\n\
    inflow_s$ID$.setConst(0)\n\
    copyVec3ToReal(source=uv_s$ID$[0], targetX=texture_u_s$ID$, targetY=texture_v_s$ID$, targetZ=texture_w_s$ID$)\n\
    copyVec3ToReal(source=uv_s$ID$[1], targetX=texture_u2_s$ID$, targetY=texture_v2_s$ID$, targetZ=texture_w2_s$ID$)\n";

//////////////////////////////////////////////////////////////////////
// STEP FUNCTIONS
//////////////////////////////////////////////////////////////////////

const std::string smoke_adaptive_step_low = "\n\
def fluid_step_low_$ID$(framenr):\n\
    s$ID$.frame = framenr\n\
    s$ID$.timeTotal = s$ID$.frame * dt0_s$ID$\n\
    last_frame_s$ID$ = s$ID$.frame\n\
    \n\
    smoke_pre_step_low_$ID$()\n\
    \n\
    while s$ID$.frame == last_frame_s$ID$:\n\
        \n\
        applyGasInflow(emission=emissionIn_s$ID$, density=density_s$ID$, densityIn=densityIn_s$ID$, heat=heat_s$ID$, heatIn=heatIn_s$ID$, fuel=fuel_s$ID$, fuelIn=fuelIn_s$ID$, react=react_s$ID$, red=color_r_s$ID$, redIn=colorIn_r_s$ID$, green=color_g_s$ID$, greenIn=colorIn_g_s$ID$, blue=color_b_s$ID$, blueIn=colorIn_b_s$ID$, flowType=flowType_s$ID$, absoluteFlow=absoluteFlow_s$ID$)\n\
        \n\
        mantaMsg('s.frame is ' + str(s$ID$.frame))\n\
        if using_obstacle_s$ID$: # TODO (sebbas): allow outflow objects when no obstacle set\n\
            phiObs_s$ID$.join(phiObsIn_s$ID$)\n\
        setObstacleFlags(flags=flags_s$ID$, phiObs=phiObs_s$ID$, phiOut=phiOutIn_s$ID$)\n\
        flags_s$ID$.fillGrid()\n\
        \n\
        fluid_adapt_time_step_low()\n\
        mantaMsg('Low step / s$ID$.frame: ' + str(s$ID$.frame))\n\
        if using_fire_s$ID$:\n\
            process_burn_low_$ID$()\n\
        step_low_$ID$()\n\
        if using_fire_s$ID$:\n\
            update_flame_low_$ID$()\n\
        \n\
        s$ID$.step()\n\
    \n\
    smoke_post_step_low_$ID$()\n";

const std::string smoke_adaptive_step_high = "\n\
def fluid_step_high_$ID$(framenr):\n\
    xl$ID$.frame = framenr\n\
    xl$ID$.timeTotal = xl$ID$.frame * dt0_s$ID$\n\
    last_frame_s$ID$ = xl$ID$.frame\n\
    \n\
    smoke_pre_step_high_$ID$()\n\
    \n\
    while xl$ID$.frame == last_frame_s$ID$:\n\
        \n\
        mantaMsg('xl.frame is ' + str(xl$ID$.frame))\n\
        setObstacleFlags(flags=flags_xl$ID$, phiObs=phiObs_xl$ID$, phiOut=phiOut_xl$ID$)\n\
        flags_xl$ID$.fillGrid()\n\
        \n\
        fluid_adapt_time_step_high()\n\
        mantaMsg('High step / xl$ID$.frame: ' + str(xl$ID$.frame))\n\
        if using_fire_s$ID$:\n\
            process_burn_high_$ID$()\n\
        step_high_$ID$()\n\
        if using_fire_s$ID$:\n\
            update_flame_high_$ID$()\n\
        \n\
        xl$ID$.step()\n\
    \n\
    smoke_post_step_high_$ID$()\n";

const std::string smoke_step_low = "\n\
def step_low_$ID$():\n\
    mantaMsg('Smoke step low')\n\
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
    # Create interpolated version of original phi grids for later use in (optional) high-res step\n\
    if using_obstacle_s$ID$ and using_highres_s$ID$:\n\
        interpolateGrid(target=phiOut_xl$ID$, source=phiOutIn_s$ID$)\n\
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
        extrapolateVec3Simple(vel=guidevelC_s$ID$, phi=phiGuideIn_s$ID$, distance=4, inside=False)\n\
        resampleVec3ToMac(source=guidevelC_s$ID$, target=guidevel_s$ID$)\n\
    \n\
    # add initial velocity\n\
    if using_invel_s$ID$:\n\
        setInitialVelocity(flags=flags_s$ID$, vel=vel_s$ID$, invel=invel_s$ID$)\n\
    \n\
    mantaMsg('Walls')\n\
    setWallBcs(flags=flags_s$ID$, vel=vel_s$ID$, obvel=obvel_s$ID$ if using_obstacle_s$ID$ else 0)\n\
    \n\
    if using_guiding_s$ID$:\n\
        mantaMsg('Guiding and pressure')\n\
        weightGuide_s$ID$.addConst(alpha_s$ID$)\n\
        PD_fluid_guiding(vel=vel_s$ID$, velT=guidevel_s$ID$, flags=flags_s$ID$, weight=weightGuide_s$ID$, blurRadius=beta_s$ID$, pressure=pressure_s$ID$, tau=tau_s$ID$, sigma=sigma_s$ID$, theta=theta_s$ID$, preconditioner=preconditioner_s$ID$, zeroPressureFixing=not doOpen_s$ID$)\n\
    else:\n\
        mantaMsg('Pressure')\n\
        solvePressure(flags=flags_s$ID$, vel=vel_s$ID$, pressure=pressure_s$ID$, preconditioner=preconditioner_s$ID$, zeroPressureFixing=not doOpen_s$ID$) # closed domains require pressure fixing\n\
\n\
def process_burn_low_$ID$():\n\
    mantaMsg('Process burn low')\n\
    processBurn(fuel=fuel_s$ID$, density=density_s$ID$, react=react_s$ID$, red=color_r_s$ID$ if using_colors_s$ID$ else 0, green=color_g_s$ID$ if using_colors_s$ID$ else 0, blue=color_b_s$ID$ if using_colors_s$ID$ else 0, heat=heat_s$ID$ if using_heat_s$ID$ else 0, burningRate=$BURNING_RATE$, flameSmoke=$FLAME_SMOKE$, ignitionTemp=$IGNITION_TEMP$, maxTemp=$MAX_TEMP$, flameSmokeColor=vec3($FLAME_SMOKE_COLOR_X$,$FLAME_SMOKE_COLOR_Y$,$FLAME_SMOKE_COLOR_Z$))\n\
\n\
def update_flame_low_$ID$():\n\
    mantaMsg('Update flame low')\n\
    updateFlame(react=react_s$ID$, flame=flame_s$ID$)\n";

const std::string smoke_step_high = "\n\
def step_high_$ID$():\n\
    mantaMsg('Smoke step high')\n\
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
    processBurn(fuel=fuel_xl$ID$, density=density_xl$ID$, react=react_xl$ID$, red=color_r_xl$ID$ if using_colors_s$ID$ else 0, green=color_g_xl$ID$ if using_colors_s$ID$ else 0, blue=color_b_xl$ID$ if using_colors_s$ID$ else 0, burningRate=$BURNING_RATE$, flameSmoke=$FLAME_SMOKE$, ignitionTemp=$IGNITION_TEMP$, maxTemp=$MAX_TEMP$, dt=dt0_s$ID$, flameSmokeColor=vec3($FLAME_SMOKE_COLOR_X$,$FLAME_SMOKE_COLOR_Y$,$FLAME_SMOKE_COLOR_Z$))\n\
\n\
def update_flame_high_$ID$():\n\
    mantaMsg('Update flame high')\n\
    updateFlame(react=react_xl$ID$, flame=flame_xl$ID$)\n";

//////////////////////////////////////////////////////////////////////
// IMPORT
//////////////////////////////////////////////////////////////////////

const std::string smoke_load_geometry_low = "\n\
def smoke_load_geometry_low_$ID$(path, framenr):\n\
    mantaMsg('Smoke load geometry')\n\
    framenr = fluid_cache_get_framenr_formatted_$ID$(framenr)\n\
    densityIn_s$ID$.load(os.path.join(path, 'densityIn_' + framenr + '.uni'))\n\
    if using_colors_s$ID$:\n\
        colorIn_r_s$ID$.load(os.path.join(path, 'colorIn_r_' + framenr + '.uni'))\n\
        colorIn_g_s$ID$.load(os.path.join(path, 'colorIn_g_' + framenr + '.uni'))\n\
        colorIn_b_s$ID$.load(os.path.join(path, 'colorIn_b_' + framenr + '.uni'))\n\
    if using_heat_s$ID$:\n\
        heatIn_s$ID$.load(os.path.join(path, 'heatIn_' + framenr + '.uni'))\n\
    if using_fire_s$ID$:\n\
        fuelIn_s$ID$.load(os.path.join(path, 'fuelIn_' + framenr + '.uni'))\n\
    emissionIn_s$ID$.load(os.path.join(path, 'emissionIn_' + framenr + '.uni'))\n";

const std::string smoke_load_data_low = "\n\
def smoke_load_data_low_$ID$(path, framenr):\n\
    mantaMsg('Smoke load data low')\n\
    framenr = fluid_cache_get_framenr_formatted_$ID$(framenr)\n\
    density_s$ID$.load(os.path.join(path, 'density_' + framenr + '.uni'))\n\
    shadow_s$ID$.load(os.path.join(path, 'shadow_' + framenr + '.uni'))\n\
    if using_colors_s$ID$:\n\
        color_r_s$ID$.load(os.path.join(path, 'color_r_' + framenr + '.uni'))\n\
        color_g_s$ID$.load(os.path.join(path, 'color_g_' + framenr + '.uni'))\n\
        color_b_s$ID$.load(os.path.join(path, 'color_b_' + framenr + '.uni'))\n\
    if using_heat_s$ID$:\n\
        heat_s$ID$.load(os.path.join(path, 'heat_' + framenr + '.uni'))\n\
    if using_fire_s$ID$:\n\
        flame_s$ID$.load(os.path.join(path, 'flame_' + framenr + '.uni'))\n\
        fuel_s$ID$.load(os.path.join(path, 'fuel_' + framenr + '.uni'))\n\
        react_s$ID$.load(os.path.join(path, 'react_' + framenr + '.uni'))\n";

const std::string smoke_load_data_high = "\n\
def smoke_load_data_high_$ID$(path, framenr):\n\
    mantaMsg('Smoke load data high')\n\
    framenr = fluid_cache_get_framenr_formatted_$ID$(framenr)\n\
    density_xl$ID$.load(os.path.join(path, 'density_xl_' + framenr + '.uni'))\n\
    texture_u_s$ID$.load(os.path.join(path, 'texture_u_' + framenr + '.uni'))\n\
    texture_v_s$ID$.load(os.path.join(path, 'texture_v_' + framenr + '.uni'))\n\
    texture_w_s$ID$.load(os.path.join(path, 'texture_w_' + framenr + '.uni'))\n\
    texture_u2_s$ID$.load(os.path.join(path, 'texture_u2_' + framenr + '.uni'))\n\
    texture_v2_s$ID$.load(os.path.join(path, 'texture_v2_' + framenr + '.uni'))\n\
    texture_w2_s$ID$.load(os.path.join(path, 'texture_w2_' + framenr + '.uni'))\n\
    if using_colors_s$ID$:\n\
        color_r_xl$ID$.load(os.path.join(path, 'color_r_xl_' + framenr + '.uni'))\n\
        color_g_xl$ID$.load(os.path.join(path, 'color_g_xl_' + framenr + '.uni'))\n\
        color_b_xl$ID$.load(os.path.join(path, 'color_b_xl_' + framenr + '.uni'))\n\
    if using_fire_s$ID$:\n\
        flame_xl$ID$.load(os.path.join(path, 'flame_xl_' + framenr + '.uni'))\n\
        fuel_xl$ID$.load(os.path.join(path, 'fuel_xl_' + framenr + '.uni'))\n\
        react_xl$ID$.load(os.path.join(path, 'react_xl_' + framenr + '.uni'))\n";

//////////////////////////////////////////////////////////////////////
// EXPORT
//////////////////////////////////////////////////////////////////////

const std::string smoke_save_geometry_low = "\n\
def smoke_save_geometry_low_$ID$(path, framenr):\n\
    mantaMsg('Smoke save geometry')\n\
    framenr = fluid_cache_get_framenr_formatted_$ID$(framenr)\n\
    densityIn_s$ID$.save(os.path.join(path, 'densityIn_' + framenr + '.uni'))\n\
    if using_colors_s$ID$:\n\
        colorIn_r_s$ID$.save(os.path.join(path, 'colorIn_r_' + framenr + '.uni'))\n\
        colorIn_g_s$ID$.save(os.path.join(path, 'colorIn_g_' + framenr + '.uni'))\n\
        colorIn_b_s$ID$.save(os.path.join(path, 'colorIn_b_' + framenr + '.uni'))\n\
    if using_heat_s$ID$:\n\
        heatIn_s$ID$.save(os.path.join(path, 'heatIn_' + framenr + '.uni'))\n\
    if using_fire_s$ID$:\n\
        fuelIn_s$ID$.save(os.path.join(path, 'fuelIn_' + framenr + '.uni'))\n\
    emissionIn_s$ID$.save(os.path.join(path, 'emissionIn_' + framenr + '.uni'))\n";

const std::string smoke_save_data_low = "\n\
def smoke_save_data_low_$ID$(path, framenr):\n\
    mantaMsg('Smoke save data low')\n\
    framenr = fluid_cache_get_framenr_formatted_$ID$(framenr)\n\
    density_s$ID$.save(os.path.join(path, 'density_' + framenr + '.uni'))\n\
    if using_colors_s$ID$:\n\
        color_r_s$ID$.save(os.path.join(path, 'color_r_' + framenr + '.uni'))\n\
        color_g_s$ID$.save(os.path.join(path, 'color_g_' + framenr + '.uni'))\n\
        color_b_s$ID$.save(os.path.join(path, 'color_b_' + framenr + '.uni'))\n\
    if using_heat_s$ID$:\n\
        heat_s$ID$.save(os.path.join(path, 'heat_' + framenr + '.uni'))\n\
    if using_fire_s$ID$:\n\
        flame_s$ID$.save(os.path.join(path, 'flame_' + framenr + '.uni'))\n\
        fuel_s$ID$.save(os.path.join(path, 'fuel_' + framenr + '.uni'))\n\
        react_s$ID$.save(os.path.join(path, 'react_' + framenr + '.uni'))\n\
\n\
def smoke_save_shadow_$ID$(path, framenr):\n\
    mantaMsg('Smoke save shadow')\n\
    framenr = fluid_cache_get_framenr_formatted_$ID$(framenr)\n\
    shadow_s$ID$.save(os.path.join(path, 'shadow_' + framenr + '.uni'))\n";

const std::string smoke_save_data_high = "\n\
def smoke_save_data_high_$ID$(path, framenr):\n\
    mantaMsg('Smoke save data high')\n\
    framenr = fluid_cache_get_framenr_formatted_$ID$(framenr)\n\
    density_xl$ID$.save(os.path.join(path, 'density_xl_' + framenr + '.uni'))\n\
    texture_u_s$ID$.save(os.path.join(path, 'texture_u_' + framenr + '.uni'))\n\
    texture_v_s$ID$.save(os.path.join(path, 'texture_v_' + framenr + '.uni'))\n\
    texture_w_s$ID$.save(os.path.join(path, 'texture_w_' + framenr + '.uni'))\n\
    texture_u2_s$ID$.save(os.path.join(path, 'texture_u2_' + framenr + '.uni'))\n\
    texture_v2_s$ID$.save(os.path.join(path, 'texture_v2_' + framenr + '.uni'))\n\
    texture_w2_s$ID$.save(os.path.join(path, 'texture_w2_' + framenr + '.uni')\n\
    if using_colors_s$ID$:\n\
        color_r_xl$ID$.save(os.path.join(path, 'color_r_xl_' + framenr + '.uni'))\n\
        color_g_xl$ID$.save(os.path.join(path, 'color_g_xl_' + framenr + '.uni'))\n\
        color_b_xl$ID$.save(os.path.join(path, 'color_b_xl_' + framenr + '.uni'))\n\
    if using_fire_s$ID$:\n\
        flame_xl$ID$.save(os.path.join(path, 'flame_xl_' + framenr + '.uni'))\n\
        fuel_xl$ID$.save(os.path.join(path, 'fuel_xl_' + framenr + '.uni'))\n\
        react_xl$ID$.save(os.path.join(path, 'react_xl_' + framenr + '.uni'))\n";

//////////////////////////////////////////////////////////////////////
// OTHER SETUPS
//////////////////////////////////////////////////////////////////////

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
path_prefix_$ID$ = '$MANTA_EXPORT_PATH$'\n\
load_smoke_data_low_$ID$(path_prefix_$ID$)\n\
if using_highres_s$ID$:\n\
    load_smoke_data_high_$ID$(path_prefix_$ID$)\n";
