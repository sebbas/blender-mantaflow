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

/** \file mantaflow/intern/strings/liquid.h
 *  \ingroup mantaflow
 */

#include <string>

//////////////////////////////////////////////////////////////////////
// VARIABLES
//////////////////////////////////////////////////////////////////////

const std::string liquid_variables_low = "\n\
mantaMsg('Liquid variables low')\n\
narrowBandWidth_s$ID$         = 3\n\
combineBandWidth_s$ID$        = narrowBandWidth_s$ID$ - 1\n\
adjustedNarrowBandWidth_s$ID$ = $PARTICLE_BAND_WIDTH$ # only used in adjustNumber to control band width\n\
\n\
particleNumber_s$ID$ = $PARTICLE_NUMBER$\n\
minParticles_s$ID$   = $PARTICLE_MINIMUM$\n\
maxParticles_s$ID$   = $PARTICLE_MAXIMUM$\n\
radiusFactor_s$ID$   = $PARTICLE_RADIUS$\n\
randomness_s$ID$     = $PARTICLE_RANDOMNESS$\n\
maxVel_s$ID$         = 1 # just declared here, do not set\n";

const std::string liquid_variables_high = "\n\
mantaMsg('Liquid variables high')\n";

//////////////////////////////////////////////////////////////////////
// GRIDS & MESH & PARTICLESYSTEM
//////////////////////////////////////////////////////////////////////

const std::string liquid_alloc_low = "\n\
mantaMsg('Liquid alloc low')\n\
flags_s$ID$      = s$ID$.create(FlagGrid)\n\
phiParts_s$ID$   = s$ID$.create(LevelsetGrid)\n\
phi_s$ID$        = s$ID$.create(LevelsetGrid)\n\
phiIn_s$ID$      = s$ID$.create(LevelsetGrid)\n\
phiOut_s$ID$     = s$ID$.create(LevelsetGrid)\n\
phiOutIn_s$ID$   = s$ID$.create(LevelsetGrid)\n\
pressure_s$ID$   = s$ID$.create(RealGrid)\n\
\n\
phiObs_s$ID$     = s$ID$.create(LevelsetGrid)\n\
fractions_s$ID$  = 0 # s$ID$.create(MACGrid) # TODO (sebbas): disabling fractions for now - not fracwallbcs not supporting obvels yet\n\
\n\
vel_s$ID$        = s$ID$.create(MACGrid)\n\
x_vel_s$ID$      = s$ID$.create(RealGrid)\n\
y_vel_s$ID$      = s$ID$.create(RealGrid)\n\
z_vel_s$ID$      = s$ID$.create(RealGrid)\n\
\n\
velOld_s$ID$     = s$ID$.create(MACGrid)\n\
velParts_s$ID$   = s$ID$.create(MACGrid)\n\
mapWeights_s$ID$ = s$ID$.create(MACGrid)\n\
\n\
pp_s$ID$         = s$ID$.create(BasicParticleSystem)\n\
pVel_pp$ID$      = pp_s$ID$.create(PdataVec3)\n\
mesh_s$ID$       = s$ID$.create(Mesh)\n\
\n\
# Acceleration data for particle nbs\n\
pindex_s$ID$     = s$ID$.create(ParticleIndexSystem)\n\
gpi_s$ID$        = s$ID$.create(IntGrid)\n\
\n\
forces_s$ID$     = s$ID$.create(MACGrid)\n\
x_force_s$ID$    = s$ID$.create(RealGrid)\n\
y_force_s$ID$    = s$ID$.create(RealGrid)\n\
z_force_s$ID$    = s$ID$.create(RealGrid)\n";

const std::string liquid_alloc_high = "\n\
mantaMsg('Liquid alloc high')\n\
flags_xl$ID$    = xl$ID$.create(FlagGrid)\n\
phiParts_xl$ID$ = xl$ID$.create(LevelsetGrid)\n\
phi_xl$ID$      = xl$ID$.create(LevelsetGrid)\n\
pp_xl$ID$       = xl$ID$.create(BasicParticleSystem)\n\
mesh_xl$ID$     = xl$ID$.create(Mesh)\n\
\n\
# Acceleration data for particle nbs\n\
pindex_xl$ID$  = xl$ID$.create(ParticleIndexSystem)\n\
gpi_xl$ID$     = xl$ID$.create(IntGrid)\n";

const std::string liquid_init_phi = "\n\
phi_s$ID$.initFromFlags(flags_s$ID$)\n\
phiIn_s$ID$.initFromFlags(flags_s$ID$)\n";

//////////////////////////////////////////////////////////////////////
// PRE / POST STEP
//////////////////////////////////////////////////////////////////////

const std::string liquid_pre_step_low = "\n\
def liquid_pre_step_low_$ID$():\n\
    # translate obvels (world space) to grid space\n\
    if using_obstacle_s$ID$:\n\
        x_obvel_s$ID$.multConst(gs_s$ID$.x)\n\
        y_obvel_s$ID$.multConst(gs_s$ID$.y)\n\
        z_obvel_s$ID$.multConst(gs_s$ID$.z)\n\
        copyRealToVec3(sourceX=x_obvel_s$ID$, sourceY=y_obvel_s$ID$, sourceZ=z_obvel_s$ID$, target=obvelC_s$ID$)\n\
    \n\
    # translate guiding vels (world space) to grid space\n\
    if using_guiding_s$ID$:\n\
        x_guidevel_s$ID$.multConst(gs_s$ID$.x)\n\
        y_guidevel_s$ID$.multConst(gs_s$ID$.y)\n\
        z_guidevel_s$ID$.multConst(gs_s$ID$.z)\n\
        copyRealToVec3(sourceX=x_guidevel_s$ID$, sourceY=y_guidevel_s$ID$, sourceZ=z_guidevel_s$ID$, target=guidevelC_s$ID$)\n\
    \n\
    # translate invels (world space) to grid space\n\
    if using_invel_s$ID$:\n\
        x_invel_s$ID$.multConst(gs_s$ID$.x)\n\
        y_invel_s$ID$.multConst(gs_s$ID$.y)\n\
        z_invel_s$ID$.multConst(gs_s$ID$.z)\n\
        copyRealToVec3(sourceX=x_invel_s$ID$, sourceY=y_invel_s$ID$, sourceZ=z_invel_s$ID$, target=invel_s$ID$)\n\
    \n\
    copyRealToVec3(sourceX=x_vel_s$ID$, sourceY=y_vel_s$ID$, sourceZ=z_vel_s$ID$, target=vel_s$ID$)\n\
    copyRealToVec3(sourceX=x_force_s$ID$, sourceY=y_force_s$ID$, sourceZ=z_force_s$ID$, target=forces_s$ID$)\n";

const std::string liquid_post_step_low = "\n\
def liquid_post_step_low_$ID$():\n\
    forces_s$ID$.clear()\n\
    if using_guiding_s$ID$:\n\
        weightGuide_s$ID$.clear()\n\
    if using_invel_s$ID$:\n\
        invel_s$ID$.clear()\n\
    \n\
    phiOut_s$ID$.setConst(9999)\n\
    copyVec3ToReal(source=vel_s$ID$, targetX=x_vel_s$ID$, targetY=y_vel_s$ID$, targetZ=z_vel_s$ID$)\n";

//////////////////////////////////////////////////////////////////////
// STEP FUNCTIONS
//////////////////////////////////////////////////////////////////////

const std::string liquid_adaptive_step = "\n\
def manta_step_$ID$(framenr):\n\
    mantaMsg('Manta step, frame ' + str(framenr))\n\
    \n\
    s$ID$.frame = framenr\n\
    s$ID$.timeTotal = s$ID$.frame * dt0_s$ID$\n\
    last_frame_s$ID$ = s$ID$.frame\n\
    \n\
    liquid_pre_step_low_$ID$()\n\
    \n\
    while s$ID$.frame == last_frame_s$ID$:\n\
        \n\
        mantaMsg('s.frame is ' + str(s$ID$.frame))\n\
        flags_s$ID$.initDomain(boundaryWidth=boundaryWidth_s$ID$, phiWalls=phiObs_s$ID$, outflow=boundConditions_s$ID$)\n\
        \n\
        # extrapolate before joining inflow levelsets\n\
        extrapolateLsSimple(phi=phiIn_s$ID$, distance=narrowBandWidth_s$ID$+2, inside=True)\n\
        extrapolateLsSimple(phi=phiIn_s$ID$, distance=3)\n\
        if using_obstacle_s$ID$:\n\
            extrapolateLsSimple(phi=phiObsIn_s$ID$, distance=9, inside=True)\n\
            extrapolateLsSimple(phi=phiObsIn_s$ID$, distance=9)\n\
        \n\
        if using_obstacle_s$ID$:\n\
            phiObs_s$ID$.join(phiObsIn_s$ID$)\n\
        phi_s$ID$.join(phiIn_s$ID$)\n\
        if using_obstacle_s$ID$:\n\
            phi_s$ID$.subtract(phiObsIn_s$ID$)\n\
        \n\
        phiOut_s$ID$.join(phiOutIn_s$ID$)\n\
        \n\
        #updateFractions(flags=flags_s$ID$, phiObs=phiObs_s$ID$, fractions=fractions_s$ID$, boundaryWidth=boundaryWidth_s$ID$) # TODO (sebbas): uncomment for fraction support\n\
        setObstacleFlags(flags=flags_s$ID$, phiObs=phiObs_s$ID$, phiOut=phiOut_s$ID$, fractions=fractions_s$ID$)\n\
        \n\
        # add initial velocity: set invel as source grid to ensure const vels in inflow region, sampling makes use of this\n\
        mapWeights_s$ID$.clear() # mis-use mapWeights\n\
        if using_invel_s$ID$:\n\
            resampleVec3ToMac(source=invel_s$ID$, target=mapWeights_s$ID$)\n\
        pVel_pp$ID$.setSource(mapWeights_s$ID$, isMAC=True)\n\
        \n\
        sampleLevelsetWithParticles(phi=phiIn_s$ID$, flags=flags_s$ID$, parts=pp_s$ID$, discretization=particleNumber_s$ID$, randomness=randomness_s$ID$, refillEmpty=True)\n\
        flags_s$ID$.updateFromLevelset(phi_s$ID$, phiObs_s$ID$)\n\
        mapWeights_s$ID$.clear() # clean up, mapweights grid used later again\n\
        \n\
        maxVel_s$ID$ = vel_s$ID$.getMaxValue()\n\
        if using_adaptTime_s$ID$:\n\
            mantaMsg('Adapt timestep')\n\
            s$ID$.adaptTimestep(maxVel_s$ID$)\n\
        \n\
        mantaMsg('Low step / s$ID$.frame: ' + str(s$ID$.frame))\n\
        liquid_step_$ID$()\n\
        \n\
        if using_highres_s$ID$:\n\
            xl$ID$.timestep = s$ID$.timestep\n\
            mantaMsg('High step / s$ID$.frame: ' + str(s$ID$.frame))\n\
            liquid_step_high_$ID$()\n\
        s$ID$.step()\n\
    \n\
    liquid_post_step_low_$ID$()\n";

const std::string liquid_step_low = "\n\
def liquid_step_$ID$():\n\
    mantaMsg('Liquid step low')\n\
    \n\
    if using_drops_s$ID$:\n\
        mantaMsg('Sampling drop particles')\n\
        sampleSndParts(type=PtypeDroplet, constraint=$SNDPARTICLE_VEL_THRESH$, phi=phi_s$ID$, flags=flags_s$ID$, vel=vel_s$ID$, parts=ppSnd_s$ID$)\n\
    \n\
    if using_floats_s$ID$:\n\
        mantaMsg('Sampling float particles')\n\
        sampleSndParts(type=PtypeFloater, constraint=$SNDPARTICLE_FLOAT_AMOUNT$, phi=phiIn_s$ID$, flags=flags_s$ID$, vel=vel_s$ID$, parts=ppSnd_s$ID$)\n\
    \n\
    if using_tracers_s$ID$:\n\
        mantaMsg('Sampling tracer particles')\n\
        sampleSndParts(type=PtypeTracer, constraint=$SNDPARTICLE_TRACER_AMOUNT$, phi=phiIn_s$ID$, flags=flags_s$ID$, vel=vel_s$ID$, parts=ppSnd_s$ID$)\n\
    \n\
    if using_sndparts_s$ID$:\n\
        mantaMsg('Updating snd particle data (velocity, life count)')\n\
        updateSndParts(phi=phi_s$ID$, flags=flags_s$ID$, vel=vel_s$ID$, gravity=gravity_s$ID$, parts=ppSnd_s$ID$, partVel=pVelSnd_pp$ID$, partLife=pLifeSnd_pp$ID$)\n\
        \n\
        mantaMsg('Adjusting snd particles')\n\
        pushOutofObs(parts=ppSnd_s$ID$, flags=flags_s$ID$, phiObs=phiObs_s$ID$, shift=1.0)\n\
        adjustSndParts(parts=ppSnd_s$ID$, flags=flags_s$ID$, phi=phi_s$ID$, partVel=pVelSnd_pp$ID$)\n\
    \n\
    mantaMsg('Advecting particles')\n\
    pp_s$ID$.advectInGrid(flags=flags_s$ID$, vel=vel_s$ID$, integrationMode=IntRK4, deleteInObstacle=using_obstacle_s$ID$, stopInObstacle=False)\n\
    \n\
    mantaMsg('Pushing particles out of obstacles')\n\
    pushOutofObs(parts=pp_s$ID$, flags=flags_s$ID$, phiObs=phiObs_s$ID$)\n\
    \n\
    mantaMsg('Advecting phi')\n\
    advectSemiLagrange(flags=flags_s$ID$, vel=vel_s$ID$, grid=phi_s$ID$, order=1) # first order is usually enough\n\
    mantaMsg('Advecting velocity')\n\
    advectSemiLagrange(flags=flags_s$ID$, vel=vel_s$ID$, grid=vel_s$ID$, order=2, openBounds=doOpen_s$ID$, boundaryWidth=boundaryWidth_s$ID$)\n\
    \n\
    # create level set of particles\n\
    gridParticleIndex(parts=pp_s$ID$, flags=flags_s$ID$, indexSys=pindex_s$ID$, index=gpi_s$ID$)\n\
    unionParticleLevelset(pp_s$ID$, pindex_s$ID$, flags_s$ID$, gpi_s$ID$, phiParts_s$ID$)\n\
    \n\
    # combine level set of particles with grid level set\n\
    phi_s$ID$.addConst(1.) # shrink slightly\n\
    phi_s$ID$.join(phiParts_s$ID$)\n\
    extrapolateLsSimple(phi=phi_s$ID$, distance=narrowBandWidth_s$ID$+2, inside=True)\n\
    extrapolateLsSimple(phi=phi_s$ID$, distance=3)\n\
    phi_s$ID$.setBoundNeumann(boundaryWidth_s$ID$) # make sure no particles are placed at outer boundary\n\
    \n\
    if doOpen_s$ID$:\n\
        resetOutflow(flags=flags_s$ID$, phi=phi_s$ID$, parts=pp_s$ID$, index=gpi_s$ID$, indexSys=pindex_s$ID$)\n\
        if using_sndparts_s$ID$:\n\
            resetOutflow(flags=flags_s$ID$, parts=ppSnd_s$ID$)\n\
    flags_s$ID$.updateFromLevelset(phi_s$ID$)\n\
    \n\
    # combine particles velocities with advected grid velocities\n\
    mapPartsToMAC(vel=velParts_s$ID$, flags=flags_s$ID$, velOld=velOld_s$ID$, parts=pp_s$ID$, partVel=pVel_pp$ID$, weight=mapWeights_s$ID$)\n\
    extrapolateMACFromWeight(vel=velParts_s$ID$, distance=2, weight=mapWeights_s$ID$)\n\
    combineGridVel(vel=velParts_s$ID$, weight=mapWeights_s$ID$, combineVel=vel_s$ID$, phi=phi_s$ID$, narrowBand=combineBandWidth_s$ID$, thresh=0)\n\
    velOld_s$ID$.copyFrom(vel_s$ID$)\n\
    \n\
    # forces & pressure solve\n\
    addGravity(flags=flags_s$ID$, vel=vel_s$ID$, gravity=gravity_s$ID$)\n\
    addForceField(flags=flags_s$ID$, vel=vel_s$ID$, force=forces_s$ID$)\n\
    \n\
    if using_obstacle_s$ID$:\n\
        mantaMsg('Extrapolating object velocity')\n\
        # ensure velocities inside of obs object, slightly add obvels outside of obs object\n\
        extrapolateVec3Simple(vel=obvelC_s$ID$, phi=phiObsIn_s$ID$, distance=int(res_s$ID$/2), inside=True)\n\
        extrapolateVec3Simple(vel=obvelC_s$ID$, phi=phiObsIn_s$ID$, distance=3, inside=False)\n\
        resampleVec3ToMac(source=obvelC_s$ID$, target=obvel_s$ID$)\n\
    \n\
    if using_guiding_s$ID$:\n\
        mantaMsg('Extrapolating guiding velocity')\n\
        # ensure velocities inside of guiding object, slightly add guiding vels outside of object too\n\
        extrapolateVec3Simple(vel=guidevelC_s$ID$, phi=phiGuideIn_s$ID$, distance=int(res_s$ID$/2), inside=True)\n\
        extrapolateVec3Simple(vel=guidevelC_s$ID$, phi=phiGuideIn_s$ID$, distance=4, inside=False)\n\
        resampleVec3ToMac(source=guidevelC_s$ID$, target=guidevel_s$ID$)\n\
    \n\
    extrapolateMACSimple(flags=flags_s$ID$, vel=vel_s$ID$, distance=2, phiObs=phiObs_s$ID$, intoObs=True)\n\
    setWallBcs(flags=flags_s$ID$, vel=vel_s$ID$, obvel=obvel_s$ID$ if using_obstacle_s$ID$ else 0, phiObs=phiObs_s$ID$, fractions=fractions_s$ID$)\n\
    \n\
    if using_guiding_s$ID$:\n\
        mantaMsg('Guiding and pressure')\n\
        weightGuide_s$ID$.addConst(alpha_s$ID$)\n\
        PD_fluid_guiding(vel=vel_s$ID$, velT=guidevel_s$ID$, flags=flags_s$ID$, phi=phi_s$ID$, fractions=fractions_s$ID$, weight=weightGuide_s$ID$, blurRadius=beta_s$ID$, pressure=pressure_s$ID$, tau=tau_s$ID$, sigma=sigma_s$ID$, theta=theta_s$ID$, zeroPressureFixing=not doOpen_s$ID$)\n\
    else:\n\
        mantaMsg('Pressure')\n\
        solvePressure(flags=flags_s$ID$, vel=vel_s$ID$, pressure=pressure_s$ID$, phi=phi_s$ID$, fractions=fractions_s$ID$)\n\
    \n\
    extrapolateMACSimple(flags=flags_s$ID$, vel=vel_s$ID$, distance=4, phiObs=phiObs_s$ID$, intoObs=True)\n\
    setWallBcs(flags=flags_s$ID$, vel=vel_s$ID$, obvel=obvel_s$ID$ if using_obstacle_s$ID$ else 0, phiObs=phiObs_s$ID$, fractions=fractions_s$ID$)\n\
    \n\
    if (dim_s$ID$==3):\n\
        # mis-use phiParts as temp grid to close the mesh\n\
        phiParts_s$ID$.copyFrom(phi_s$ID$)\n\
        phiParts_s$ID$.setBound(0.5,0)\n\
        phiParts_s$ID$.createMesh(mesh_s$ID$)\n\
    \n\
    # Create interpolated version of original phi grid for later use in (optional) high-res step\n\
    if using_highres_s$ID$:\n\
        interpolateGrid(target=phi_xl$ID$, source=phiParts_s$ID$)\n\
    \n\
    extrapolateMACSimple(flags=flags_s$ID$, vel=vel_s$ID$, distance=(int(maxVel_s$ID$*1.25 )) ) # TODO (sebbas): extrapolation because of no fractions\n\
    # set source grids for resampling, used in adjustNumber!\n\
    pVel_pp$ID$.setSource(vel_s$ID$, isMAC=True)\n\
    adjustNumber(parts=pp_s$ID$, vel=vel_s$ID$, flags=flags_s$ID$, minParticles=minParticles_s$ID$, maxParticles=maxParticles_s$ID$, phi=phi_s$ID$, exclude=phiObs_s$ID$, radiusFactor=radiusFactor_s$ID$, narrowBand=adjustedNarrowBandWidth_s$ID$)\n\
    flipVelocityUpdate(vel=vel_s$ID$, velOld=velOld_s$ID$, flags=flags_s$ID$, parts=pp_s$ID$, partVel=pVel_pp$ID$, flipRatio=0.97)\n";

const std::string liquid_step_high = "\n\
def liquid_step_high_$ID$():\n\
    mantaMsg('Liquid step high')\n\
    pp_xl$ID$.readParticles(pp_s$ID$)\n\
    \n\
    # create surface\n\
    gridParticleIndex(parts=pp_xl$ID$, flags=flags_xl$ID$, indexSys=pindex_xl$ID$, index=gpi_xl$ID$)\n\
    averagedParticleLevelset(pp_xl$ID$, pindex_xl$ID$, flags_xl$ID$, gpi_xl$ID$, phiParts_xl$ID$, radiusFactor_s$ID$ , 1, 1)\n\
    phi_xl$ID$.join(phiParts_xl$ID$)\n\
    \n\
    phi_xl$ID$.createMesh(mesh_xl$ID$)\n";

//////////////////////////////////////////////////////////////////////
// IMPORT / EXPORT
//////////////////////////////////////////////////////////////////////

const std::string liquid_save_mesh_low = "\n\
def save_mesh_low_$ID$(path):\n\
    mesh_s$ID$.save(path)\n";

const std::string liquid_save_mesh_high = "\n\
def save_mesh_high_$ID$(path):\n\
    mesh_xl$ID$.save(path)\n";

const std::string liquid_save_particles_low = "\n\
def save_particles_low_$ID$(path):\n\
    pp_s$ID$.save(path)\n";

const std::string liquid_save_particle_velocities = "\n\
def save_particles_velocities_$ID$(path):\n\
    pVel_pp$ID$.save(path)\n";

const std::string liquid_import_low = "\n\
def load_liquid_data_low_$ID$(path):\n\
    flags_s$ID$.load(path + '_flags.uni')\n\
    \n\
    phiParts_s$ID$.load(path + '_phiParts.uni')\n\
    phi_s$ID$.load(path + '_phi.uni')\n\
    phiIn_s$ID$.load(path + '_phiIn.uni')\n\
    phiObs_s$ID$.load(path + '_phiObs.uni')\n\
    phiOut_s$ID$.load(path + '_phiOut.uni')\n\
    phiOutIn_s$ID$.load(path + '_phiOutIn.uni')\n\
    if fractions_s$ID$: fractions_s$ID$.load(path + '_fractions.uni') # TODO (sebbas)\n\
    pressure_s$ID$.load(path + '_pressure.uni')\n\
    \n\
    vel_s$ID$.load(path + '_vel.uni')\n\
    velOld_s$ID$.load(path + '_velOld.uni')\n\
    velParts_s$ID$.load(path + '_velParts.uni')\n\
    mapWeights_s$ID$.load(path + '_mapWeights.uni')\n\
    \n\
    x_vel_s$ID$.load(path + '_x_vel.uni')\n\
    y_vel_s$ID$.load(path + '_y_vel.uni')\n\
    z_vel_s$ID$.load(path + '_z_vel.uni')\n\
    \n\
    pp_s$ID$.load(path + '_pp.uni')\n\
    pVel_pp$ID$.load(path + '_pVel.uni')\n\
    \n\
    gpi_s$ID$.load(path + '_gpi.uni')\n";

const std::string liquid_import_high = "\n\
def load_liquid_data_high_$ID$(path):\n\
    flags_xl$ID$.load(path + '_flags_xl.uni')\n\
    \n\
    phiParts_xl$ID$.load(path + '_phiParts_xl.uni')\n\
    phi_xl$ID$.load(path + '_phi_xl.uni')\n\
    \n\
    pp_xl$ID$.load(path + '_pp_xl.uni')\n";

const std::string liquid_export_low = "\n\
def save_liquid_data_low_$ID$(path):\n\
    flags_s$ID$.save(path + '_flags.uni')\n\
    \n\
    phiParts_s$ID$.save(path + '_phiParts.uni')\n\
    phi_s$ID$.save(path + '_phi.uni')\n\
    phiIn_s$ID$.save(path + '_phiIn.uni')\n\
    phiObs_s$ID$.save(path + '_phiObs.uni')\n\
    phiOut_s$ID$.save(path + '_phiOut.uni')\n\
    phiOutIn_s$ID$.save(path + '_phiOutIn.uni')\n\
    if fractions_s$ID$: fractions_s$ID$.save(path + '_fractions.uni') # TODO (sebbas)\n\
    pressure_s$ID$.save(path + '_pressure.uni')\n\
    \n\
    vel_s$ID$.save(path + '_vel.uni')\n\
    velOld_s$ID$.save(path + '_velOld.uni')\n\
    velParts_s$ID$.save(path + '_velParts.uni')\n\
    mapWeights_s$ID$.save(path + '_mapWeights.uni')\n\
    \n\
    x_vel_s$ID$.save(path + '_x_vel.uni')\n\
    y_vel_s$ID$.save(path + '_y_vel.uni')\n\
    z_vel_s$ID$.save(path + '_z_vel.uni')\n\
    \n\
    pp_s$ID$.save(path + '_pp.uni')\n\
    pVel_pp$ID$.save(path + '_pVel.uni')\n\
    \n\
    gpi_s$ID$.save(path + '_gpi.uni')\n";

const std::string liquid_export_high = "\n\
def save_liquid_data_high_$ID$(path):\n\
    flags_xl$ID$.save(path + '_flags_xl.uni')\n\
    \n\
    phiParts_xl$ID$.save(path + '_phiParts_xl.uni')\n\
    phi_xl$ID$.save(path + '_phi_xl.uni')\n\
    \n\
    pp_xl$ID$.save(path + '_pp_xl.uni')\n";

//////////////////////////////////////////////////////////////////////
// DESTRUCTION
//////////////////////////////////////////////////////////////////////

const std::string liquid_delete_grids_low = "\n\
mantaMsg('Deleting lowres grids, mesh, particlesystem')\n\
if 'flags_s$ID$'      in globals() : del flags_s$ID$\n\
if 'phiParts_s$ID$'   in globals() : del phiParts_s$ID$\n\
if 'phi_s$ID$'        in globals() : del phi_s$ID$\n\
if 'phiIn_s$ID$'      in globals() : del phiIn_s$ID$\n\
if 'phiOut_s$ID$'     in globals() : del phiOut_s$ID$\n\
if 'phiOutIn_s$ID$'   in globals() : del phiOutIn_s$ID$\n\
if 'pressure_s$ID$'   in globals() : del pressure_s$ID$\n\
if 'vel_s$ID$'        in globals() : del vel_s$ID$\n\
if 'x_vel_s$ID$'      in globals() : del x_vel_s$ID$\n\
if 'y_vel_s$ID$'      in globals() : del y_vel_s$ID$\n\
if 'z_vel_s$ID$'      in globals() : del z_vel_s$ID$\n\
if 'velOld_s$ID$'     in globals() : del velOld_s$ID$\n\
if 'velParts_s$ID$'   in globals() : del velParts_s$ID$\n\
if 'mapWeights_s$ID$' in globals() : del mapWeights_s$ID$\n\
if 'pp_s$ID$'         in globals() : del pp_s$ID$\n\
if 'pVel_pp$ID$'      in globals() : del pVel_pp$ID$\n\
if 'mesh_s$ID$'       in globals() : del mesh_s$ID$\n\
if 'pindex_s$ID$'     in globals() : del pindex_s$ID$\n\
if 'gpi_s$ID$'        in globals() : del gpi_s$ID$\n\
if 'forces_s$ID$'     in globals() : del forces_s$ID$\n\
if 'x_force_s$ID$'    in globals() : del x_force_s$ID$\n\
if 'y_force_s$ID$'    in globals() : del y_force_s$ID$\n\
if 'z_force_s$ID$'    in globals() : del z_force_s$ID$\n\
if 'phiObs_s$ID$'     in globals() : del phiObs_s$ID$\n\
if 'fractions_s$ID$'  in globals() : del fractions_s$ID$\n";

const std::string liquid_delete_grids_high = "\n\
mantaMsg('Deleting highres grids, mesh, particlesystem')\n\
if 'flags_xl$ID$'    in globals() : del flags_xl$ID$\n\
if 'phiParts_xl$ID$' in globals() : del phiParts_xl$ID$\n\
if 'phi_xl$ID$'      in globals() : del phi_xl$ID$\n\
if 'pp_xl$ID$'       in globals() : del pp_xl$ID$\n\
if 'mesh_xl$ID$'     in globals() : del mesh_xl$ID$\n\
if 'pindex_xl$ID$'   in globals() : del pindex_xl$ID$\n\
if 'gpi_xl$ID$'      in globals() : del gpi_xl$ID$\n";

const std::string liquid_delete_variables_low = "\n\
mantaMsg('Deleting lowres liquid variables')\n\
if 'narrowBandWidth_s$ID$'  in globals() : del narrowBandWidth_s$ID$\n\
if 'combineBandWidth_s$ID$' in globals() : del combineBandWidth_s$ID$\n\
if 'minParticles_s$ID$'     in globals() : del minParticles_s$ID$\n\
if 'maxParticles_s$ID$'     in globals() : del maxParticles_s$ID$\n\
if 'particleNumber_s$ID$'   in globals() : del particleNumber_s$ID$\n\
if 'maxVel_s$ID$'           in globals() : del maxVel_s$ID$\n";

const std::string liquid_delete_variables_high = "\n\
mantaMsg('Deleting highres liquid variables')\n";

//////////////////////////////////////////////////////////////////////
// STANDALONE MODE
//////////////////////////////////////////////////////////////////////

const std::string liquid_standalone_load = "\n\
# import *.uni files\n\
path_prefix_$ID$ = '$MANTA_EXPORT_PATH$'\n\
load_liquid_data_low_$ID$(path_prefix_$ID$)\n\
if using_highres_s$ID$:\n\
    load_liquid_data_high_$ID$(path_prefix_$ID$)\n";
