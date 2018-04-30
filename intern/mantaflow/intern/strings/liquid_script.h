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

const std::string liquid_variables = "\n\
mantaMsg('Liquid variables low')\n\
narrowBandWidth_s$ID$         = 3\n\
combineBandWidth_s$ID$        = narrowBandWidth_s$ID$ - 1\n\
adjustedNarrowBandWidth_s$ID$ = $PARTICLE_BAND_WIDTH$ # only used in adjustNumber to control band width\n\
\n\
particleNumber_s$ID$ = $PARTICLE_NUMBER$\n\
minParticles_s$ID$   = $PARTICLE_MINIMUM$\n\
maxParticles_s$ID$   = $PARTICLE_MAXIMUM$\n\
radiusFactor_s$ID$   = $PARTICLE_RADIUS$\n\
smoothenUpper_s$ID$  = $MESH_SMOOTHEN_UPPER$\n\
smoothenLower_s$ID$  = $MESH_SMOOTHEN_LOWER$\n\
randomness_s$ID$     = $PARTICLE_RANDOMNESS$\n\
surfaceTension_s$ID$ = $LIQUID_SURFACE_TENSION$\n";

const std::string liquid_variables_mesh = "\n\
mantaMsg('Liquid variables high')\n";

//////////////////////////////////////////////////////////////////////
// GRIDS & MESH & PARTICLESYSTEM
//////////////////////////////////////////////////////////////////////

const std::string liquid_alloc = "\n\
mantaMsg('Liquid alloc low')\n\
phiParts_s$ID$   = s$ID$.create(LevelsetGrid)\n\
phi_s$ID$        = s$ID$.create(LevelsetGrid)\n\
phiTmp_s$ID$     = s$ID$.create(LevelsetGrid)\n\
phiIn_s$ID$      = s$ID$.create(LevelsetGrid)\n\
curvature_s$ID$  = s$ID$.create(RealGrid)\n\
\n\
fractions_s$ID$  = s$ID$.create(MACGrid) # TODO (sebbas): disabling fractions for now - not fracwallbcs not supporting obvels yet\n\
\n\
velOld_s$ID$     = s$ID$.create(MACGrid)\n\
velParts_s$ID$   = s$ID$.create(MACGrid)\n\
mapWeights_s$ID$ = s$ID$.create(MACGrid)\n\
\n\
pp_s$ID$         = s$ID$.create(BasicParticleSystem)\n\
pVel_pp$ID$      = pp_s$ID$.create(PdataVec3)\n\
\n\
# Acceleration data for particle nbs\n\
pindex_s$ID$     = s$ID$.create(ParticleIndexSystem)\n\
gpi_s$ID$        = s$ID$.create(IntGrid)\n";

const std::string liquid_alloc_mesh = "\n\
mantaMsg('Liquid alloc high')\n\
phiParts_xl$ID$ = sm$ID$.create(LevelsetGrid)\n\
phi_xl$ID$      = sm$ID$.create(LevelsetGrid)\n\
pp_xl$ID$       = sm$ID$.create(BasicParticleSystem)\n\
flags_xl$ID$    = sm$ID$.create(FlagGrid)\n\
mesh_xl$ID$     = sm$ID$.create(Mesh)\n\
\n\
# Acceleration data for particle nbs\n\
pindex_xl$ID$  = sm$ID$.create(ParticleIndexSystem)\n\
gpi_xl$ID$     = sm$ID$.create(IntGrid)\n";

const std::string liquid_init_phi = "\n\
phi_s$ID$.initFromFlags(flags_s$ID$)\n\
phiIn_s$ID$.initFromFlags(flags_s$ID$)\n";

//////////////////////////////////////////////////////////////////////
// PRE / POST STEP
//////////////////////////////////////////////////////////////////////

const std::string liquid_pre_step = "\n\
def liquid_pre_step_$ID$():\n\
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
    copyRealToVec3(sourceX=x_force_s$ID$, sourceY=y_force_s$ID$, sourceZ=z_force_s$ID$, target=forces_s$ID$)\n";

const std::string liquid_post_step = "\n\
def liquid_post_step_$ID$():\n\
    forces_s$ID$.clear()\n\
    if using_guiding_s$ID$:\n\
        weightGuide_s$ID$.clear()\n\
    if using_invel_s$ID$:\n\
        invel_s$ID$.clear()\n\
    \n\
#    phiIn_s$ID$.setConst(9999)\n\
#    phiOutIn_s$ID$.setConst(9999)\n\
    \n\
    #copyVec3ToReal(source=vel_s$ID$, targetX=x_vel_s$ID$, targetY=y_vel_s$ID$, targetZ=z_vel_s$ID$)\n\
    #x_vel_s$ID$.multConst( 1.0/Real(gs_s$ID$.x) )\n\
    #y_vel_s$ID$.multConst( 1.0/Real(gs_s$ID$.y) )\n\
    #z_vel_s$ID$.multConst( 1.0/Real(gs_s$ID$.z) )\n";

//////////////////////////////////////////////////////////////////////
// STEP FUNCTIONS
//////////////////////////////////////////////////////////////////////

const std::string liquid_adaptive_step = "\n\
def liquid_adaptive_step_$ID$(framenr):\n\
    mantaMsg('Manta step, frame ' + str(framenr))\n\
    \n\
    # time params are animatable\n\
    s$ID$.frameLength = dt0_s$ID$ \n\
    s$ID$.cfl = cfl_cond_s$ID$\n\
    \n\
    liquid_pre_step_$ID$()\n\
    \n\
    flags_s$ID$.initDomain(boundaryWidth=boundaryWidth_s$ID$, phiWalls=phiObs_s$ID$, outflow=boundConditions_s$ID$)\n\
    \n\
    if using_obstacle_s$ID$:\n\
        phiObs_s$ID$.join(phiObsIn_s$ID$)\n\
    phi_s$ID$.join(phiIn_s$ID$)\n\
    if using_obstacle_s$ID$:\n\
        phi_s$ID$.subtract(phiObsIn_s$ID$)\n\
    \n\
    phiOut_s$ID$.join(phiOutIn_s$ID$)\n\
    \n\
    updateFractions(flags=flags_s$ID$, phiObs=phiObs_s$ID$, fractions=fractions_s$ID$, boundaryWidth=boundaryWidth_s$ID$) # TODO (sebbas): uncomment for fraction support\n\
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
    mantaMsg('Low step / s$ID$.frame: ' + str(s$ID$.frame))\n\
    liquid_step_$ID$()\n\
    \n\
    s$ID$.step()\n\
    \n\
    liquid_post_step_$ID$()\n";

const std::string liquid_step = "\n\
def liquid_step_$ID$():\n\
    mantaMsg('Liquid step low')\n\
    \n\
    mantaMsg('Advecting particles')\n\
    pp_s$ID$.advectInGrid(flags=flags_s$ID$, vel=vel_s$ID$, integrationMode=IntRK4, deleteInObstacle=False, stopInObstacle=False)\n\
    \n\
    mantaMsg('Pushing particles out of obstacles')\n\
    pushOutofObs(parts=pp_s$ID$, flags=flags_s$ID$, phiObs=phiObs_s$ID$)\n\
    \n\
    mantaMsg('Advecting phi')\n\
    advectSemiLagrange(flags=flags_s$ID$, vel=vel_s$ID$, grid=phi_s$ID$, order=1) # first order is usually enough\n\
    mantaMsg('Advecting velocity')\n\
    advectSemiLagrange(flags=flags_s$ID$, vel=vel_s$ID$, grid=vel_s$ID$, order=2, openBounds=doOpen_s$ID$, boundaryWidth=boundaryWidth_s$ID$)\n\
    \n\
    phiTmp_s$ID$.copyFrom(phi_s$ID$) # save original phi for later use in mesh creation\n\
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
            resetOutflow(flags=flags_s$ID$, parts=ppSnd_sp$ID$)\n\
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
    \n\
    # vel diffusion / viscosity!\n\
    if viscosity_s$ID$ > 0.:\n\
        mantaMsg('Viscosity')\n\
        # diffusion param for solve = const * dt / dx^2\n\
        alphaV = viscosity_s$ID$ * s$ID$.timestep * float(res_s$ID$*res_s$ID$)\n\
        setWallBcs(flags=flags_s$ID$, vel=vel_s$ID$, obvel=obvel_s$ID$ if using_obstacle_s$ID$ else 0, phiObs=phiObs_s$ID$, fractions=fractions_s$ID$)\n\
        cgSolveDiffusion(flags_s$ID$, vel_s$ID$, alphaV)\n\
    \n\
    setWallBcs(flags=flags_s$ID$, vel=vel_s$ID$, obvel=obvel_s$ID$ if using_obstacle_s$ID$ else 0, phiObs=phiObs_s$ID$, fractions=fractions_s$ID$)\n\
    \n\
    mantaMsg('Calculating curvature')\n\
    getLaplacian(laplacian=curvature_s$ID$, grid=phi_s$ID$)\n\
    \n\
    if using_guiding_s$ID$:\n\
        mantaMsg('Guiding and pressure')\n\
        weightGuide_s$ID$.addConst(alpha_s$ID$)\n\
        PD_fluid_guiding(vel=vel_s$ID$, velT=guidevel_s$ID$, flags=flags_s$ID$, phi=phi_s$ID$, curv=curvature_s$ID$, surfTens=surfaceTension_s$ID$, fractions=fractions_s$ID$, weight=weightGuide_s$ID$, blurRadius=beta_s$ID$, pressure=pressure_s$ID$, tau=tau_s$ID$, sigma=sigma_s$ID$, theta=theta_s$ID$, zeroPressureFixing=not doOpen_s$ID$)\n\
    else:\n\
        mantaMsg('Pressure')\n\
        solvePressure(flags=flags_s$ID$, vel=vel_s$ID$, pressure=pressure_s$ID$, phi=phi_s$ID$, curv=curvature_s$ID$, surfTens=surfaceTension_s$ID$, fractions=fractions_s$ID$)\n\
    \n\
    extrapolateMACSimple(flags=flags_s$ID$, vel=vel_s$ID$, distance=4, phiObs=phiObs_s$ID$, intoObs=True)\n\
    setWallBcs(flags=flags_s$ID$, vel=vel_s$ID$, obvel=obvel_s$ID$ if using_obstacle_s$ID$ else 0, phiObs=phiObs_s$ID$, fractions=fractions_s$ID$)\n\
    \n\
    #extrapolateMACSimple(flags=flags_s$ID$, vel=vel_s$ID$, distance=(int(maxVel_s$ID$*1.25 )) ) # TODO (sebbas): extrapolation because of no fractions\n\
    # set source grids for resampling, used in adjustNumber!\n\
    pVel_pp$ID$.setSource(vel_s$ID$, isMAC=True)\n\
    adjustNumber(parts=pp_s$ID$, vel=vel_s$ID$, flags=flags_s$ID$, minParticles=minParticles_s$ID$, maxParticles=maxParticles_s$ID$, phi=phi_s$ID$, exclude=phiObs_s$ID$, radiusFactor=1., narrowBand=adjustedNarrowBandWidth_s$ID$)\n\
    flipVelocityUpdate(vel=vel_s$ID$, velOld=velOld_s$ID$, flags=flags_s$ID$, parts=pp_s$ID$, partVel=pVel_pp$ID$, flipRatio=0.97)\n";

const std::string liquid_step_particles = "\n\
def liquid_step_particles_$ID$():\n\
    mantaMsg('Sampling snd particles')\n\
    interpolateMACGrid(target=vel_sp$ID$, source=vel_s$ID$)\n\
    interpolateGrid(target=phi_sp$ID$, source=phi_s$ID$)\n\
    interpolateGrid(target=phiIn_sp$ID$, source=phiIn_s$ID$)\n\
    \n\
    # Manually recreate phiObs levelset to ensure correct values at border. Then use that to create upscaled flags grid\n\
    flags_sp$ID$.initDomain(boundaryWidth=boundaryWidth_s$ID$, phiWalls=phiObs_sp$ID$, outflow=boundConditions_s$ID$)\n\
    if using_obstacle_s$ID$:\n\
        interpolateGrid(target=phiObsIn_sp$ID$, source=phiObsIn_s$ID$)\n\
        phiObs_sp$ID$.join(phiObsIn_sp$ID$)\n\
    setObstacleFlags(flags=flags_sp$ID$, phiObs=phiObs_sp$ID$) #, phiOut=phiOut_sp$ID$) # TODO (sebbas): outflow for snd parts\n\
    \n\
    sampleSndParts(phi=phi_sp$ID$, phiIn=phiIn_sp$ID$, flags=flags_sp$ID$, vel=vel_sp$ID$, parts=ppSnd_sp$ID$, type=$SNDPARTICLE_TYPES$, amountDroplet=$SNDPARTICLE_DROPLET_AMOUNT$, amountFloater=$SNDPARTICLE_FLOATER_AMOUNT$, amountTracer=$SNDPARTICLE_TRACER_AMOUNT$, thresholdDroplet=$SNDPARTICLE_DROPLET_THRESH$)\n\
    mantaMsg('Updating snd particle data (velocity, life count)')\n\
    updateSndParts(phi=phi_sp$ID$, flags=flags_sp$ID$, vel=vel_sp$ID$, gravity=gravity_s$ID$, parts=ppSnd_sp$ID$, partVel=pVelSnd_pp$ID$, partLife=pLifeSnd_pp$ID$, riseBubble=$SNDPARTICLE_BUBBLE_RISE$, lifeDroplet=$SNDPARTICLE_DROPLET_LIFE$, lifeBubble=$SNDPARTICLE_BUBBLE_LIFE$, lifeFloater=$SNDPARTICLE_FLOATER_LIFE$, lifeTracer=$SNDPARTICLE_TRACER_LIFE$)\n\
    mantaMsg('Adjusting snd particles')\n\
    pushOutofObs(parts=ppSnd_sp$ID$, flags=flags_sp$ID$, phiObs=phiObs_sp$ID$, shift=1.0)\n\
    adjustSndParts(parts=ppSnd_sp$ID$, flags=flags_sp$ID$, phi=phi_sp$ID$, partVel=pVelSnd_pp$ID$, partLife=pLifeSnd_pp$ID$, maxDroplet=$SNDPARTICLE_DROPLET_MAX$, maxBubble=$SNDPARTICLE_BUBBLE_MAX$, maxFloater=$SNDPARTICLE_FLOATER_MAX$, maxTracer=$SNDPARTICLE_TRACER_MAX$)\n";

//////////////////////////////////////////////////////////////////////
// IMPORT
//////////////////////////////////////////////////////////////////////

const std::string liquid_load_data = "\n\
def liquid_load_data_$ID$(path, framenr, withParticles):\n\
    mantaMsg('Liquid load data')\n\
    framenr = fluid_cache_get_framenr_formatted_$ID$(framenr)\n\
    phi_s$ID$.load(os.path.join(path, 'phi_' + framenr +'.uni'))\n\
    phiTmp_s$ID$.load(os.path.join(path, 'phiTmp_' + framenr +'.uni'))\n\
    phiIn_s$ID$.load(os.path.join(path, 'phiIn_' + framenr + '.uni'))\n\
    phiParts_s$ID$.load(os.path.join(path, 'phiParts_' + framenr + '.uni'))\n\
    if withParticles:\n\
        pp_s$ID$.load(os.path.join(path, 'pp_' + framenr +'.uni'))\n\
        pVel_pp$ID$.load(os.path.join(path, 'pVel_' + framenr +'.uni'))\n";

const std::string liquid_load_mesh = "\n\
def liquid_load_mesh_$ID$(path, framenr):\n\
    mantaMsg('Liquid load mesh')\n\
    framenr = fluid_cache_get_framenr_formatted_$ID$(framenr)\n\
    mesh_xl$ID$.load(os.path.join(path, 'liquid_mesh_' + framenr + '.bobj.gz'))\n\
    pp_xl$ID$.load(os.path.join(path, 'pp_xl_' + framenr +'.uni'))\n";

const std::string liquid_load_particles = "\n\
def liquid_load_particles_$ID$(path, framenr):\n\
    mantaMsg('Liquid load particles')\n\
    framenr = fluid_cache_get_framenr_formatted_$ID$(framenr)\n\
    ppSnd_sp$ID$.load(os.path.join(path, 'ppSnd_' + framenr + '.uni'))\n\
    pVelSnd_pp$ID$.load(os.path.join(path, 'pVelSnd_' + framenr + '.uni'))\n\
    pLifeSnd_pp$ID$.load(os.path.join(path, 'pLifeSnd_' + framenr + '.uni'))\n";

//////////////////////////////////////////////////////////////////////
// EXPORT
//////////////////////////////////////////////////////////////////////

const std::string liquid_save_data = "\n\
def liquid_save_data_$ID$(path, framenr):\n\
    mantaMsg('Liquid save data low')\n\
    framenr = fluid_cache_get_framenr_formatted_$ID$(framenr)\n\
    phi_s$ID$.save(os.path.join(path, 'phi_' + framenr +'.uni'))\n\
    phiTmp_s$ID$.save(os.path.join(path, 'phiTmp_' + framenr +'.uni'))\n\
    phiIn_s$ID$.save(os.path.join(path, 'phiIn_' + framenr + '.uni'))\n\
    phiParts_s$ID$.save(os.path.join(path, 'phiParts_' + framenr + '.uni'))\n\
    pp_s$ID$.save(os.path.join(path, 'pp_' + framenr +'.uni'))\n\
    pVel_pp$ID$.save(os.path.join(path, 'pVel_' + framenr +'.uni'))\n";

const std::string liquid_save_mesh = "\n\
def liquid_save_mesh_$ID$(path, framenr):\n\
    mantaMsg('Liquid save mesh high')\n\
    framenr = fluid_cache_get_framenr_formatted_$ID$(framenr)\n\
    \n\
    interpolateGrid(target=phi_xl$ID$, source=phiTmp_s$ID$) # mis-use phiParts as temp grid\n\
    \n\
    # create surface\n\
    pp_xl$ID$.readParticles(pp_s$ID$)\n\
    gridParticleIndex(parts=pp_xl$ID$, flags=flags_xl$ID$, indexSys=pindex_xl$ID$, index=gpi_xl$ID$)\n\
    averagedParticleLevelset(pp_xl$ID$, pindex_xl$ID$, flags_xl$ID$, gpi_xl$ID$, phiParts_xl$ID$, radiusFactor_s$ID$, 1, 1)\n\
#    unionParticleLevelset(pp_xl$ID$, pindex_xl$ID$, flags_xl$ID$, gpi_xl$ID$, phiParts_xl$ID$, radiusFactor_s$ID$)\n\
    \n\
    phi_xl$ID$.addConst(1.) # shrink slightly\n\
    phi_xl$ID$.join(phiParts_xl$ID$)\n\
    extrapolateLsSimple(phi=phi_xl$ID$, distance=narrowBandWidth_s$ID$+2, inside=True)\n\
    extrapolateLsSimple(phi=phi_xl$ID$, distance=3)\n\
    phi_xl$ID$.setBoundNeumann(boundaryWidth_s$ID$) # make sure no particles are placed at outer boundary\n\
    \n\
    phi_xl$ID$.setBound(0.5,int(((upres_sm$ID$)*2)-2) )\n\
    phi_xl$ID$.createMesh(mesh_xl$ID$)\n\
    mesh_xl$ID$.save(os.path.join(path, 'liquid_mesh_' + framenr + '.bobj.gz'))\n\
    pp_xl$ID$.save(os.path.join(path, 'pp_xl_' + framenr +'.uni'))\n";

const std::string liquid_save_particles = "\n\
def liquid_save_particles_$ID$(path, framenr):\n\
    mantaMsg('Liquid save particles low')\n\
    framenr = fluid_cache_get_framenr_formatted_$ID$(framenr)\n\
    ppSnd_sp$ID$.save(os.path.join(path, 'ppSnd_' + framenr + '.uni'))\n\
    pVelSnd_pp$ID$.save(os.path.join(path, 'pVelSnd_' + framenr + '.uni'))\n\
    pLifeSnd_pp$ID$.save(os.path.join(path, 'pLifeSnd_' + framenr + '.uni'))\n";

const std::string liquid_save_particles_high = "\n\
def liquid_save_particles_high$ID$(path, framenr):\n\
    mantaMsg('Liquid save particles high')\n\
    # Nothing to do here yet!\n";

//////////////////////////////////////////////////////////////////////
// STANDALONE MODE
//////////////////////////////////////////////////////////////////////

const std::string liquid_standalone_load = "\n\
# import *.uni files\n\
path_prefix_$ID$ = '$MANTA_EXPORT_PATH$'\n\
load_liquid_data_$ID$(path_prefix_$ID$)\n\
if using_highres_s$ID$:\n\
    load_liquid_data_high_$ID$(path_prefix_$ID$)\n";
