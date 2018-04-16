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
surfaceTension_s$ID$ = $LIQUID_SURFACE_TENSION$\n\
maxVel_s$ID$         = 1 # just declared here, do not set\n\
tauMin_wc_s$ID$      = $SNDPARTICLE_TAU_MIN_WC$\n\
tauMax_wc_s$ID$      = $SNDPARTICLE_TAU_MAX_WC$\n\
tauMin_ta_s$ID$      = $SNDPARTICLE_TAU_MIN_TA$\n\
tauMax_ta_s$ID$      = $SNDPARTICLE_TAU_MAX_TA$\n\
tauMin_k_s$ID$       = $SNDPARTICLE_TAU_MIN_K$\n\
tauMax_k_s$ID$       = $SNDPARTICLE_TAU_MAX_K$\n\
k_wc_s$ID$           = $SNDPARTICLE_K_WC$\n\
k_ta_s$ID$           = $SNDPARTICLE_K_TA$\n\
k_b_s$ID$            = $SNDPARTICLE_K_B$\n\
k_d_s$ID$            = $SNDPARTICLE_K_D$\n\
lMin_s$ID$           = $SNDPARTICLE_L_MIN$\n\
lMax_s$ID$           = $SNDPARTICLE_L_MAX$\n\
c_s_s$ID$            = 0.4   # classification constant for snd parts\n\
c_b_s$ID$            = 0.77  # classification constant for snd parts\n\
scaleFromManta_s$ID$ = $FLUID_DOMAIN_SIZE$ / float(res_s$ID$) # resize factor for snd parts\n";

const std::string liquid_variables_high = "\n\
mantaMsg('Liquid variables high')\n";

//////////////////////////////////////////////////////////////////////
// GRIDS & MESH & PARTICLESYSTEM
//////////////////////////////////////////////////////////////////////

const std::string liquid_alloc_low = "\n\
mantaMsg('Liquid alloc low')\n\
phiParts_s$ID$   = s$ID$.create(LevelsetGrid)\n\
phi_s$ID$        = s$ID$.create(LevelsetGrid)\n\
phiIn_s$ID$      = s$ID$.create(LevelsetGrid)\n\
phiOut_s$ID$     = s$ID$.create(LevelsetGrid)\n\
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
mesh_s$ID$       = s$ID$.create(Mesh)\n\
\n\
# Acceleration data for particle nbs\n\
pindex_s$ID$     = s$ID$.create(ParticleIndexSystem)\n\
gpi_s$ID$        = s$ID$.create(IntGrid)\n\
\n\
normal_s$ID$ = s$ID$.create(VecGrid)\n\
neighborRatio_s$ID$ = s$ID$.create(RealGrid)\n\
trappedAir_s$ID$ = s$ID$.create(RealGrid)\n\
waveCrest_s$ID$ = s$ID$.create(RealGrid)\n\
kineticEnergy_s$ID$ = s$ID$.create(RealGrid)\n";

const std::string liquid_alloc_high = "\n\
mantaMsg('Liquid alloc high')\n\
phiParts_xl$ID$ = xl$ID$.create(LevelsetGrid)\n\
phi_xl$ID$      = xl$ID$.create(LevelsetGrid)\n\
pp_xl$ID$       = xl$ID$.create(BasicParticleSystem)\n\
mesh_xl$ID$     = xl$ID$.create(Mesh)\n\
\n\
# Acceleration data for particle nbs\n\
pindex_xl$ID$  = xl$ID$.create(ParticleIndexSystem)\n\
gpi_xl$ID$     = xl$ID$.create(IntGrid)\n\
\n\
normal_xl$ID$ = s$ID$.create(VecGrid)\n\
neighborRatio_xl$ID$ = s$ID$.create(RealGrid)\n\
trappedAir_xl$ID$ = s$ID$.create(RealGrid)\n\
waveCrest_xl$ID$ = s$ID$.create(RealGrid)\n\
kineticEnergy_xl$ID$ = s$ID$.create(RealGrid)\n";

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

const std::string liquid_post_step_low = "\n\
def liquid_post_step_low_$ID$():\n\
    forces_s$ID$.clear()\n\
    if using_guiding_s$ID$:\n\
        weightGuide_s$ID$.clear()\n\
    if using_invel_s$ID$:\n\
        invel_s$ID$.clear()\n\
    \n\
    phiIn_s$ID$.setConst(9999)\n\
    phiOut_s$ID$.setConst(9999)\n\
    #copyVec3ToReal(source=vel_s$ID$, targetX=x_vel_s$ID$, targetY=y_vel_s$ID$, targetZ=z_vel_s$ID$)\n\
    #x_vel_s$ID$.multConst( 1.0/Real(gs_s$ID$.x) )\n\
    #y_vel_s$ID$.multConst( 1.0/Real(gs_s$ID$.y) )\n\
    #z_vel_s$ID$.multConst( 1.0/Real(gs_s$ID$.z) )\n";

//////////////////////////////////////////////////////////////////////
// STEP FUNCTIONS
//////////////////////////////////////////////////////////////////////

const std::string liquid_adaptive_step_low = "\n\
def liquid_adaptive_step_low_$ID$(framenr):\n\
    mantaMsg('Manta step low, frame ' + str(framenr))\n\
    \n\
    # time params are animatable\n\
    s$ID$.frameLength = dt0_s$ID$ \n\
    s$ID$.cfl = cfl_cond_s$ID$\n\
    \n\
    mantaMsg('s.frame is ' + str(s$ID$.frame))\n\
    mantaMsg('s.timestep is ' + str(s$ID$.timestep))\n\
    mantaMsg('s.cfl is ' + str(s$ID$.cfl))\n\
    mantaMsg('s.frameLength is ' + str(s$ID$.frameLength))\n\
    \n\
    liquid_pre_step_low_$ID$()\n\
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
    liquid_post_step_low_$ID$()\n";

const std::string liquid_adaptive_step_high = "\n\
def liquid_adaptive_step_high_$ID$(framenr):\n\
    mantaMsg('Manta step high, frame ' + str(framenr))\n\
    \n\
    xl$ID$.frame = framenr\n\
    xl$ID$.timeTotal = xl$ID$.frame * dt0_s$ID$\n\
    last_frame_s$ID$ = xl$ID$.frame\n\
    \n\
    liquid_pre_step_low_$ID$()\n\
    \n\
    while xl$ID$.frame == last_frame_s$ID$:\n\
        \n\
        mantaMsg('xl.frame is ' + str(xl$ID$.frame))\n\
        \n\
        fluid_adapt_time_step_high_$ID$()\n\
        mantaMsg('High step / xl$ID$.frame: ' + str(xl$ID$.frame))\n\
        liquid_step_high_$ID$()\n\
        xl$ID$.step()\n";

const std::string liquid_step_low = "\n\
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
    adjustNumber(parts=pp_s$ID$, vel=vel_s$ID$, flags=flags_s$ID$, minParticles=minParticles_s$ID$, maxParticles=maxParticles_s$ID$, phi=phi_s$ID$, exclude=phiObs_s$ID$, radiusFactor=radiusFactor_s$ID$, narrowBand=adjustedNarrowBandWidth_s$ID$)\n\
    flipVelocityUpdate(vel=vel_s$ID$, velOld=velOld_s$ID$, flags=flags_s$ID$, parts=pp_s$ID$, partVel=pVel_pp$ID$, flipRatio=0.97)\n";

const std::string liquid_step_high = "\n\
def liquid_step_high_$ID$():\n\
    mantaMsg('Liquid step high')\n\
    # interpolateGrid(target=flags_xl$ID$, source=flags_s$ID$) # TODO (sebbas): needs interpolation?\n\
    pp_xl$ID$.readParticles(pp_s$ID$)\n\
    \n\
    # create surface\n\
    gridParticleIndex(parts=pp_xl$ID$, flags=flags_xl$ID$, indexSys=pindex_xl$ID$, index=gpi_xl$ID$)\n\
    averagedParticleLevelset(pp_xl$ID$, pindex_xl$ID$, flags_xl$ID$, gpi_xl$ID$, phiParts_xl$ID$, radiusFactor_s$ID$ , 1, 1)\n";

const std::string liquid_step_particles_low = "\n\
def liquid_step_particles_low_$ID$():\n\
    mantaMsg('Secondary particles step low')\n\
    if using_sndparts_s$ID$:\n\
        radius = 2 if $SNDPARTICLE_POTENTIAL_QUALITY_HIGH$ else 1\n\
        flipComputeSecondaryParticlePotentials(potTA=trappedAir_s$ID$, potWC=waveCrest_s$ID$, potKE=kineticEnergy_s$ID$, neighborRatio=neighborRatio_s$ID$, flags=flags_s$ID$, v=vel_s$ID$, normal=normal_s$ID$, phi=phi_s$ID$, radius=radius, tauMinTA=tauMin_ta_s$ID$, tauMaxTA=tauMax_ta_s$ID$, tauMinWC=tauMin_wc_s$ID$, tauMaxWC=tauMax_wc_s$ID$, tauMinKE=tauMin_k_s$ID$, tauMaxKE=tauMax_k_s$ID$, scaleFromManta=scaleFromManta_s$ID$)\n\
        flipSampleSecondaryParticles(mode='single', flags=flags_s$ID$, v=vel_s$ID$, pts_sec=ppSnd_s$ID$, v_sec=pVelSnd_pp$ID$, l_sec=pLifeSnd_pp$ID$, lMin=lMin_s$ID$, lMax=lMax_s$ID$, potTA=trappedAir_s$ID$, potWC=waveCrest_s$ID$, potKE=kineticEnergy_s$ID$, neighborRatio=neighborRatio_s$ID$, c_s=c_s_s$ID$, c_b=c_b_s$ID$, k_ta=k_ta_s$ID$, k_wc=k_wc_s$ID$, dt=s$ID$.frameLength)\n\
        flipUpdateSecondaryParticles(mode='linear', pts_sec=ppSnd_s$ID$, v_sec=pVelSnd_pp$ID$, l_sec=pLifeSnd_pp$ID$, f_sec=pForceSnd_pp$ID$, flags=flags_s$ID$, v=vel_s$ID$, neighborRatio=neighborRatio_s$ID$, radius=1, gravity=gravity_s$ID$, k_b=k_b_s$ID$, k_d=k_d_s$ID$, c_s=c_s_s$ID$, c_b=c_b_s$ID$, dt=s$ID$.frameLength)\n\
        if $SNDPARTICLE_BOUNDARY_PUSHOUT$:\n\
            pushOutofObs(parts=ppSnd_s$ID$, flags=flags_s$ID$, phiObs=phiObs_s$ID$, shift=1.0)\n\
        flipDeleteParticlesInObstacle(pts=ppSnd_s$ID$, flags=flags_s$ID$)\n\
        #debugGridInfo(flags=flags_s$ID$, grid=trappedAir_s$ID$, name='Trapped Air')\n\
        #debugGridInfo(flags=flags_s$ID$, grid=waveCrest_s$ID$, name='Wave Crest')\n\
        #debugGridInfo(flags=flags_s$ID$, grid=kineticEnergy_s$ID$, name='Kinetic Energy')\n";

const std::string liquid_step_particles_high = "\n\
def liquid_step_particles_high_$ID$():\n\
    mantaMsg('Secondary particles step high (WORK IN PROGRESS)')\n\
    mantaMsg('      Grid size: ' + str(xl$ID$.getGridSize()))\n\
    if using_sndparts_s$ID$:\n\
        setFlagsFromLevelset(flags=flags_xl$ID$, phi=phi_xl$ID$)\n\
        #flags_xl$ID$.printGrid(5, false)\n\
        pVel_pp$ID$.multConst(Vec3(2,2,2))\n\
        mapPartsToMAC(vel=vel_xl$ID$, flags=flags_xl$ID$, velOld=vel_xl$ID$, parts=pp_s$ID$, partVel=pVel_pp$ID$)\n\
        radius = 3 if $SNDPARTICLE_POTENTIAL_QUALITY_HIGH$ else 2\n\
        flipComputeSecondaryParticlePotentials(potTA=trappedAir_xl$ID$, potWC=waveCrest_xl$ID$, potKE=kineticEnergy_xl$ID$, neighborRatio=neighborRatio_xl$ID$, flags=flags_xl$ID$, v=vel_xl$ID$, normal=normal_xl$ID$, phi=phi_xl$ID$, radius=radius, tauMinTA=tauMin_ta_s$ID$, tauMaxTA=tauMax_ta_s$ID$, tauMinWC=tauMin_wc_s$ID$, tauMaxWC=tauMax_wc_s$ID$, tauMinKE=tauMin_k_s$ID$, tauMaxKE=tauMax_k_s$ID$, scaleFromManta=scaleFromManta_s$ID$)\n\
        flipSampleSecondaryParticles(mode='single', flags=flags_xl$ID$, v=vel_xl$ID$, pts_sec=ppSnd_s$ID$, v_sec=pVelSnd_pp$ID$, l_sec=pLifeSnd_pp$ID$, lMin=lMin_s$ID$, lMax=lMax_s$ID$, potTA=trappedAir_xl$ID$, potWC=waveCrest_xl$ID$, potKE=kineticEnergy_xl$ID$, neighborRatio=neighborRatio_xl$ID$, c_s=c_s_s$ID$, c_b=c_b_s$ID$, k_ta=k_ta_s$ID$, k_wc=k_wc_s$ID$, dt=s$ID$.frameLength)\n\
        flipUpdateSecondaryParticles(mode='linear', pts_sec=ppSnd_s$ID$, v_sec=pVelSnd_pp$ID$, l_sec=pLifeSnd_pp$ID$, f_sec=pForceSnd_pp$ID$, flags=flags_xl$ID$, v=vel_xl$ID$, neighborRatio=neighborRatio_xl$ID$, radius=1, gravity=gravity_s$ID$, k_b=k_b_s$ID$, k_d=k_d_s$ID$, c_s=c_s_s$ID$, c_b=c_b_s$ID$, dt=s$ID$.frameLength)\n\
        if $SNDPARTICLE_BOUNDARY_PUSHOUT$:\n\
            interpolateGrid(target=phiOut_xl$ID$, source=phiOutIn_s$ID$)\n\
            interpolateGrid(target=phiObs_xl$ID$, source=phiObs_s$ID$)\n\
            setObstacleFlags(flags=flags_xl$ID$, phiObs=phiObs_xl$ID$, phiOut=phiOut_xl$ID$)\n\
            pushOutofObs(parts=ppSnd_s$ID$, flags=flags_xl$ID$, phiObs=phiObs_xl$ID$, shift=1.0)\n\
        flipDeleteParticlesInObstacle(pts=ppSnd_s$ID$, flags=flags_s$ID$)\n\
        debugGridInfo(flags=flags_xl$ID$, grid=trappedAir_xl$ID$, name='Trapped Air')\n\
        debugGridInfo(flags=flags_xl$ID$, grid=waveCrest_xl$ID$, name='Wave Crest')\n\
        debugGridInfo(flags=flags_xl$ID$, grid=kineticEnergy_xl$ID$, name='Kinetic Energy')\n";

//////////////////////////////////////////////////////////////////////
// IMPORT
//////////////////////////////////////////////////////////////////////

const std::string liquid_load_geometry_low = "\n\
def liquid_load_geometry_low_$ID$(path, framenr):\n\
    mantaMsg('Liquid load geometry low')\n\
    framenr = fluid_cache_get_framenr_formatted_$ID$(framenr)\n\
    phiIn_s$ID$.load(os.path.join(path, 'phiIn_' + framenr + '.uni'))\n";

const std::string liquid_load_data_low = "\n\
def liquid_load_data_low_$ID$(path, framenr, withParticles):\n\
    mantaMsg('Liquid load data low')\n\
    framenr = fluid_cache_get_framenr_formatted_$ID$(framenr)\n\
    phi_s$ID$.load(os.path.join(path, 'phi_' + framenr +'.uni'))\n\
    phiIn_s$ID$.load(os.path.join(path, 'phiIn_' + framenr + '.uni'))\n\
    pressure_s$ID$.load(os.path.join(path, 'pressure_' + framenr + '.uni'))\n\
    if withParticles:\n\
        pp_s$ID$.load(os.path.join(path, 'pp_' + framenr +'.uni'))\n\
        pVel_pp$ID$.load(os.path.join(path, 'pVel_' + framenr +'.uni'))\n";

const std::string liquid_load_data_high = "\n\
def liquid_load_data_high_$ID$(path, framenr, withParticles):\n\
    mantaMsg('Liquid load data high')\n\
    framenr = fluid_cache_get_framenr_formatted_$ID$(framenr)\n\
    phi_xl$ID$.load(os.path.join(path, 'phi_xl_' + framenr +'.uni'))\n\
    if withParticles:\n\
        pp_xl$ID$.load(os.path.join(path, 'pp_xl_' + framenr +'.uni'))\n";

const std::string liquid_load_mesh_low = "\n\
def liquid_load_mesh_low_$ID$(path, framenr):\n\
    mantaMsg('Liquid load mesh low')\n\
    framenr = fluid_cache_get_framenr_formatted_$ID$(framenr)\n\
    mesh_s$ID$.load(os.path.join(path, 'mesh_low_' + framenr + '.bobj.gz'))\n";

const std::string liquid_load_mesh_high = "\n\
def liquid_load_mesh_high_$ID$(path, framenr):\n\
    mantaMsg('Liquid load mesh high')\n\
    framenr = fluid_cache_get_framenr_formatted_$ID$(framenr)\n\
    mesh_xl$ID$.load(os.path.join(path, 'mesh_high_' + framenr + '.bobj.gz'))\n";

const std::string liquid_load_particles_low = "\n\
def liquid_load_particles_low_$ID$(path, framenr):\n\
    mantaMsg('Liquid load particles low')\n\
    framenr = fluid_cache_get_framenr_formatted_$ID$(framenr)\n\
    ppSnd_s$ID$.load(os.path.join(path, 'ppSnd_' + framenr + '.uni'))\n\
    pVelSnd_pp$ID$.load(os.path.join(path, 'pVelSnd_' + framenr + '.uni'))\n\
    pLifeSnd_pp$ID$.load(os.path.join(path, 'pLifeSnd_' + framenr + '.uni'))\n\
    if os.path.isfile(os.path.join(path, 'sndTrappedAir_' + framenr + '.uni')):\n\
        trappedAir_s$ID$.load(os.path.join(path, 'sndTrappedAir_' + framenr + '.uni'))\n\
        waveCrest_s$ID$.load(os.path.join(path, 'sndWaveCrest_' + framenr + '.uni'))\n\
        kineticEnergy_s$ID$.load(os.path.join(path, 'sndKineticEnergy_' + framenr + '.uni'))\n";

const std::string liquid_load_particles_high = "\n\
def liquid_load_particles_high_$ID$(path, framenr):\n\
    mantaMsg('Liquid load particles high (WORK IN PROGRESS)')\n\
    framenr = fluid_cache_get_framenr_formatted_$ID$(framenr)\n\
    ppSnd_s$ID$.load(os.path.join(path, 'ppSnd_xl_' + framenr + '.uni'))\n\
    pVelSnd_pp$ID$.load(os.path.join(path, 'pVelSnd_xl_' + framenr + '.uni'))\n\
    pLifeSnd_pp$ID$.load(os.path.join(path, 'pLifeSnd_xl_' + framenr + '.uni'))\n\
    if os.path.isfile(os.path.join(path, 'sndTrappedAir_xl_' + framenr + '.uni')):\n\
        trappedAir_xl$ID$.load(os.path.join(path, 'sndTrappedAir_xl_' + framenr + '.uni'))\n\
        waveCrest_xl$ID$.load(os.path.join(path, 'sndWaveCrest_xl_' + framenr + '.uni'))\n\
        kineticEnergy_xl$ID$.load(os.path.join(path, 'sndKineticEnergy_xl_' + framenr + '.uni'))\n";


//////////////////////////////////////////////////////////////////////
// EXPORT
//////////////////////////////////////////////////////////////////////

const std::string liquid_save_geometry_low = "\n\
def liquid_save_geometry_low_$ID$(path, framenr):\n\
    mantaMsg('Liquid save geometry')\n\
    framenr = fluid_cache_get_framenr_formatted_$ID$(framenr)\n\
    phiIn_s$ID$.save(os.path.join(path, 'phiIn_' + framenr + '.uni'))\n";

const std::string liquid_save_data_low = "\n\
def liquid_save_data_low_$ID$(path, framenr):\n\
    mantaMsg('Liquid save data low')\n\
    framenr = fluid_cache_get_framenr_formatted_$ID$(framenr)\n\
    phi_s$ID$.save(os.path.join(path, 'phi_' + framenr +'.uni'))\n\
    phiIn_s$ID$.save(os.path.join(path, 'phiIn_' + framenr + '.uni'))\n\
    pp_s$ID$.save(os.path.join(path, 'pp_' + framenr +'.uni'))\n\
    pVel_pp$ID$.save(os.path.join(path, 'pVel_' + framenr +'.uni'))\n\
    pressure_s$ID$.save(os.path.join(path, 'pressure_' + framenr + '.uni'))\n";

const std::string liquid_save_data_high = "\n\
def liquid_save_data_high_$ID$(path, framenr):\n\
    mantaMsg('Liquid save data high')\n\
    framenr = fluid_cache_get_framenr_formatted_$ID$(framenr)\n\
    phi_xl$ID$.save(os.path.join(path, 'phi_xl_' + framenr +'.uni'))\n\
    pp_s$ID$.save(os.path.join(path, 'pp_xl_' + framenr +'.uni'))\n";

const std::string liquid_save_mesh_low = "\n\
def liquid_save_mesh_low_$ID$(path, framenr):\n\
    mantaMsg('Liquid save mesh low')\n\
    framenr = fluid_cache_get_framenr_formatted_$ID$(framenr)\n\
    phiParts_s$ID$.copyFrom(phi_s$ID$) # mis-use phiParts as temp grid to close the mesh\n\
    phiParts_s$ID$.setBound(0.5,0)\n\
    phiParts_s$ID$.createMesh(mesh_s$ID$)\n\
    mesh_s$ID$.save(os.path.join(path, 'mesh_low_' + framenr + '.bobj.gz'))\n";

const std::string liquid_save_mesh_high = "\n\
def liquid_save_mesh_high_$ID$(path, framenr):\n\
    mantaMsg('Liquid save mesh high')\n\
    framenr = fluid_cache_get_framenr_formatted_$ID$(framenr)\n\
    phiParts_s$ID$.copyFrom(phi_s$ID$) # mis-use phiParts as temp grid to close the mesh\n\
    phiParts_s$ID$.setBound(0.5,0)\n\
    interpolateGrid(target=phi_xl$ID$, source=phiParts_s$ID$)\n\
    phi_xl$ID$.join(phiParts_xl$ID$)\n\
    phi_xl$ID$.createMesh(mesh_xl$ID$)\n\
    mesh_xl$ID$.save(os.path.join(path, 'mesh_high_' + framenr + '.bobj.gz'))\n";

const std::string liquid_save_particles_low = "\n\
def liquid_save_particles_low_$ID$(path, framenr):\n\
    mantaMsg('Liquid save particles low')\n\
    framenr = fluid_cache_get_framenr_formatted_$ID$(framenr)\n\
    ppSnd_s$ID$.save(os.path.join(path, 'ppSnd_' + framenr + '.uni'))\n\
    pVelSnd_pp$ID$.save(os.path.join(path, 'pVelSnd_' + framenr + '.uni'))\n\
    pLifeSnd_pp$ID$.save(os.path.join(path, 'pLifeSnd_' + framenr + '.uni'))\n\
    if $SNDPARTICLE_POTENTIAL_GRID_SAVE_LOW$ or $SNDPARTICLE_POTENTIAL_GRID_SAVE_ON$:\n\
        trappedAir_s$ID$.save(os.path.join(path, 'sndTrappedAir_' + framenr + '.uni'))\n\
        waveCrest_s$ID$.save(os.path.join(path, 'sndWaveCrest_' + framenr + '.uni'))\n\
        kineticEnergy_s$ID$.save(os.path.join(path, 'sndKineticEnergy_' + framenr + '.uni'))\n";

const std::string liquid_save_particles_high = "\n\
def liquid_save_particles_high_$ID$(path, framenr):\n\
    mantaMsg('Liquid save particles high (WORK IN PROGRESS)')\n\
    framenr = fluid_cache_get_framenr_formatted_$ID$(framenr)\n\
    ppSnd_s$ID$.save(os.path.join(path, 'ppSnd_xl_' + framenr + '.uni'))\n\
    pVelSnd_pp$ID$.save(os.path.join(path, 'pVelSnd_xl_' + framenr + '.uni'))\n\
    pLifeSnd_pp$ID$.save(os.path.join(path, 'pLifeSnd_xl_' + framenr + '.uni'))\n\
    if $SNDPARTICLE_POTENTIAL_GRID_SAVE_ON$:\n\
        trappedAir_xl$ID$.save(os.path.join(path, 'sndTrappedAir_xl_' + framenr + '.uni'))\n\
        waveCrest_xl$ID$.save(os.path.join(path, 'sndWaveCrest_xl_' + framenr + '.uni'))\n\
        kineticEnergy_xl$ID$.save(os.path.join(path, 'sndKineticEnergy_xl_' + framenr + '.uni'))\n";

//////////////////////////////////////////////////////////////////////
// STANDALONE MODE
//////////////////////////////////////////////////////////////////////

const std::string liquid_standalone_load = "\n\
# import *.uni files\n\
path_prefix_$ID$ = '$MANTA_EXPORT_PATH$'\n\
load_liquid_data_low_$ID$(path_prefix_$ID$)\n\
if using_highres_s$ID$:\n\
    load_liquid_data_high_$ID$(path_prefix_$ID$)\n";
