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
mantaMsg('Liquid variables')\n\
narrowBandWidth_s$ID$         = 3\n\
combineBandWidth_s$ID$        = narrowBandWidth_s$ID$ - 1\n\
adjustedNarrowBandWidth_s$ID$ = $PARTICLE_BAND_WIDTH$ # only used in adjustNumber to control band width\n\
particleNumber_s$ID$   = $PARTICLE_NUMBER$\n\
minParticles_s$ID$     = $PARTICLE_MINIMUM$\n\
maxParticles_s$ID$     = $PARTICLE_MAXIMUM$\n\
radiusFactor_s$ID$     = $PARTICLE_RADIUS$\n\
using_mesh_s$ID$       = $USING_MESH$\n\
using_final_mesh_s$ID$ = $USING_IMPROVED_MESH$\n\
smoothenUpper_s$ID$    = $MESH_SMOOTHEN_UPPER$\n\
smoothenLower_s$ID$    = $MESH_SMOOTHEN_LOWER$\n\
smoothenPos_s$ID$      = $MESH_SMOOTHEN_POS$\n\
smoothenNeg_s$ID$      = $MESH_SMOOTHEN_NEG$\n\
randomness_s$ID$       = $PARTICLE_RANDOMNESS$\n\
surfaceTension_s$ID$   = $LIQUID_SURFACE_TENSION$\n";

//////////////////////////////////////////////////////////////////////
// GRIDS & MESH & PARTICLESYSTEM
//////////////////////////////////////////////////////////////////////

const std::string liquid_alloc = "\n\
mantaMsg('Liquid alloc')\n\
phiParts_s$ID$   = s$ID$.create(LevelsetGrid)\n\
phi_s$ID$        = s$ID$.create(LevelsetGrid)\n\
phiTmp_s$ID$     = s$ID$.create(LevelsetGrid)\n\
phiIn_s$ID$      = s$ID$.create(LevelsetGrid)\n\
curvature_s$ID$  = s$ID$.create(RealGrid)\n\
fractions_s$ID$  = 0# s$ID$.create(MACGrid) # TODO (sebbas): disabling fractions for now - not fracwallbcs not supporting obvels yet\n\
velOld_s$ID$     = s$ID$.create(MACGrid)\n\
velParts_s$ID$   = s$ID$.create(MACGrid)\n\
mapWeights_s$ID$ = s$ID$.create(MACGrid)\n\
\n\
pp_s$ID$         = s$ID$.create(BasicParticleSystem)\n\
pVel_pp$ID$      = pp_s$ID$.create(PdataVec3)\n\
\n\
# Acceleration data for particle nbs\n\
pindex_s$ID$     = s$ID$.create(ParticleIndexSystem)\n\
gpi_s$ID$        = s$ID$.create(IntGrid)\n\
\n\
# Keep track of important objects in dict to load them later on\n\
liquid_data_dict_s$ID$ = dict(phiParts=phiParts_s$ID$, phi=phi_s$ID$, phiIn=phiIn_s$ID$, phiTmp=phiTmp_s$ID$)\n\
liquid_flip_dict_s$ID$ = dict(pp=pp_s$ID$, pVel=pVel_pp$ID$)\n";

const std::string liquid_alloc_mesh = "\n\
mantaMsg('Liquid alloc mesh')\n\
phiParts_sm$ID$ = sm$ID$.create(LevelsetGrid)\n\
phi_sm$ID$      = sm$ID$.create(LevelsetGrid)\n\
pp_sm$ID$       = sm$ID$.create(BasicParticleSystem)\n\
flags_sm$ID$    = sm$ID$.create(FlagGrid)\n\
mesh_sm$ID$     = sm$ID$.create(Mesh)\n\
\n\
if using_speedvectors_s$ID$:\n\
    mVel_mesh$ID$ = mesh_sm$ID$.create(MdataVec3)\n\
    vel_sm$ID$    = sm$ID$.create(MACGrid)\n\
\n\
# Acceleration data for particle nbs\n\
pindex_sm$ID$  = sm$ID$.create(ParticleIndexSystem)\n\
gpi_sm$ID$     = sm$ID$.create(IntGrid)\n\
\n\
# Keep track of important objects in dict to load them later on\n\
liquid_mesh_dict_s$ID$ = dict(lMesh=mesh_sm$ID$)\n\
\n\
if using_speedvectors_s$ID$:\n\
    liquid_meshvel_dict_s$ID$ = dict(lVelMesh=mVel_mesh$ID$)\n";

const std::string liquid_init_phi = "\n\
# Prepare domain\n\
phi_s$ID$.initFromFlags(flags_s$ID$)\n\
phiIn_s$ID$.initFromFlags(flags_s$ID$)\n";

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
    fluid_pre_step_$ID$()\n\
    \n\
    flags_s$ID$.initDomain(boundaryWidth=0, phiWalls=phiObs_s$ID$, outflow=boundConditions_s$ID$)\n\
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
    if using_invel_s$ID$:\n\
        resampleVec3ToMac(source=invelC_s$ID$, target=invel_s$ID$)\n\
        pVel_pp$ID$.setSource(invel_s$ID$, isMAC=True)\n\
    \n\
    sampleLevelsetWithParticles(phi=phiIn_s$ID$, flags=flags_s$ID$, parts=pp_s$ID$, discretization=particleNumber_s$ID$, randomness=randomness_s$ID$, refillEmpty=True)\n\
    flags_s$ID$.updateFromLevelset(phi_s$ID$, phiObs_s$ID$)\n\
    \n\
    mantaMsg('Liquid step / s$ID$.frame: ' + str(s$ID$.frame))\n\
    liquid_step_$ID$()\n\
    \n\
    s$ID$.step()\n\
    \n\
    fluid_post_step_$ID$()\n";

const std::string liquid_step = "\n\
def liquid_step_$ID$():\n\
    mantaMsg('Liquid step')\n\
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
    phi_s$ID$.setBoundNeumann(0) # make sure no particles are placed at outer boundary\n\
    \n\
    if doOpen_s$ID$:\n\
        resetOutflow(flags=flags_s$ID$, phi=phi_s$ID$, parts=pp_s$ID$, index=gpi_s$ID$, indexSys=pindex_s$ID$)\n\
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
    \n\
    mantaMsg('Adding external forces')\n\
    addForceField(flags=flags_s$ID$, vel=vel_s$ID$, force=forces_s$ID$)\n\
    \n\
    if using_obstacle_s$ID$:\n\
        mantaMsg('Extrapolating object velocity')\n\
        # ensure velocities inside of obs object, slightly add obvels outside of obs object\n\
        extrapolateVec3Simple(vel=obvelC_s$ID$, phi=phiObsIn_s$ID$, distance=int(res_s$ID$/2), inside=True)\n\
        extrapolateVec3Simple(vel=obvelC_s$ID$, phi=phiObsIn_s$ID$, distance=3, inside=False)\n\
        resampleVec3ToMac(source=obvelC_s$ID$, target=obvel_s$ID$)\n\
    \n\
    extrapolateMACSimple(flags=flags_s$ID$, vel=vel_s$ID$, distance=2)#, intoObs=True) # TODO (sebbas): uncomment for fraction support\n\
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
        PD_fluid_guiding(vel=vel_s$ID$, velT=velT_s$ID$, flags=flags_s$ID$, phi=phi_s$ID$, curv=curvature_s$ID$, surfTens=surfaceTension_s$ID$, fractions=fractions_s$ID$, weight=weightGuide_s$ID$, blurRadius=beta_sg$ID$, pressure=pressure_s$ID$, tau=tau_sg$ID$, sigma=sigma_sg$ID$, theta=theta_sg$ID$, zeroPressureFixing=not doOpen_s$ID$)\n\
    else:\n\
        mantaMsg('Pressure')\n\
        solvePressure(flags=flags_s$ID$, vel=vel_s$ID$, pressure=pressure_s$ID$, phi=phi_s$ID$, curv=curvature_s$ID$, surfTens=surfaceTension_s$ID$, fractions=fractions_s$ID$)\n\
    \n\
    extrapolateMACSimple(flags=flags_s$ID$, vel=vel_s$ID$, distance=4)#, intoObs=True) # TODO (sebbas): uncomment for fraction support\n\
    setWallBcs(flags=flags_s$ID$, vel=vel_s$ID$, obvel=obvel_s$ID$ if using_obstacle_s$ID$ else 0, phiObs=phiObs_s$ID$, fractions=fractions_s$ID$)\n\
    \n\
    extrapolateMACSimple(flags=flags_s$ID$, vel=vel_s$ID$, distance=(int(maxVel_s$ID$*1.25 )) ) # TODO (sebbas): extrapolation because of no fractions\n\
    # set source grids for resampling, used in adjustNumber!\n\
    pVel_pp$ID$.setSource(vel_s$ID$, isMAC=True)\n\
    adjustNumber(parts=pp_s$ID$, vel=vel_s$ID$, flags=flags_s$ID$, minParticles=minParticles_s$ID$, maxParticles=maxParticles_s$ID$, phi=phi_s$ID$, exclude=phiObs_s$ID$, radiusFactor=1., narrowBand=adjustedNarrowBandWidth_s$ID$)\n\
    flipVelocityUpdate(vel=vel_s$ID$, velOld=velOld_s$ID$, flags=flags_s$ID$, parts=pp_s$ID$, partVel=pVel_pp$ID$, flipRatio=0.97)\n\
    \n\
    sampleLevelsetWithParticles(phi=phiIn_s$ID$, flags=flags_s$ID$, parts=pp_s$ID$, discretization=particleNumber_s$ID$, randomness=randomness_s$ID$, refillEmpty=True)\n";

const std::string liquid_step_mesh = "\n\
def liquid_step_mesh_$ID$():\n\
    mantaMsg('Liquid step mesh')\n\
    \n\
    interpolateGrid(target=phi_sm$ID$, source=phiTmp_s$ID$) # mis-use phiParts as temp grid\n\
    \n\
    # create surface\n\
    pp_sm$ID$.readParticles(pp_s$ID$)\n\
    gridParticleIndex(parts=pp_sm$ID$, flags=flags_sm$ID$, indexSys=pindex_sm$ID$, index=gpi_sm$ID$)\n\
    \n\
    if using_final_mesh_s$ID$:\n\
        mantaMsg('Liquid using improved particle levelset')\n\
        improvedParticleLevelset(pp_sm$ID$, pindex_sm$ID$, flags_sm$ID$, gpi_sm$ID$, phiParts_sm$ID$, radiusFactor_s$ID$, smoothenPos_s$ID$, smoothenNeg_s$ID$, smoothenLower_s$ID$, smoothenUpper_s$ID$)\n\
    else:\n\
        mantaMsg('Liquid using union particle levelset')\n\
        unionParticleLevelset(pp_sm$ID$, pindex_sm$ID$, flags_sm$ID$, gpi_sm$ID$, phiParts_sm$ID$, radiusFactor_s$ID$)\n\
    \n\
    phi_sm$ID$.addConst(1.) # shrink slightly\n\
    phi_sm$ID$.join(phiParts_sm$ID$)\n\
    extrapolateLsSimple(phi=phi_sm$ID$, distance=narrowBandWidth_s$ID$+2, inside=True)\n\
    extrapolateLsSimple(phi=phi_sm$ID$, distance=3)\n\
    phi_sm$ID$.setBoundNeumann(boundaryWidth_s$ID$) # make sure no particles are placed at outer boundary\n\
    \n\
    # Vert vel vector needs to pull data from vel grid with correct dim\n\
    if using_speedvectors_s$ID$:\n\
        interpolateMACGrid(target=vel_sm$ID$, source=vel_s$ID$)\n\
        mVel_mesh$ID$.setSource(vel_sm$ID$, isMAC=True)\n\
    \n\
    phi_sm$ID$.setBound(0.5,int(((upres_sm$ID$)*2)-2) )\n\
    phi_sm$ID$.createMesh(mesh_sm$ID$)\n";

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
    if doOpen_s$ID$:\n\
        resetOutflow(flags=flags_s$ID$, parts=ppSnd_sp$ID$)\n\
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
def liquid_load_data_$ID$(path, framenr, file_format):\n\
    mantaMsg('Liquid load data')\n\
    fluid_file_import_s$ID$(dict=liquid_data_dict_s$ID$, path=path, framenr=framenr, file_format=file_format)\n";

const std::string liquid_load_flip = "\n\
def liquid_load_flip_$ID$(path, framenr, file_format):\n\
    mantaMsg('Liquid load flip')\n\
    fluid_file_import_s$ID$(dict=liquid_flip_dict_s$ID$, path=path, framenr=framenr, file_format=file_format)\n";

const std::string liquid_load_mesh = "\n\
def liquid_load_mesh_$ID$(path, framenr, file_format):\n\
    mantaMsg('Liquid load mesh')\n\
    fluid_file_import_s$ID$(dict=liquid_mesh_dict_s$ID$, path=path, framenr=framenr, file_format=file_format)\n";

const std::string liquid_load_meshvel = "\n\
def liquid_load_meshvel_$ID$(path, framenr, file_format):\n\
    mantaMsg('Liquid load meshvel')\n\
    fluid_file_import_s$ID$(dict=liquid_meshvel_dict_s$ID$, path=path, framenr=framenr, file_format=file_format)\n";

//////////////////////////////////////////////////////////////////////
// EXPORT
//////////////////////////////////////////////////////////////////////

const std::string liquid_save_data = "\n\
def liquid_save_data_$ID$(path, framenr, file_format):\n\
    mantaMsg('Liquid save data')\n\
    fluid_file_export_s$ID$(dict=liquid_data_dict_s$ID$, path=path, framenr=framenr, file_format=file_format)\n";

const std::string liquid_save_flip = "\n\
def liquid_save_flip_$ID$(path, framenr, file_format):\n\
    mantaMsg('Liquid save flip')\n\
    fluid_file_export_s$ID$(dict=liquid_flip_dict_s$ID$, path=path, framenr=framenr, file_format=file_format)\n";

const std::string liquid_save_mesh = "\n\
def liquid_save_mesh_$ID$(path, framenr, file_format):\n\
    mantaMsg('Liquid save mesh')\n\
    fluid_file_export_s$ID$(dict=liquid_mesh_dict_s$ID$, path=path, framenr=framenr, file_format=file_format)\n";

const std::string liquid_save_meshvel = "\n\
def liquid_save_meshvel_$ID$(path, framenr, file_format):\n\
    mantaMsg('Liquid save mesh vel')\n\
    fluid_file_export_s$ID$(dict=liquid_meshvel_dict_s$ID$, path=path, framenr=framenr, file_format=file_format)\n";

//////////////////////////////////////////////////////////////////////
// STANDALONE MODE
//////////////////////////////////////////////////////////////////////

const std::string liquid_standalone = "\n\
# Helper function to call cache load functions\n\
def load(frame):\n\
    fluid_load_data_$ID$(os.path.join(cache_dir, 'data'), frame, file_format_data)\n\
    liquid_load_data_$ID$(os.path.join(cache_dir, 'data'), frame, file_format_data)\n\
    liquid_load_flip_$ID$(os.path.join(cache_dir, 'data'), frame, file_format_particles)\n\
    if using_sndparts_s$ID$:\n\
        fluid_load_particles_$ID$(os.path.join(cache_dir, 'particles'), frame, file_format_particles)\n\
    if using_mesh_s$ID$:\n\
        liquid_load_mesh_$ID$(os.path.join(cache_dir, 'mesh'), frame, file_format_mesh)\n\
    if using_guiding_s$ID$:\n\
        fluid_load_guiding_$ID$(os.path.join(cache_dir, 'guiding'), frame, file_format_data)\n\
\n\
# Helper function to call step functions\n\
def step(frame):\n\
    liquid_adaptive_step_$ID$(frame)\n\
    if using_mesh_s$ID$:\n\
        liquid_step_mesh_$ID$()\n\
    if using_sndparts_s$ID$:\n\
        liquid_step_particles_$ID$()\n";



