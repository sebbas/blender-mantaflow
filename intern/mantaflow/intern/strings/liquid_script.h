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
// BOUNDS
//////////////////////////////////////////////////////////////////////

const std::string liquid_bounds_low = "\n\
# prepare domain low\n\
mantaMsg('Liquid domain low')\n\
flags_s$ID$.initDomain(boundaryWidth=boundaryWidth_s$ID$, phiWalls=phiObs_s$ID$)\n\
if doOpen_s$ID$:\n\
    setOpenBound(flags=flags_s$ID$, bWidth=boundaryWidth_s$ID$, openBound=boundConditions_s$ID$, type=FlagOutflow|FlagEmpty)\n";

const std::string liquid_bounds_high = "\n\
# prepare domain high\n\
mantaMsg('Liquid domain high')\n\
flags_xl$ID$.initDomain(boundaryWidth=boundaryWidth_s$ID$)\n\
if doOpen_s$ID$:\n\
    setOpenBound(flags=flags_xl$ID$, bWidth=boundaryWidth_s$ID$, openBound=boundConditions_s$ID$, type=FlagOutflow|FlagEmpty)\n";

//////////////////////////////////////////////////////////////////////
// VARIABLES
//////////////////////////////////////////////////////////////////////

const std::string liquid_variables_low = "\n\
mantaMsg('Liquid variables low')\n\
narrowBandWidth_s$ID$  = 3\n\
combineBandWidth_s$ID$ = narrowBandWidth_s$ID$ - 1\n\
\n\
particleNumber_s$ID$ = $PARTICLE_NUMBER$\n\
minParticles_s$ID$   = pow(particleNumber_s$ID$, dim_s$ID$)\n\
radiusFactor_s$ID$   = $PARTICLE_RADIUS$\n\
randomness_s$ID$     = $PARTICLE_RANDOMNESS$\n";

const std::string liquid_variables_high = "\n\
mantaMsg('Liquid variables high')\n";

//////////////////////////////////////////////////////////////////////
// GRIDS & MESH & PARTICLESYSTEM
//////////////////////////////////////////////////////////////////////

const std::string liquid_alloc_low = "\n\
mantaMsg('Liquid alloc low')\n\
flags_s$ID$      = s$ID$.create(FlagGrid)\n\
numObs_s$ID$     = s$ID$.create(IntGrid)\n\
phiParts_s$ID$   = s$ID$.create(LevelsetGrid)\n\
phi_s$ID$        = s$ID$.create(LevelsetGrid)\n\
phiIn_s$ID$      = s$ID$.create(LevelsetGrid)\n\
phiOut_s$ID$     = s$ID$.create(LevelsetGrid)\n\
pressure_s$ID$   = s$ID$.create(RealGrid)\n\
\n\
phiObs_s$ID$     = s$ID$.create(LevelsetGrid)\n\
phiObsIn_s$ID$   = s$ID$.create(LevelsetGrid)\n\
fractions_s$ID$  = s$ID$.create(MACGrid)\n\
\n\
vel_s$ID$        = s$ID$.create(MACGrid)\n\
x_vel_s$ID$      = s$ID$.create(RealGrid)\n\
y_vel_s$ID$      = s$ID$.create(RealGrid)\n\
z_vel_s$ID$      = s$ID$.create(RealGrid)\n\
obvel_s$ID$      = s$ID$.create(MACGrid)\n\
x_obvel_s$ID$    = s$ID$.create(RealGrid)\n\
y_obvel_s$ID$    = s$ID$.create(RealGrid)\n\
z_obvel_s$ID$    = s$ID$.create(RealGrid)\n\
velOld_s$ID$     = s$ID$.create(MACGrid)\n\
velParts_s$ID$   = s$ID$.create(MACGrid)\n\
mapWeights_s$ID$ = s$ID$.create(MACGrid)\n\
\n\
pp_s$ID$         = s$ID$.create(BasicParticleSystem)\n\
pVel_s$ID$       = pp_s$ID$.create(PdataVec3)\n\
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
    copyRealToVec3(sourceX=x_vel_s$ID$, sourceY=y_vel_s$ID$, sourceZ=z_vel_s$ID$, target=vel_s$ID$)\n\
    copyRealToVec3(sourceX=x_obvel_s$ID$, sourceY=y_obvel_s$ID$, sourceZ=z_obvel_s$ID$, target=obvel_s$ID$)\n\
    copyRealToVec3(sourceX=x_force_s$ID$, sourceY=y_force_s$ID$, sourceZ=z_force_s$ID$, target=forces_s$ID$)\n\
    \n\
    clearInObstacle(flags=flags_s$ID$, grid=phi_s$ID$)\n";

const std::string liquid_pre_step_high = "\n\
def liquid_pre_step_high_s$ID$():\n\
    clearInObstacle(flags=flags_xl$ID$, grid=phi_xl$ID$)\n";

const std::string liquid_post_step_low = "\n\
def liquid_post_step_low_$ID$():\n\
    forces_s$ID$.clear()\n\
    obvel_s$ID$.clear()\n\
    \n\
    phiIn_s$ID$.setConst(0.5)\n\
    phiObs_s$ID$.setConst(0.5)\n\
    phiObsIn_s$ID$.setConst(0)\n\
    \n\
    copyVec3ToReal(source=vel_s$ID$, targetX=x_vel_s$ID$, targetY=y_vel_s$ID$, targetZ=z_vel_s$ID$)\n";

const std::string liquid_post_step_high = "\n\
def liquid_post_step_high_$ID$():\n\
    clearInObstacle(flags=flags_xl$ID$, grid=phi_xl$ID$)\n";

//////////////////////////////////////////////////////////////////////
// STEP FUNCTIONS
//////////////////////////////////////////////////////////////////////

const std::string liquid_adaptive_step = "\n\
def manta_step_$ID$(start_frame):\n\
    s$ID$.frame = start_frame\n\
    s$ID$.timeTotal = s$ID$.frame * dt0_s$ID$\n\
    last_frame_s$ID$ = s$ID$.frame\n\
    \n\
    liquid_pre_step_low_$ID$()\n\
    if using_highres_s$ID$:\n\
        liquid_pre_step_high_s$ID$()\n\
    \n\
    while s$ID$.frame == last_frame_s$ID$:\n\
        \n\
        flags_s$ID$.initDomain(boundaryWidth=boundaryWidth_s$ID$, phiWalls=phiObs_s$ID$)\n\
        if doOpen_s$ID$:\n\
            setOpenBound(flags=flags_s$ID$, bWidth=boundaryWidth_s$ID$, openBound=boundConditions_s$ID$, type=FlagOutflow|FlagEmpty)\n\
        \n\
        phiObs_s$ID$.join(phiObsIn_s$ID$)\n\
        phiIn_s$ID$.subtract(phiObs_s$ID$)\n\
        phi_s$ID$.join(phiIn_s$ID$)\n\
        \n\
        updateFractions(flags=flags_s$ID$, phiObs=phiObs_s$ID$, fractions=fractions_s$ID$, boundaryWidth=boundaryWidth_s$ID$)\n\
        setObstacleFlags(flags=flags_s$ID$, phiObs=phiObs_s$ID$, fractions=fractions_s$ID$, phiOut=phiOut_s$ID$)\n\
        \n\
        sampleLevelsetWithParticles(phi=phiIn_s$ID$, flags=flags_s$ID$, parts=pp_s$ID$, discretization=particleNumber_s$ID$, randomness=randomness_s$ID$, refillEmpty=True)\n\
        flags_s$ID$.updateFromLevelset(phi_s$ID$)\n\
        \n\
        if using_adaptTime_s$ID$:\n\
            mantaMsg('Adapt timestep')\n\
            maxvel_s$ID$ = vel_s$ID$.getMaxValue()\n\
            s$ID$.adaptTimestep(maxvel_s$ID$)\n\
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
    liquid_post_step_low_$ID$()\n\
    if using_highres_s$ID$:\n\
        liquid_post_step_high_$ID$()\n";

const std::string liquid_step_low = "\n\
def liquid_step_$ID$():\n\
    mantaMsg('Liquid step low')\n\
    # FLIP\n\
    # Create interpolated version of original phi grid for later use in (optional) high-res step\n\
    if using_highres_s$ID$:\n\
        interpolateGrid(target=phi_xl$ID$, source=phi_s$ID$)\n\
    \n\
    pp_s$ID$.advectInGrid(flags=flags_s$ID$, vel=vel_s$ID$, integrationMode=IntRK4, deleteInObstacle=False, stopInObstacle=False)\n\
    pushOutofObs(parts=pp_s$ID$, flags=flags_s$ID$, phiObs=phiObs_s$ID$)\n\
    \n\
    advectSemiLagrange(flags=flags_s$ID$, vel=vel_s$ID$, grid=phi_s$ID$, order=1, openBounds=doOpen_s$ID$, boundaryWidth=boundaryWidth_s$ID$) # first order is usually enough\n\
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
    #if doOpen:\n\
    resetOutflow(flags=flags_s$ID$, phi=phi_s$ID$, parts=pp_s$ID$, index=gpi_s$ID$, indexSys=pindex_s$ID$) # open boundaries\n\
    flags_s$ID$.updateFromLevelset(phi_s$ID$)\n\
    \n\
    # combine particles velocities with advected grid velocities\n\
    mapPartsToMAC(vel=velParts_s$ID$, flags=flags_s$ID$, velOld=velOld_s$ID$, parts=pp_s$ID$, partVel=pVel_s$ID$, weight=mapWeights_s$ID$)\n\
    extrapolateMACFromWeight(vel=velParts_s$ID$, distance=2, weight=mapWeights_s$ID$)\n\
    combineGridVel(vel=velParts_s$ID$, weight=mapWeights_s$ID$, combineVel=vel_s$ID$, phi=phi_s$ID$, narrowBand=combineBandWidth_s$ID$, thresh=0)\n\
    velOld_s$ID$.copyFrom(vel_s$ID$)\n\
    \n\
    # forces & pressure solve\n\
    addGravity(flags=flags_s$ID$, vel=vel_s$ID$, gravity=gravity_s$ID$)\n\
    addForceField(flags=flags_s$ID$, vel=vel_s$ID$, force=forces_s$ID$)\n\
    \n\
    extrapolateMACSimple(flags=flags_s$ID$, vel=vel_s$ID$, distance=2, intoObs=True)\n\
    setWallBcs(flags=flags_s$ID$, vel=vel_s$ID$, fractions=fractions_s$ID$, phiObs=phiObs_s$ID$)\n\
    \n\
    solvePressure(flags=flags_s$ID$, vel=vel_s$ID$, pressure=pressure_s$ID$, phi=phi_s$ID$, fractions=fractions_s$ID$)\n\
    \n\
    extrapolateMACSimple(flags=flags_s$ID$, vel=vel_s$ID$, distance=4, intoObs=True)\n\
    setWallBcs(flags=flags_s$ID$, vel=vel_s$ID$, fractions=fractions_s$ID$, phiObs=phiObs_s$ID$)\n\
    \n\
    clearInObstacle(flags=flags_s$ID$, grid=phi_s$ID$)\n\
    clearInObstacle(flags=flags_s$ID$, grid=phiParts_s$ID$)\n\
    pushOutofObs(parts=pp_s$ID$, flags=flags_s$ID$, phiObs=phiObs_s$ID$)\n\
    \n\
    if (dim_s$ID$==3):\n\
        # mis-use phiParts as temp grid to close the mesh\n\
        phiParts_s$ID$.copyFrom(phi_s$ID$)\n\
        phiParts_s$ID$.setBound(0.5,0)\n\
        phiParts_s$ID$.createMesh(mesh_s$ID$)\n\
    \n\
    # set source grids for resampling, used in adjustNumber!\n\
    pVel_s$ID$.setSource(vel_s$ID$, isMAC=True)\n\
    adjustNumber(parts=pp_s$ID$, vel=vel_s$ID$, flags=flags_s$ID$, minParticles=1*minParticles_s$ID$, maxParticles=2*minParticles_s$ID$, phi=phi_s$ID$, exclude=phiObs_s$ID$, radiusFactor=radiusFactor_s$ID$, narrowBand=narrowBandWidth_s$ID$)\n\
    flipVelocityUpdate(vel=vel_s$ID$, velOld=velOld_s$ID$, flags=flags_s$ID$, parts=pp_s$ID$, partVel=pVel_s$ID$, flipRatio=0.97)\n";

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

const std::string liquid_import_low = "\n\
def load_liquid_data_low_$ID$(path):\n\
    flags_s$ID$.load(os.path.join(path, 'flags_s$ID$.uni'))\n\
    \n\
    phiParts_s$ID$.load(os.path.join(path, 'phiParts_s$ID$.uni'))\n\
    phi_s$ID$.load(os.path.join(path, 'phi_s$ID$.uni'))\n\
    phiIn_s$ID$.load(os.path.join(path, 'phiIn_s$ID$.uni'))\n\
    phiObs_s$ID$.load(os.path.join(path, 'phiObs_s$ID$.uni'))\n\
    phiObsIn_s$ID$.load(os.path.join(path, 'phiObsIn_s$ID$.uni'))\n\
    phiOut_s$ID$.load(os.path.join(path, 'phiOut_s$ID$.uni'))\n\
    fractions_s$ID$.load(os.path.join(path, 'fractions_s$ID$.uni'))\n\
    pressure_s$ID$.load(os.path.join(path, 'pressure_s$ID$.uni'))\n\
    \n\
    vel_s$ID$.load(os.path.join(path, 'vel_s$ID$.uni'))\n\
    velOld_s$ID$.load(os.path.join(path, 'velOld_s$ID$.uni'))\n\
    velParts_s$ID$.load(os.path.join(path, 'velParts_s$ID$.uni'))\n\
    mapWeights_s$ID$.load(os.path.join(path, 'mapWeights_s$ID$.uni'))\n\
    \n\
    x_vel_s$ID$.load(os.path.join(path, 'x_vel_s$ID$.uni'))\n\
    y_vel_s$ID$.load(os.path.join(path, 'y_vel_s$ID$.uni'))\n\
    z_vel_s$ID$.load(os.path.join(path, 'z_vel_s$ID$.uni'))\n\
    x_obvel_s$ID$.load(os.path.join(path, 'x_obvel_s$ID$.uni'))\n\
    y_obvel_s$ID$.load(os.path.join(path, 'y_obvel_s$ID$.uni'))\n\
    z_obvel_s$ID$.load(os.path.join(path, 'z_obvel_s$ID$.uni'))\n\
    \n\
    pp_s$ID$.load(os.path.join(path, 'pp_s$ID$.uni'))\n\
    pVel_s$ID$.load(os.path.join(path, 'pVel_s$ID$.uni'))\n\
    \n\
    gpi_s$ID$.load(os.path.join(path, 'gpi_s$ID$.uni'))\n";

const std::string liquid_import_high = "\n\
def load_liquid_data_high_$ID$(path):\n\
    flags_xl$ID$.load(os.path.join(path, 'flags_xl$ID$.uni'))\n\
    \n\
    phiParts_xl$ID$.load(os.path.join(path, 'phiParts_xl$ID$.uni'))\n\
    phi_xl$ID$.load(os.path.join(path, 'phi_xl$ID$.uni'))\n\
    \n\
    pp_xl$ID$.load(os.path.join(path, 'pp_xl$ID$.uni'))\n";

const std::string liquid_export_low = "\n\
def save_liquid_data_low_$ID$(path):\n\
    flags_s$ID$.save(os.path.join(path, 'flags_s$ID$.uni'))\n\
    \n\
    phiParts_s$ID$.save(os.path.join(path, 'phiParts_s$ID$.uni'))\n\
    phi_s$ID$.save(os.path.join(path, 'phi_s$ID$.uni'))\n\
    phiIn_s$ID$.save(os.path.join(path, 'phiIn_s$ID$.uni'))\n\
    phiObs_s$ID$.save(os.path.join(path, 'phiObs_s$ID$.uni'))\n\
    phiObsIn_s$ID$.save(os.path.join(path, 'phiObsIn_s$ID$.uni'))\n\
    phiOut_s$ID$.save(os.path.join(path, 'phiOut_s$ID$.uni'))\n\
    fractions_s$ID$.save(os.path.join(path, 'fractions_s$ID$.uni'))\n\
    pressure_s$ID$.save(os.path.join(path, 'pressure_s$ID$.uni'))\n\
    \n\
    vel_s$ID$.save(os.path.join(path, 'vel_s$ID$.uni'))\n\
    velOld_s$ID$.save(os.path.join(path, 'velOld_s$ID$.uni'))\n\
    velParts_s$ID$.save(os.path.join(path, 'velParts_s$ID$.uni'))\n\
    mapWeights_s$ID$.save(os.path.join(path, 'mapWeights_s$ID$.uni'))\n\
    \n\
    x_vel_s$ID$.save(os.path.join(path, 'x_vel_s$ID$.uni'))\n\
    y_vel_s$ID$.save(os.path.join(path, 'y_vel_s$ID$.uni'))\n\
    z_vel_s$ID$.save(os.path.join(path, 'z_vel_s$ID$.uni'))\n\
    x_obvel_s$ID$.save(os.path.join(path, 'x_obvel_s$ID$.uni'))\n\
    y_obvel_s$ID$.save(os.path.join(path, 'y_obvel_s$ID$.uni'))\n\
    z_obvel_s$ID$.save(os.path.join(path, 'z_obvel_s$ID$.uni'))\n\
    \n\
    pp_s$ID$.save(os.path.join(path, 'pp_s$ID$.uni'))\n\
    pVel_s$ID$.save(os.path.join(path, 'pVel_s$ID$.uni'))\n\
    \n\
    gpi_s$ID$.save(os.path.join(path, 'gpi_s$ID$.uni'))\n";

const std::string liquid_export_high = "\n\
def save_liquid_data_high_$ID$(path):\n\
    flags_xl$ID$.save(os.path.join(path, 'flags_xl$ID$.uni'))\n\
    \n\
    phiParts_xl$ID$.save(os.path.join(path, 'phiParts_xl$ID$.uni'))\n\
    phi_xl$ID$.save(os.path.join(path, 'phi_xl$ID$.uni'))\n\
    \n\
    pp_xl$ID$.save(os.path.join(path, 'pp_xl$ID$.uni'))\n";

//////////////////////////////////////////////////////////////////////
// DESTRUCTION
//////////////////////////////////////////////////////////////////////

const std::string liquid_delete_grids_low = "\n\
mantaMsg('Deleting lowres grids, mesh, particlesystem')\n\
if 'flags_s$ID$'      in globals() : del flags_s$ID$\n\
if 'numObs_s$ID$'     in globals() : del numObs_s$ID$\n\
if 'phiParts_s$ID$'   in globals() : del phiParts_s$ID$\n\
if 'phi_s$ID$'        in globals() : del phi_s$ID$\n\
if 'phiIn_s$ID$'      in globals() : del phiIn_s$ID$\n\
if 'phiOut_s$ID$'     in globals() : del phiOut_s$ID$\n\
if 'pressure_s$ID$'   in globals() : del pressure_s$ID$\n\
if 'vel_s$ID$'        in globals() : del vel_s$ID$\n\
if 'x_vel_s$ID$'      in globals() : del x_vel_s$ID$\n\
if 'y_vel_s$ID$'      in globals() : del y_vel_s$ID$\n\
if 'z_vel_s$ID$'      in globals() : del z_vel_s$ID$\n\
if 'obvel_s$ID$'      in globals() : del obvel_s$ID$\n\
if 'x_obvel_s$ID$'    in globals() : del x_obvel_s$ID$\n\
if 'y_obvel_s$ID$'    in globals() : del y_obvel_s$ID$\n\
if 'z_obvel_s$ID$'    in globals() : del z_obvel_s$ID$\n\
if 'velOld_s$ID$'     in globals() : del velOld_s$ID$\n\
if 'velParts_s$ID$'   in globals() : del velParts_s$ID$\n\
if 'mapWeights_s$ID$' in globals() : del mapWeights_s$ID$\n\
if 'pp_s$ID$'         in globals() : del pp_s$ID$\n\
if 'pVel_s$ID$'       in globals() : del pVel_s$ID$\n\
if 'mesh_s$ID$'       in globals() : del mesh_s$ID$\n\
if 'pindex_s$ID$'     in globals() : del pindex_s$ID$\n\
if 'gpi_s$ID$'        in globals() : del gpi_s$ID$\n\
if 'forces_s$ID$'     in globals() : del forces_s$ID$\n\
if 'x_force_s$ID$'    in globals() : del x_force_s$ID$\n\
if 'y_force_s$ID$'    in globals() : del y_force_s$ID$\n\
if 'z_force_s$ID$'    in globals() : del z_force_s$ID$\n\
if 'phiObs_s$ID$'     in globals() : del phiObs_s$ID$\n\
if 'phiObsIn_s$ID$'   in globals() : del phiObsIn_s$ID$\n\
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
