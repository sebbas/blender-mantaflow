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
flags.initDomain(boundaryWidth=boundaryWidth)\n\
if doOpen:\n\
    setOpenBound(flags=flags, bWidth=boundaryWidth, openBound=boundConditions, type=FlagOutflow|FlagEmpty)\n";

const std::string liquid_bounds_high = "\n\
# prepare domain high\n\
mantaMsg('Liquid domain high')\n\
xl_flags.initDomain(boundaryWidth=boundaryWidth)\n\
if doOpen:\n\
    setOpenBound(flags=xl_flags, bWidth=boundaryWidth, openBound=boundConditions, type=FlagOutflow|FlagEmpty)\n";

//////////////////////////////////////////////////////////////////////
// VARIABLES
//////////////////////////////////////////////////////////////////////

const std::string liquid_variables_low = "\n\
narrowBand       = True\n\
narrowBandWidth  = 3\n\
combineBandWidth = narrowBandWidth - 1\n\
\n\
minParticles   = pow(2,dim)\n\
particleNumber = 2\n\
radiusFactor   = 1.0\n\
randomness     = $RANDOMNESS$\n\
\n\
step    = -1\n\
maxVel  = 0\n\
\n\
using_highres = $USE_WAVELETS$\n";

const std::string liquid_variables_high = "\n\
scale = 0.5\n\
xl_radiusFactor = 2.5\n";

//////////////////////////////////////////////////////////////////////
// GRIDS & MESH & PARTICLESYSTEM
//////////////////////////////////////////////////////////////////////

const std::string liquid_alloc_low = "\n\
flags      = s.create(FlagGrid)\n\
\n\
phiParts   = s.create(LevelsetGrid)\n\
phi        = s.create(LevelsetGrid)\n\
phiInit    = s.create(LevelsetGrid)\n\
pressure   = s.create(RealGrid)\n\
# TODO (sebbas): remove density grid. only here because of stand in\n\
density    = s.create(RealGrid)\n\
\n\
vel        = s.create(MACGrid)\n\
velOld     = s.create(MACGrid)\n\
velParts   = s.create(MACGrid)\n\
mapWeights = s.create(MACGrid)\n\
\n\
pp         = s.create(BasicParticleSystem)\n\
pVel       = pp.create(PdataVec3)\n\
mesh       = s.create(Mesh)\n\
\n\
# Acceleration data for particle nbs\n\
pindex     = s.create(ParticleIndexSystem)\n\
gpi        = s.create(IntGrid)\n";

const std::string liquid_alloc_high = "\n\
xl_flags   = xl.create(FlagGrid)\n\
xl_phi     = xl.create(LevelsetGrid)\n\
xl_pp      = xl.create(BasicParticleSystem)\n\
xl_mesh    = xl.create(Mesh)\n\
\n\
# Acceleration data for particle nbs\n\
xl_pindex  = xl.create(ParticleIndexSystem)\n\
xl_gpi     = xl.create(IntGrid)\n";

const std::string liquid_init_phi = "\n\
phi.initFromFlags(flags)\n\
phiInit.initFromFlags(flags)\n";

//////////////////////////////////////////////////////////////////////
// STEP FUNCTIONS
//////////////////////////////////////////////////////////////////////

const std::string liquid_adaptive_step = "\n\
def manta_step(start_frame):\n\
    s.frame = start_frame\n\
    s.timeTotal = s.frame * dt0\n\
    last_frame = s.frame\n\
    \n\
    sampleLevelsetWithParticles( phi=phiInit, flags=flags, parts=pp, discretization=2, randomness=randomness, refillEmpty=True )\n\
    mapGridToPartsVec3(source=vel, parts=pp, target=pVel )\n\
    phi.join(phiInit)\n\
    flags.updateFromLevelset(phi)\n\
    \n\
    while s.frame == last_frame:\n\
        mantaMsg('Adapt timestep')\n\
        maxvel = vel.getMaxValue()\n\
        s.adaptTimestep(maxvel)\n\
        \n\
        mantaMsg('Low step / s.frame: ' + str(s.frame))\n\
        liquid_step()\n\
        \n\
        # TODO (sebbas)\n\
        if using_highres:\n\
            xl.timestep = s.timestep\n\
            mantaMsg('High step / s.frame: ' + str(s.frame))\n\
            liquid_step_high()\n\
        \n\
        s.step()\n";

const std::string liquid_step_low = "\n\
def liquid_step():\n\
    # Advect particles and grid phi\n\
    # Note: Grid velocities are extrapolated at the end of each step\n\
    pp.advectInGrid(flags=flags, vel=vel, integrationMode=IntRK4, deleteInObstacle=False)\n\
    advectSemiLagrange(flags=flags, vel=vel, grid=phi, order=1, openBounds=doOpen, boundaryWidth=boundaryWidth)\n\
    flags.updateFromLevelset(phi)\n\
    \n\
    # Advect grid velocity\n\
    if narrowBand:\n\
        advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2, openBounds=doOpen, boundaryWidth=boundaryWidth)\n\
    \n\
    # Create level set of particles\n\
    gridParticleIndex(parts=pp , flags=flags, indexSys=pindex, index=gpi)\n\
    unionParticleLevelset(pp, pindex, flags, gpi, phiParts)\n\
    \n\
    if doOpen:\n\
        resetOutflow(flags=flags, phi=phi, parts=pp, index=gpi, indexSys=pindex)\n\
        flags.updateFromLevelset(phi)\n\
    \n\
    if narrowBand:\n\
        # Combine level set of particles with grid level set\n\
        phi.addConst(1.) # shrink slightly\n\
        phi.join(phiParts)\n\
        extrapolateLsSimple(phi=phi, distance=narrowBandWidth+2, inside=True )\n\
    else:\n\
        # Overwrite grid level set with level set of particles\n\
        phi.copyFrom(phiParts)\n\
        extrapolateLsSimple(phi=phi, distance=4, inside=True )\n\
    \n\
    extrapolateLsSimple(phi=phi, distance=3)\n\
    flags.updateFromLevelset(phi)\n\
    \n\
    # Make sure we have velocities throught liquid region\n\
    if narrowBand:\n\
        # Combine particles velocities with advected grid velocities\n\
        mapPartsToMAC(vel=velParts, flags=flags, velOld=velOld, parts=pp, partVel=pVel, weight=mapWeights)\n\
        extrapolateMACFromWeight(vel=velParts , distance=2, weight=mapWeights)\n\
        combineGridVel(vel=velParts, weight=mapWeights , combineVel=vel, phi=phi, narrowBand=combineBandWidth, thresh=0)\n\
        velOld.copyFrom(vel)\n\
    else:\n\
        # Map particle velocities to grid\n\
        mapPartsToMAC(vel=vel, flags=flags, velOld=velOld, parts=pp, partVel=pVel, weight=mapWeights)\n\
        extrapolateMACFromWeight(vel=vel , distance=2, weight=mapWeights)\n\
    \n\
    # Forces & pressure solve\n\
    addGravity(flags=flags, vel=vel, gravity=gravity)\n\
    setWallBcs(flags=flags, vel=vel)\n\
    solvePressure(flags=flags, vel=vel, pressure=pressure, phi=phi)\n\
    setWallBcs(flags=flags, vel=vel)\n\
    \n\
    # Extrapolate velocities\n\
    extrapolateMACSimple(flags=flags, vel=vel, distance=(int(maxVel*1.25 + 2.)))\n\
    \n\
    # Update particle velocities\n\
    flipVelocityUpdate(vel=vel, velOld=velOld, flags=flags, parts=pp, partVel=pVel, flipRatio=0.95)\n\
    \n\
    if dim==3:\n\
        phi.createMesh(mesh)\n\
    \n\
    # Resample particles\n\
    pVel.setSource(vel, isMAC=True) # Set source grids for resampling, used in adjustNumber!\n\
    if narrowBand:\n\
        phi.setBoundNeumann(boundaryWidth) # make sure no particles are placed at outer boundary\n\
        adjustNumber(parts=pp, vel=vel, flags=flags, minParticles=1*minParticles, maxParticles=2*minParticles, phi=phi, radiusFactor=radiusFactor, narrowBand=narrowBandWidth)\n\
    else:\n\
        adjustNumber(parts=pp, vel=vel, flags=flags, minParticles=1*minParticles, maxParticles=2*minParticles, phi=phi, radiusFactor=radiusFactor)\n\
    \n\
    # TODO (sebbas): HACK - saving particle system for highres step\n\
    #if using_highres:\n\
        #pp.save('/tmp/partfile.uni')\n\
    \n\
    # reset inflow grid\n\
    phiInit.setConst(0.5)\n";

const std::string liquid_step_high = "\n\
def liquid_step_high():\n\
    xl_phi.setBound(value=0., boundaryWidth=1)\n\
    xl_pp.load('/tmp/partfile.uni')\n\
    \n\
    # create surface\n\
    gridParticleIndex( parts=xl_pp , flags=xl_flags, indexSys=xl_pindex, index=xl_gpi )\n\
    unionParticleLevelset( xl_pp, xl_pindex, xl_flags, xl_gpi, xl_phi , xl_radiusFactor )\n\
    #averagedParticleLevelset( xl_pp, xl_pindex, xl_flags, xl_gpi, xl_phi , xl_radiusFactor , 1, 1 )\n\
    \n\
    xl_phi.setBound(value=0., boundaryWidth=1)\n\
    xl_phi.createMesh(xl_mesh)\n\
    \n\
    # beautify mesh, too slow right now!\n\
    #subdivideMesh(mesh=xl_mesh, minAngle=0.01, minLength=scale, maxLength=3*scale, cutTubes=False)\n\
    # perform smoothing\n\
    #for iters in range(10):\n\
        #smoothMesh(mesh=xl_mesh, strength=1e-3, steps=10)\n\
        #subdivideMesh(mesh=xl_mesh, minAngle=0.01, minLength=scale, maxLength=3*scale, cutTubes=True)\n";

//////////////////////////////////////////////////////////////////////
// IMPORT / EXPORT
//////////////////////////////////////////////////////////////////////

const std::string liquid_save_mesh = "\n\
def save_mesh(path):\n\
    mesh.save(path)\n\
    # TODO (sebbas)\n\
	#if using_highres:\n\
        #xl_mesh.save(path)\n";

const std::string liquid_import_low = "\n\
def load_liquid_data(path):\n\
    flags.load(path + str('flags.uni'))\n\
    \n\
    phiParts.load(path + str('phiParts.uni'))\n\
    phi.load(path + str('phi.uni'))\n\
    phiInit.load(path + str('phiInit.uni'))\n\
    pressure.load(path + str('pressure.uni'))\n\
    \n\
    vel.load(path + str('vel.uni'))\n\
    velOld.load(path + str('velOld.uni'))\n\
    velParts.load(path + str('velParts.uni'))\n\
    mapWeights.load(path + str('mapWeights.uni'))\n\
    \n\
    pp.load(path + str('pp.uni'))\n\
    pVel.load(path + str('pVel.uni'))\n\
    \n\
    gpi.load(path + str('gpi.uni'))\n";

const std::string liquid_export_low = "\n\
def save_liquid_data(path):\n\
    flags.save(path + str('flags.uni'))\n\
    \n\
    phiParts.save(path + str('phiParts.uni'))\n\
    phi.save(path + str('phi.uni'))\n\
    phiInit.save(path + str('phiInit.uni'))\n\
    pressure.save(path + str('pressure.uni'))\n\
    \n\
    vel.save(path + str('vel.uni'))\n\
    velOld.save(path + str('velOld.uni'))\n\
    velParts.save(path + str('velParts.uni'))\n\
    mapWeights.save(path + str('mapWeights.uni'))\n\
    \n\
    pp.save(path + str('pp.uni'))\n\
    pVel.save(path + str('pVel.uni'))\n\
    \n\
    gpi.save(path + str('gpi.uni'))\n";

//////////////////////////////////////////////////////////////////////
// DESTRUCTION
//////////////////////////////////////////////////////////////////////

const std::string liquid_delete_grids_low = "\n\
mantaMsg('Deleting lowres grids, mesh, particlesystem')\n\
if 'flags'      in globals() : del flags\n\
if 'phiParts'   in globals() : del phiParts\n\
if 'phi'        in globals() : del phi\n\
if 'phiInit'    in globals() : del phiInit\n\
if 'pressure'   in globals() : del pressure\n\
if 'density'    in globals() : del density\n\
if 'vel'        in globals() : del vel\n\
if 'velOld'     in globals() : del velOld\n\
if 'velParts'   in globals() : del velParts\n\
if 'mapWeights' in globals() : del mapWeights\n\
if 'pp'         in globals() : del pp\n\
if 'pVel'       in globals() : del pVel\n\
if 'mesh'       in globals() : del mesh\n\
if 'pindex'     in globals() : del pindex\n\
if 'gpi'        in globals() : del gpi\n";

const std::string liquid_delete_grids_high = "\n\
mantaMsg('Deleting highres grids, mesh, particlesystem')\n\
if 'xl_flags'   in globals() : del xl_flags\n\
if 'xl_phi'     in globals() : del xl_phi\n\
if 'xl_pp'      in globals() : del xl_pp\n\
if 'xl_mesh'    in globals() : del xl_mesh\n\
if 'xl_pindex'  in globals() : del xl_pindex\n\
if 'xl_gpi'     in globals() : del xl_gpi\n";

const std::string liquid_delete_variables_low = "\n\
mantaMsg('Deleting lowres liquid variables')\n\
if 'narrowBand'       in globals() : del narrowBand\n\
if 'narrowBandWidth'  in globals() : del narrowBandWidth\n\
if 'combineBandWidth' in globals() : del combineBandWidth\n\
if 'minParticles'     in globals() : del minParticles\n\
if 'particleNumber'   in globals() : del particleNumber\n\
if 'step'             in globals() : del step\n\
if 'maxVel'           in globals() : del maxVel\n";

const std::string liquid_delete_variables_high = "\n\
mantaMsg('Deleting highres liquid variables')\n\
if 'scale'            in globals() : del scale\n\
if 'xl_radiusFactor'  in globals() : del xl_radiusFactor\n";

