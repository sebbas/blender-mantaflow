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
// GENERAL SETUP
//////////////////////////////////////////////////////////////////////

const std::string liquid_variables = "\n\
narrowBand       = True\n\
narrowBandWidth  = 3\n\
combineBandWidth = narrowBandWidth - 1\n\
\n\
minParticles   = pow(2,dim)\n\
particleNumber = 2\n\
\n\
gravity = (0,0,-1)\n\
step    = -1\n\
maxVel  = 0\n";

//////////////////////////////////////////////////////////////////////
// GRIDS & MESH & PARTICLESYSTEM
//////////////////////////////////////////////////////////////////////

const std::string alloc_liquid = "\n\
flags      = s.create(FlagGrid)\n\
\n\
phiParts   = s.create(LevelsetGrid)\n\
phi        = s.create(LevelsetGrid)\n\
phiTemp    = s.create(LevelsetGrid)\n\
pressure   = s.create(RealGrid)\n\
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

const std::string prep_domain = "\n\
flags.initDomain(boundaryWidth=0)\n\
phi.initFromFlags(flags)\n\
phiTemp.initFromFlags(flags)\n";

//////////////////////////////////////////////////////////////////////
// ADAPTIVE STEP
//////////////////////////////////////////////////////////////////////

const std::string adaptive_step_liquid = "\n\
def manta_step(start_frame):\n\
    s.frame = start_frame\n\
    s.timeTotal = s.frame * dt0\n\
    last_frame = s.frame\n\
    \n\
    # Sample particles on first frame\n\
    if (start_frame == 1):\n\
        phi.copyFrom(phiTemp)\n\
        flags.updateFromLevelset(phi)\n\
        sampleLevelsetWithParticles( phi=phi, flags=flags, parts=pp, discretization=2, randomness=0.01 )\n\
        mapGridToPartsVec3(source=vel, parts=pp, target=pVel )\n\
        phi.save('/Users/sbarschkis/Desktop/phi.uni')\n\
    \n\
    #for i in range(int(gs.z)):\n\
        #phiTemp.printGrid(zSlice=int(i))\n\
    while s.frame == last_frame:\n\
        global step\n\
        step = step + 1\n\
        mantaMsg('Adapt timestep')\n\
        maxvel = vel.getMaxValue()\n\
        s.adaptTimestep(maxvel)\n\
        \n\
        mantaMsg('Liquid step / s.frame: ' + str(s.frame))\n\
        liquid_step()\n\
        s.step()\n";

//////////////////////////////////////////////////////////////////////
// STEP FUNCTIONS
//////////////////////////////////////////////////////////////////////

const std::string liquid_step = "\n\
def liquid_step():\n\
    # Advect particles and grid phi\n\
    # Note: Grid velocities are extrapolated at the end of each step\n\
    pp.advectInGrid(flags=flags, vel=vel, integrationMode=IntRK4, deleteInObstacle=False )\n\
    advectSemiLagrange(flags=flags, vel=vel, grid=phi, order=1)\n\
    flags.updateFromLevelset(phi)\n\
    \n\
    # Advect grid velocity\n\
    if narrowBand:\n\
        advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2)\n\
    \n\
    # Create level set of particles\n\
    gridParticleIndex( parts=pp , flags=flags, indexSys=pindex, index=gpi )\n\
    unionParticleLevelset( pp, pindex, flags, gpi, phiParts )\n\
    \n\
    if narrowBand:\n\
        # Combine level set of particles with grid level set\n\
        phi.addConst(1.); # shrink slightly\n\
        phi.join( phiParts );\n\
        extrapolateLsSimple(phi=phi, distance=narrowBandWidth+2, inside=True )\n\
    else:\n\
        # Overwrite grid level set with level set of particles\n\
        phi.copyFrom( phiParts );\n\
        extrapolateLsSimple(phi=phi, distance=4, inside=True )\n\
    \n\
    extrapolateLsSimple(phi=phi, distance=3 )\n\
    flags.updateFromLevelset(phi)\n\
    \n\
    # Make sure we have velocities throught liquid region\n\
    if narrowBand:\n\
        # Combine particles velocities with advected grid velocities\n\
        mapPartsToMAC(vel=velParts, flags=flags, velOld=velOld, parts=pp, partVel=pVel, weight=mapWeights)\n\
        extrapolateMACFromWeight( vel=velParts , distance=2, weight=mapWeights )\n\
        combineGridVel(vel=velParts, weight=mapWeights , combineVel=vel, phi=phi, narrowBand=combineBandWidth, thresh=0)\n\
        velOld.copyFrom(vel)\n\
    else:\n\
        # Map particle velocities to grid\n\
        mapPartsToMAC(vel=vel, flags=flags, velOld=velOld, parts=pp, partVel=pVel, weight=mapWeights)\n\
        extrapolateMACFromWeight( vel=vel , distance=2, weight=mapWeights )\n\
    \n\
    # Forces & pressure solve\n\
    addGravity(flags=flags, vel=vel, gravity=gravity)\n\
    setWallBcs(flags=flags, vel=vel)\n\
    solvePressure(flags=flags, vel=vel, pressure=pressure, phi=phi)\n\
    setWallBcs(flags=flags, vel=vel)\n\
    \n\
    # Extrapolate velocities\n\
    extrapolateMACSimple( flags=flags, vel=vel, distance=(int(maxVel*1.25 + 2.)) )\n\
    \n\
    # Update particle velocities\n\
    flipVelocityUpdate(vel=vel, velOld=velOld, flags=flags, parts=pp, partVel=pVel, flipRatio=0.95 )\n\
    \n\
    if dim==3:\n\
        phi.createMesh(mesh)\n\
        mesh.save('/Users/sbarschkis/Desktop/surface/fluidsurface_final_%04d.bobj.gz' % step)\n\
    \n\
    # Resample particles\n\
    pVel.setSource( vel, isMAC=True ) # Set source grids for resampling, used in adjustNumber!\n\
    if narrowBand:\n\
        phi.setBoundNeumann(0) # make sure no particles are placed at outer boundary\n\
        adjustNumber( parts=pp, vel=vel, flags=flags, minParticles=1*minParticles, maxParticles=2*minParticles, phi=phi, narrowBand=narrowBandWidth )\n\
    else:\n\
        adjustNumber( parts=pp, vel=vel, flags=flags, minParticles=1*minParticles, maxParticles=2*minParticles, phi=phi )\n";

//////////////////////////////////////////////////////////////////////
// DESTRUCTION
//////////////////////////////////////////////////////////////////////

const std::string del_liquid_grids = "\n\
mantaMsg('Deleting grids, mesh, particlesystem')\n\
if 'flags'      in globals() : del flags\n\
if 'phiParts'   in globals() : del phiParts\n\
if 'phi'        in globals() : del phi\n\
if 'phiTemp'    in globals() : del phiTemp\n\
if 'pressure'   in globals() : del pressure\n\
if 'vel'        in globals() : del vel\n\
if 'velOld'     in globals() : del velOld\n\
if 'velParts'   in globals() : del velParts\n\
if 'mapWeights' in globals() : del mapWeights\n\
if 'pp'         in globals() : del pp\n\
if 'pVel'       in globals() : del pVel\n\
if 'mesh'       in globals() : del mesh\n\
if 'pindex'     in globals() : del pindex\n\
if 'gpi'        in globals() : del gpi\n";

const std::string del_liquid_vars = "\n\
mantaMsg('Deleting liquid variables')\n\
if 'narrowBand'       in globals() : del narrowBand\n\
if 'narrowBandWidth'  in globals() : del narrowBandWidth\n\
if 'combineBandWidth' in globals() : del combineBandWidth\n\
if 'minParticles'     in globals() : del minParticles\n\
if 'particleNumber'   in globals() : del particleNumber\n\
if 'gravity'          in globals() : del gravity\n\
if 'step'             in globals() : del step\n\
if 'maxVel'           in globals() : del maxVel\n";

