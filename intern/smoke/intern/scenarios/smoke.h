#include <string>
using namespace std;
const string smoke_setup = "\n\
from manta import * \n\
import os, shutil, math, sys \n\
def transform_back(obj, gs):\n\
	obj.scale(gs/2)\n\
	obj.offset(gs/2)\n\
\n\
uvs = $UVS_CNT$\n\
velInflow = vec3(0, 0, 1)\n\
upres = $UPRES$\n\
wltStrength = $WLT_STR$\n\
octaves = int( math.log(upres)/ math.log(2.0) + 0.5 ) \n\
res = $RES$\n\
gs = vec3($RESX$, $RESY$, $RESZ$) \n\
s = Solver(name = 'main', gridSize = gs, dim = $SOLVER_DIM$) \n\
s.timestep = 1 \n\
noise = s.create(NoiseField, fixedSeed=256, loadFromFile=True) \n\
noise.posScale = vec3(20) \n\
noise.clamp = False \n\
noise.clampNeg = $NOISE_CN$\n\
noise.clampPos = $NOISE_CP$\n\
noise.valScale = $NOISE_VALSCALE$\n\
noise.valOffset = $NOISE_VALOFFSET$\n\
noise.timeAnim = 0.2 \n\
source = s.create(Mesh)\n\
source.load('manta_flow.obj')\n\
transform_back(source, gs)\n\
sourceVel = s.create(Mesh)\n\
sourceVel.load('manta_flow.obj')\n\
transform_back(sourceVel, gs)\n\
xl_gs = vec3(42, 42, 64) \n\
xl = Solver(name = 'larger', gridSize = xl_gs, dim = 3) \n\
xl.timestep = 0.5 \n\
xl_vel = xl.create(MACGrid) \n\
xl_density = xl.create(RealGrid) \n\
xl_flags = xl.create(FlagGrid) \n\
xl_flags.initDomain() \n\
xl_flags.fillGrid() \n\
xl_source = s.create(Mesh)\n\
xl_source.load('manta_flow.obj')\n\
transform_back(xl_source, gs)\n\
xl_noise = xl.create(NoiseField, fixedSeed=256, loadFromFile=True) \n\
xl_noise.posScale = vec3(20) \n\
xl_noise.clamp = False \n\
xl_noise.clampNeg = 0 \n\
xl_noise.clampPos = 1 \n\
xl_noise.valScale = 0 \n\
xl_noise.valOffset = 0.075 \n\
xl_noise.timeAnim = 0.4 \n\
flags = s.create(FlagGrid) \n\
flags.initDomain() \n\
flags.fillGrid() \n\
uv = [] \n\
for i in range(uvs): \n\
	uvGrid = s.create(VecGrid) \n\
	uv.append(uvGrid) \n\
	resetUvGrid( uv[i] ) \n\
vel = s.create(MACGrid) \n\
density = s.create(RealGrid) \n\
pressure = s.create(RealGrid) \n\
energy = s.create(RealGrid) \n\
tempFlag  = s.create(FlagGrid)\n\
sdf_flow  = s.create(LevelsetGrid)\n\
forces = s.create(MACGrid)\n\
source.meshSDF(source, sdf_flow, 1.1)\n\
source_shape = s.create(Cylinder, center=gs*vec3(0.5,0.1,0.5), radius=res*0.14, z=gs*vec3(0, 0.02, 0))\n\
xl_wltnoise = s.create(NoiseField, loadFromFile=True) \n\
xl_wltnoise.posScale = vec3( int(1.0*gs.x) ) * 0.5 \n\
xl_wltnoise.posScale = xl_wltnoise.posScale * 0.5\n\
xl_wltnoise.timeAnim = 0.1 \n\
def sim_step(t):\n\
	forces.load('manta_forces.uni')\n\
	addForceField(flags=flags, vel=vel,force=forces)\n\
	addBuoyancy(density=density, vel=vel, gravity=vec3(0,0,-0.1), flags=flags) \n\
	advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2) \n\
	advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2) \n\
	for i in range(uvs): \n\
		advectSemiLagrange(flags=flags, vel=vel, grid=uv[i], order=2) \n\
		updateUvWeight( resetTime=16.5 , index=i, numUvs=uvs, uv=uv[i] )\n\
	applyInflow=False\n\
	if (t>=0 and t<75):\n\
		densityInflowMesh(flags=flags, density=density, mesh=source, value=1)\n\
		applyInflow=True\n\
	setWallBcs(flags=flags, vel=vel) \n\
	vorticityConfinement( vel=vel, flags=flags, strength=0.2 ) \n\
	solvePressure(flags=flags, vel=vel, pressure=pressure, useResNorm=True, openBound='xXyYzZ', cgMaxIterFac=1, cgAccuracy=0.01) \n\
	setWallBcs(flags=flags, vel=vel) \n\
	computeEnergy(flags=flags, vel=vel, energy=energy)\n\
	tempFlag.copyFrom(flags)\n\
	extrapolateSimpleFlags( flags=flags, val=tempFlag, distance=2, flagFrom=FlagObstacle, flagTo=FlagFluid )\n\
	extrapolateSimpleFlags( flags=tempFlag, val=energy, distance=6, flagFrom=FlagFluid, flagTo=FlagObstacle )\n\
	computeWaveletCoeffs(energy)\n\
	density.save('den%04d_temp.uni' % t) \n\
	os.rename('den%04d_temp.uni' % t, 'den%04d.uni' % t) \n\
	s.step()\n\
	\n\
	interpolateMACGrid( source=vel, target=xl_vel ) \n\
	sStr = 1.0 * wltStrength  \n\
	sPos = 2.0  \n\
	for o in range(octaves): \n\
		for i in range(uvs): \n\
			uvWeight = getUvWeight(uv[i])  \n\
			applyNoiseVec3( flags=xl_flags, target=xl_vel, noise=xl_wltnoise, scale=sStr * uvWeight, scaleSpatial=sPos , weight=energy, uv=uv[i] ) \n\
		sStr *= 0.06 # magic kolmogorov factor \n\
		sPos *= 2.0 \n\
	for substep in range(upres):  \n\
		advectSemiLagrange(flags=xl_flags, vel=xl_vel, grid=xl_density, order=2)  \n\
	if (applyInflow): \n\
		densityInflowMesh(flags=xl_flags, density=xl_density, mesh=source, value=1)\n\
	xl_density.save('densityXl_%04d.uni' % t)\n\
	xl.step()\n\
";
