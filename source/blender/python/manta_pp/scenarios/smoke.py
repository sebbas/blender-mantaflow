from manta import * 
import os, shutil, math, sys 
def transform_back(obj, gs):
	obj.scale(gs/2)
	obj.offset(gs/2)

uvs = $UVS_CNT$
velInflow = vec3(0, 0, 1)
upres = $UPRES$
wltStrength = $WLT_STR$
octaves = int( math.log(upres)/ math.log(2.0) + 0.5 ) 
res = $RES$
gs = vec3($RESX$, $RESY$, $RESZ$) 
s = Solver(name = 'main', gridSize = gs, dim = $SOLVER_DIM$) 
s.timestep = 1 
noise = s.create(NoiseField, fixedSeed=256, loadFromFile=True) 
noise.posScale = vec3(20) 
noise.clamp = False 
noise.clampNeg = $NOISE_CN$
noise.clampPos = $NOISE_CP$
noise.valScale = $NOISE_VALSCALE$
noise.valOffset = $NOISE_VALOFFSET$
noise.timeAnim = 0.2 
source = s.create(Mesh)
source.load('manta_flow.obj')
transform_back(source, gs)
sourceVel = s.create(Mesh)
sourceVel.load('manta_flow.obj')
transform_back(sourceVel, gs)
xl_gs = vec3(42, 42, 64) 
xl = Solver(name = 'larger', gridSize = xl_gs, dim = 3) 
xl.timestep = 0.5 
xl_vel = xl.create(MACGrid) 
xl_density = xl.create(RealGrid) 
xl_flags = xl.create(FlagGrid) 
xl_flags.initDomain() 
xl_flags.fillGrid() 
xl_source = s.create(Mesh)
xl_source.load('manta_flow.obj')
transform_back(xl_source, gs)
xl_noise = xl.create(NoiseField, fixedSeed=256, loadFromFile=True) 
xl_noise.posScale = vec3(20) 
xl_noise.clamp = False 
xl_noise.clampNeg = 0 
xl_noise.clampPos = 1 
xl_noise.valScale = 0 
xl_noise.valOffset = 0.075 
xl_noise.timeAnim = 0.4 
flags = s.create(FlagGrid) 
flags.initDomain() 
flags.fillGrid() 
uv = [] 
for i in range(uvs): 
	uvGrid = s.create(VecGrid) 
	uv.append(uvGrid) 
	resetUvGrid( uv[i] ) 
vel = s.create(MACGrid) 
density = s.create(RealGrid) 
pressure = s.create(RealGrid) 
energy = s.create(RealGrid) 
tempFlag  = s.create(FlagGrid)
sdf_flow  = s.create(LevelsetGrid)
forces = s.create(MACGrid)
source.meshSDF(source, sdf_flow, 1.1)
source_shape = s.create(Cylinder, center=gs*vec3(0.5,0.1,0.5), radius=res*0.14, z=gs*vec3(0, 0.02, 0))
xl_wltnoise = s.create(NoiseField, loadFromFile=True) 
xl_wltnoise.posScale = vec3( int(1.0*gs.x) ) * 0.5 
xl_wltnoise.posScale = xl_wltnoise.posScale * 0.5
xl_wltnoise.timeAnim = 0.1 
def sim_step(t):
	forces.load('manta_forces.uni')
	addForceField(flags=flags, vel=vel,force=forces)
	addBuoyancy(density=density, vel=vel, gravity=vec3(0,0,-0.1), flags=flags) 
	advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2) 
	advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2) 
	for i in range(uvs): 
		advectSemiLagrange(flags=flags, vel=vel, grid=uv[i], order=2) 
		updateUvWeight( resetTime=16.5 , index=i, numUvs=uvs, uv=uv[i] )
	applyInflow=False
	if (t>=0 and t<75):
		densityInflowMesh(flags=flags, density=density, mesh=source, value=1)
		applyInflow=True
	setWallBcs(flags=flags, vel=vel) 
	vorticityConfinement( vel=vel, flags=flags, strength=0.2 ) 
	solvePressure(flags=flags, vel=vel, pressure=pressure, useResNorm=True, openBound='xXyYzZ', cgMaxIterFac=1, cgAccuracy=0.01) 
	setWallBcs(flags=flags, vel=vel) 
	computeEnergy(flags=flags, vel=vel, energy=energy)
	tempFlag.copyFrom(flags)
	extrapolateSimpleFlags( flags=flags, val=tempFlag, distance=2, flagFrom=FlagObstacle, flagTo=FlagFluid )
	extrapolateSimpleFlags( flags=tempFlag, val=energy, distance=6, flagFrom=FlagFluid, flagTo=FlagObstacle )
	computeWaveletCoeffs(energy)
	density.save('den%04d_temp.uni' % t) 
	os.rename('den%04d_temp.uni' % t, 'den%04d.uni' % t) 
	s.step()
	
	interpolateMACGrid( source=vel, target=xl_vel ) 
	sStr = 1.0 * wltStrength  
	sPos = 2.0  
	for o in range(octaves): 
		for i in range(uvs): 
			uvWeight = getUvWeight(uv[i])  
			applyNoiseVec3( flags=xl_flags, target=xl_vel, noise=xl_wltnoise, scale=sStr * uvWeight, scaleSpatial=sPos , weight=energy, uv=uv[i] ) 
		sStr *= 0.06 # magic kolmogorov factor 
		sPos *= 2.0 
	for substep in range(upres):  
		advectSemiLagrange(flags=xl_flags, vel=xl_vel, grid=xl_density, order=2)  
	if (applyInflow): 
		densityInflowMesh(flags=xl_flags, density=xl_density, mesh=source, value=1)
	xl_density.save('densityXl_%04d.uni' % t)
	xl.step()   
