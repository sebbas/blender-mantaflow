#
# Simple example scene for a 3D simulation
# Simulation of a flame with smoke
#
from manta import *

# solver params
dim = 3
res = 32
gs = vec3(res, res, res)
if (dim==2):
	gs.z=1
s = Solver(name='main', gridSize = gs, dim=dim)
s.timestep = 1.0
timings = Timings()

# prepare grids
flags = s.create(FlagGrid)
vel = s.create(MACGrid)
density = s.create(RealGrid) # smoke
react = s.create(RealGrid)
fuel = s.create(RealGrid)
heat = s.create(RealGrid)
red = s.create(RealGrid)
green = s.create(RealGrid)
blue = s.create(RealGrid)
flame = s.create(RealGrid)
pressure = s.create(RealGrid)

# noise field
noise = s.create(NoiseField, loadFromFile=True)
noise.posScale = vec3(45)
noise.clamp = True
noise.clampNeg = 0
noise.clampPos = 1
noise.valScale = 1
noise.valOffset = 0.75
noise.timeAnim = 0.2

# initialize domain with boundary
bWidth=0
flags.initDomain(boundaryWidth=bWidth)
flags.fillGrid()

#setOpenBound(flags,bWidth,'yY',FlagOutflow|FlagEmpty)

if (GUI):
	gui = Gui()
	gui.show(True)
	gui.pause()

# Cube in center of domain (x, y), standing on bottom (z=0)
boxSize = vec3(res/8, res/8, res/8)
boxCenter = vec3(res/2, boxSize.z+1, res/2)
source = s.create(Box, center=boxCenter, size=boxSize)

# main loop
for t in range(1000):
	
	print('Inflow')
	burnStep(flags=flags, density=density, heat=heat, fuel=fuel, react=react, red=red, green=green, blue=blue, flame=flame, shape=source, gridSize=gs);
	
	print('Advecting density')
	#advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2)
	
	print('Advecting react')
	#advectSemiLagrange(flags=flags, vel=vel, grid=react, order=2)

	print('Advecting fuel')
	advectSemiLagrange(flags=flags, vel=vel, grid=fuel, order=2)

	print('Advecting velocity')
	advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2, strength=1.0)

	print('Walls')
	setWallBcs(flags=flags, vel=vel)
	
	print('Buoyancy')
	addBuoyancy(density=density, vel=vel, gravity=vec3(0,-6e-4,0), flags=flags)
	
	print('Pressure')
	solvePressure(flags=flags, vel=vel, pressure=pressure)
	
	print('Walls')
	setWallBcs(flags=flags , vel=vel)

	timings.display()
	s.step()
