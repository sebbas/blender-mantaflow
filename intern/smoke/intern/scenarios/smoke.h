#include <string>
using namespace std;

//////////////////////////////////////////////////////////////////////
// GENERAL SETUP
//////////////////////////////////////////////////////////////////////

const string manta_import = "\
from manta import *\n\
import os, shutil, math, sys\n";

const string flags = "\n\
using_colors = $USING_COLORS$\n\
using_heat = $USING_HEAT$\n\
using_fire = $USING_FIRE$\n\
low_flags_updated = False\n\
using_wavelets = $USE_WAVELETS$\n";

const string uv_setup = "\n\
# create the array of uv grids\n\
uv = []\n\
for i in range(uvs):\n\
  uvGrid = s.create(VecGrid)\n\
  uv.append(uvGrid)\n\
  resetUvGrid(uv[i])\n";

//////////////////////////////////////////////////////////////////////
// LOW RESOLUTION SETUP
//////////////////////////////////////////////////////////////////////

const string solver_setup_low = "\n\
# solver low params\n\
dim = $SOLVER_DIM$\n\
doOpen = $DO_OPEN$\n\
boundConditions = '$BOUNDCONDITIONS$'\n\
res = $RES$\n\
gs = vec3($RESX$,$RESY$,$RESZ$)\n\
if dim == 2:\n\
  gs.z = 1\n\
s = FluidSolver(name='main', gridSize=gs, dim=dim)\n\
s.timestep = $TIMESTEP$\n\
timings = Timings()\n\
vorticity = $VORTICITY$\n\
uvs = $UVS_CNT$\n";

const string alloc_base_grids_low = "\n\
# prepare grids low\n\
flags = s.create(FlagGrid)\n\
vel = s.create(MACGrid)\n\
x_vel = s.create(RealGrid)\n\
y_vel = s.create(RealGrid)\n\
z_vel = s.create(RealGrid)\n\
density = s.create(LevelsetGrid)\n\
pressure = s.create(RealGrid)\n\
energy = s.create(RealGrid)\n\
forces = s.create(MACGrid)\n\
inflow_grid = s.create(LevelsetGrid)\n\
fuel_inflow = s.create(LevelsetGrid)\n";

const string noise_low = "\n\
# noise field low\n\
noise = s.create(NoiseField, loadFromFile=True)\n\
noise.posScale = vec3(45)\n\
noise.clamp = True\n\
noise.clampNeg = $NOISE_CN$\n\
noise.clampPos = $NOISE_CP$\n\
noise.valScale = $NOISE_VALSCALE$\n\
noise.valOffset = $NOISE_VALOFFSET$\n\
noise.timeAnim = $NOISE_TIMEANIM$\n";

const string prep_domain_low = "\n\
# prepare domain low\n\
flags.initDomain()\n\
flags.fillGrid()\n\
setOpenBound(flags=flags, bWidth=1, openBound=boundConditions, type=FlagOutflow|FlagEmpty)\n";

//////////////////////////////////////////////////////////////////////
// HIGH RESOLUTION SETUP
//////////////////////////////////////////////////////////////////////

const string solver_setup_high = "\n\
# solver high params\n\
upres = $UPRES$\n\
xl_gs = vec3($HRESX$, $HRESY$, $HRESZ$)\n\
if dim == 2:\n\
  xl_gs.z = 1\n\
xl = Solver(name = 'larger', gridSize = xl_gs)\n\
xl.timestep = $XL_TIMESTEP$\n\
wltStrength = $WLT_STR$\n\
octaves = 0\n\
if(upres>0):\n\
  octaves = int(math.log(upres)/ math.log(2.0) + 0.5)\n";

const string alloc_base_grids_high = "\n\
# prepare grids high\n\
xl_flags = xl.create(FlagGrid)\n\
xl_vel = xl.create(MACGrid)\n\
xl_x_vel = s.create(RealGrid)\n\
xl_y_vel = s.create(RealGrid)\n\
xl_z_vel = s.create(RealGrid)\n\
xl_density = xl.create(RealGrid)\n\
xl_weight  = xl.create(RealGrid)\n";

const string noise_high = "\n\
# noise field high\n\
xl_noise = xl.create(NoiseField, fixedSeed=256, loadFromFile=True)\n\
xl_noise.posScale = vec3(20)\n\
xl_noise.clamp = False\n\
xl_noise.clampNeg = $NOISE_CN$\n\
xl_noise.clampPos = $NOISE_CP$\n\
xl_noise.valScale = $NOISE_VALSCALE$\n\
xl_noise.valOffset = $NOISE_VALOFFSET$\n\
xl_noise.timeAnim = $NOISE_TIMEANIM$ * upres\n";

const string prep_domain_high = "\n\
# prepare domain high\n\
xl_flags.initDomain()\n\
xl_flags.fillGrid()\n";

const string wavelet_turbulence_noise = "\n\
# wavelet turbulence noise field\n\
xl_wltnoise = s.create(NoiseField, loadFromFile=True)\n\
xl_wltnoise.posScale = vec3( int(1.0*gs.x) ) * 0.5\n\
xl_wltnoise.timeAnim = 0.1\n\
if(upres>0):\n\
  xl_wltnoise.posScale = xl_wltnoise.posScale * (1./upres)\n";

//////////////////////////////////////////////////////////////////////
// ADDITIONAL GRIDS
//////////////////////////////////////////////////////////////////////

const string alloc_colors_low = "\n\
print('Allocating colors low')\n\
color_r = s.create(RealGrid)\n\
color_g = s.create(RealGrid)\n\
color_b = s.create(RealGrid)\n";

const string alloc_colors_high = "\n\
print('Allocating colors high')\n\
xl_color_r = xl.create(RealGrid)\n\
xl_color_g = xl.create(RealGrid)\n\
xl_color_b = xl.create(RealGrid)\n";

const string init_colors_low = "\n\
print('Initializing colors low')\n\
color_r.copyFrom(density) \n\
color_r.multConst(manta_color_r) \n\
color_g.copyFrom(density) \n\
color_g.multConst(manta_color_g) \n\
color_b.copyFrom(density) \n\
color_b.multConst(manta_color_b)\n";

const string init_colors_high = "\n\
print('Initializing colors high')\n\
xl_color_r.copyFrom(xl_density) \n\
xl_color_r.multConst(manta_color_r) \n\
xl_color_g.copyFrom(xl_density) \n\
xl_color_g.multConst(manta_color_g) \n\
xl_color_b.copyFrom(xl_density) \n\
xl_color_b.multConst(manta_color_b)\n";

const string alloc_heat_low = "\n\
print('Allocating heat low')\n\
heat = s.create(RealGrid)\n";

const string alloc_fire_low = "\n\
print('Allocating fire low')\n\
flame = s.create(RealGrid)\n\
fuel = s.create(RealGrid)\n\
react = s.create(RealGrid)\n";

const string alloc_fire_high = "\n\
print('Allocating fire high')\n\
xl_flame = xl.create(RealGrid)\n\
xl_fuel = xl.create(RealGrid)\n\
xl_react = xl.create(RealGrid)\n";

const string with_heat = "\n\
using_heat = True\n";

const string with_colors = "\n\
using_colors = True\n";

const string with_fire = "\n\
using_fire = True\n";

//////////////////////////////////////////////////////////////////////
// STANDALONE MODE
//////////////////////////////////////////////////////////////////////

const string standalone = "\n\
if (GUI):\n\
  gui=Gui()\n\
  gui.show()\n\
  gui.pause()\n\
\n\
import_grids_low()\n\
if using_wavelets:\n\
  import_grids_high()\n\
\n\
for step in range(1000):\n\
  apply_inflow()\n\
  \n\
  print('Low step '+ str(step))\n\
  if using_fire:\n\
    process_burn_low()\n\
  step_low()\n\
  if using_fire:\n\
    update_flame_low()\n\
  \n\
  print('High step '+ str(step))\n\
  if using_wavelets:\n\
    if using_fire:\n\
      process_burn_high()\n\
    step_high()\n\
    if using_fire:\n\
      update_flame_high()\n";

//////////////////////////////////////////////////////////////////////
// DESCTRUCTION
//////////////////////////////////////////////////////////////////////

const string del_colors_low = "\n\
print('Deleting colors low')\n\
del color_r\n\
del color_g\n\
del color_b\n";

const string del_colors_high = "\n\
print('Deleting colors high')\n\
del xl_color_r\n\
del xl_color_g\n\
del xl_color_b\n";

const string del_fire_low = "\n\
print('Deleting fire low')\n\
del flame\n\
del fuel\n\
del react\n";

const string del_fire_high = "\n\
print('Deleting fire high')\n\
del xl_flame\n\
del xl_fuel\n\
del xl_react\n";

const string del_heat_low = "\n\
print('Deleting heat low')\n\
del heat\n";

const string del_base_grids_low = "\n\
print('Deleting base grids low')\n\
del res\n\
del dim\n\
del gs\n\
del doOpen\n\
del boundConditions\n\
del s\n\
del timings\n\
del using_colors\n\
del using_heat\n\
del using_fire\n\
del flags\n\
del vel\n\
del x_vel\n\
del y_vel\n\
del z_vel\n\
del density\n\
del pressure\n\
del energy\n\
del forces\n\
del inflow_grid\n\
del fuel_inflow\n\
del noise\n";

const string del_base_grids_high = "\n\
print('Deleting base grids high')\n\
del upres\n\
del xl_gs\n\
del xl\n\
del uvs\n\
del wltStrength\n\
del octaves\n\
del xl_flags\n\
del xl_vel\n\
del xl_x_vel\n\
del xl_y_vel\n\
del xl_z_vel\n\
del xl_density\n\
del xl_weight\n\
del xl_noise\n\
del xl_wltnoise\n";

//////////////////////////////////////////////////////////////////////
// STEP FUNCTIONS LOW
//////////////////////////////////////////////////////////////////////

const string smoke_step_low = "\n\
def step_low():\n\
  print('Step low')\n\
  copyRealToMac(sourceX=x_vel, sourceY=y_vel, sourceZ=z_vel, target=vel)\n\
  if dim == 2:\n\
    density.add(inflow_grid)\n\
  \n\
  if using_colors:\n\
    print ('Advecting colors')\n\
    advectSemiLagrange(flags=flags, vel=vel, grid=color_r, order=$ADVECT_ORDER$)\n\
    advectSemiLagrange(flags=flags, vel=vel, grid=color_g, order=$ADVECT_ORDER$)\n\
    advectSemiLagrange(flags=flags, vel=vel, grid=color_b, order=$ADVECT_ORDER$)\n\
  \n\
  if using_fire:\n\
    print ('Advecting fire')\n\
    advectSemiLagrange(flags=flags, vel=vel, grid=fuel, order=$ADVECT_ORDER$)\n\
    advectSemiLagrange(flags=flags, vel=vel, grid=react, order=$ADVECT_ORDER$)\n\
  \n\
  print('Advecting density')\n\
  advectSemiLagrange(flags=flags, vel=vel, grid=density, order=$ADVECT_ORDER$)\n\
  \n\
  print('Advecting velocity')\n\
  advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=$ADVECT_ORDER$, openBounds=doOpen)\n\
  \n\
  for i in range(uvs):\n\
    print('Advecting UV and updating UVWeight')\n\
    advectSemiLagrange(flags=flags, vel=vel, grid=uv[i], order=$ADVECT_ORDER$)\n\
    updateUvWeight(resetTime=16.5 , index=i, numUvs=uvs, uv=uv[i])\n\
  \n\
  print('Walls')\n\
  setWallBcs(flags=flags, vel=vel)\n\
  \n\
  if using_heat:\n\
    print ('Adding heat buoyancy')\n\
    gravity=vec3(0,0,-0.0981) if dim==3 else vec3(0,-0.0981,0)\n\
    addBuoyancy2(flags=flags, grid=density, vel=vel, gravity=gravity, coefficient=$ALPHA$)\n\
    addBuoyancy2(flags=flags, grid=heat, vel=vel, gravity=gravity, coefficient=$BETA$*10)\n\
  else:\n\
    print ('Adding buoyancy')\n\
    gravity=vec3(0,0,-0.01 * $ALPHA$) if dim==3 else vec3(0,-0.01* $ALPHA$,0)\n\
    addBuoyancy(density=density, vel=vel, gravity=gravity, flags=flags)\n\
  \n\
  print('Vorticity')\n\
  if vorticity > 0.01:\n\
    vorticityConfinement( vel=vel, flags=flags, strength=$VORTICITY$ )\n\
  # TODO: print('Forcefield')\n\
  # TODO: addForceField(flags=flags, vel=vel, force=forces)\n\
  # TODO: forces.clear()\n\
  \n\
  print('Pressure')\n\
  solvePressure(flags=flags, vel=vel, pressure=pressure)\n\
  \n\
  print('Walls')\n\
  setWallBcs(flags=flags, vel=vel)\n\
  \n\
  print('Energy')\n\
  computeEnergy(flags=flags, vel=vel, energy=energy)\n\
  \n\
  copyMacToReal(source=vel, targetX=x_vel, targetY=y_vel, targetZ=z_vel)\n\
  s.step()\n\
\n\
def process_burn_low():\n\
  print('Process burn low')\n\
  if (using_colors):\n\
    processBurn(fuel=fuel, density=density, react=react, red=color_r, green=color_g, blue=color_b, heat=heat, burningRate=$BURNING_RATE$, flameSmoke=$FLAME_SMOKE$, ignitionTemp=$IGNITION_TEMP$, maxTemp=$MAX_TEMP$, dt=$DT$, flameSmokeColor=vec3($FLAME_SMOKE_COLOR_X$,$FLAME_SMOKE_COLOR_Y$,$FLAME_SMOKE_COLOR_Z$))\n\
  else:\n\
    processBurn(fuel=fuel, density=density, react=react, heat=heat, burningRate=$BURNING_RATE$, flameSmoke=$FLAME_SMOKE$, ignitionTemp=$IGNITION_TEMP$, maxTemp=$MAX_TEMP$, dt=$DT$, flameSmokeColor=vec3($FLAME_SMOKE_COLOR_X$,$FLAME_SMOKE_COLOR_Y$,$FLAME_SMOKE_COLOR_Z$))\n\
\n\
def update_flame_low():\n\
  print('Update flame low')\n\
  updateFlame(react=react, flame=flame)\n";

//////////////////////////////////////////////////////////////////////
// STEP FUNCTIONS HIGH
//////////////////////////////////////////////////////////////////////

const string smoke_step_high = "\n\
def step_high():\n\
  print('Step high')\n\
  interpolateMACGrid(source=vel, target=xl_vel)\n\
  sStr = 1.0 * wltStrength\n\
  sPos = 2.0\n\
  \n\
  print('Octaves')\n\
  for o in range(octaves):\n\
    for i in range(uvs):\n\
      uvWeight = getUvWeight(uv[i])\n\
      applyNoiseVec3(flags=xl_flags, target=xl_vel, noise=xl_wltnoise, scale=sStr * uvWeight, scaleSpatial=sPos , weight=energy, uv=uv[i])\n\
    sStr *= 0.06 # magic kolmogorov factor \n\
    sPos *= 2.0 \n\
  \n\
  for substep in range(upres):\n\
    if using_colors: \n\
      print ('Advecting colors high')\n\
      advectSemiLagrange(flags=xl_flags, vel=xl_vel, grid=xl_color_r, order=$ADVECT_ORDER$)\n\
      advectSemiLagrange(flags=xl_flags, vel=xl_vel, grid=xl_color_g, order=$ADVECT_ORDER$)\n\
      advectSemiLagrange(flags=xl_flags, vel=xl_vel, grid=xl_color_b, order=$ADVECT_ORDER$)\n\
    \n\
    if using_fire: \n\
      print ('Advecting fire high')\n\
      advectSemiLagrange(flags=xl_flags, vel=xl_vel, grid=xl_fuel, order=$ADVECT_ORDER$)\n\
      advectSemiLagrange(flags=xl_flags, vel=xl_vel, grid=xl_react, order=$ADVECT_ORDER$)\n\
    \n\
    print('Advecting density high')\n\
    advectSemiLagrange(flags=xl_flags, vel=xl_vel, grid=xl_density, order=$ADVECT_ORDER$)\n\
  \n\
  xl.step()\n\
\n\
def process_burn_high():\n\
  print('Process burn high')\n\
  if (using_colors):\n\
    processBurn(fuel=xl_fuel, density=xl_density, react=xl_react, red=xl_color_r, green=xl_color_g, blue=xl_color_b, burningRate=$BURNING_RATE$, flameSmoke=$FLAME_SMOKE$, ignitionTemp=$IGNITION_TEMP$, maxTemp=$MAX_TEMP$, dt=$DT$, flameSmokeColor=vec3($FLAME_SMOKE_COLOR_X$,$FLAME_SMOKE_COLOR_Y$,$FLAME_SMOKE_COLOR_Z$))\n\
  else:\n\
    processBurn(fuel=xl_fuel, density=xl_density, react=xl_react, burningRate=$BURNING_RATE$, flameSmoke=$FLAME_SMOKE$, ignitionTemp=$IGNITION_TEMP$, maxTemp=$MAX_TEMP$, dt=$DT$, flameSmokeColor=vec3($FLAME_SMOKE_COLOR_X$,$FLAME_SMOKE_COLOR_Y$,$FLAME_SMOKE_COLOR_Z$))\n\
\n\
def update_flame_high():\n\
  print('Update flame high')\n\
  updateFlame(react=xl_react, flame=xl_flame)\n";

//////////////////////////////////////////////////////////////////////
// STEP FUNCTIONS LIQUID
//////////////////////////////////////////////////////////////////////

const string liquid_step_low = "\n\
def sim_step_low(t):\n\
#update flags from density on first step\n\
  setWallBcs(flags=flags, vel=vel)\n\
  density.multConst(-1.)\n\
  print(using_colors)\n\
  global low_flags_updated\n\
  if not low_flags_updated:\n\
    print('Updating Flags from Levelset on startup!')\n\
    flags.updateFromLevelset(density)\n\
  low_flags_updated = True \n\
  setWallBcs(flags=flags, vel=vel)\n\
  density.reinitMarching(flags=flags, velTransport=vel)\n\
  advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2)\n\
  flags.updateFromLevelset(density)\n\
  \n\
  advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2)\n\
  addGravity(flags=flags, vel=vel, gravity=vec3(0,0,-0.981))\n\
  \n\
  # print current maximal velocity\n\
  maxvel = vel.getMaxValue()\n\
  print('Current max velocity %f ' % maxvel)\n\
  \n\
  # pressure solve\n\
  setWallBcs(flags=flags, vel=vel)\n\
  solvePressure(flags=flags, vel=vel, pressure=pressure, cgMaxIterFac=0.5, useResNorm=True) \n\
  setWallBcs(flags=flags, vel=vel)\n\
  s.step()\n\
  density.multConst(-1.)\n";

//////////////////////////////////////////////////////////////////////
// EXPORT GRIDS
//////////////////////////////////////////////////////////////////////

const string smoke_export_low = "\n\
import os\n\
print('Exporting grids')\n\
density.save(os.path.join('$MANTA_EXPORT_PATH$','density.uni'))\n\
flags.save(os.path.join('$MANTA_EXPORT_PATH$','flags.uni'))\n\
forces.save(os.path.join('$MANTA_EXPORT_PATH$','forces.uni'))\n\
inflow_grid.save(os.path.join('$MANTA_EXPORT_PATH$','inflow_low.uni'))\n\
fuel_inflow.save(os.path.join('$MANTA_EXPORT_PATH$','fuel_inflow.uni'))\n\
if using_colors:\n\
  color_r.save(os.path.join('$MANTA_EXPORT_PATH$','color_r.uni'))\n\
  color_g.save(os.path.join('$MANTA_EXPORT_PATH$','color_g.uni'))\n\
  color_b.save(os.path.join('$MANTA_EXPORT_PATH$','color_b.uni'))\n\
if using_heat:\n\
  heat.save(os.path.join('$MANTA_EXPORT_PATH$','heat.uni'))\n\
if using_fire:\n\
  flame.save(os.path.join('$MANTA_EXPORT_PATH$','flame.uni'))\n\
  fuel.save(os.path.join('$MANTA_EXPORT_PATH$','fuel.uni'))\n\
  react.save(os.path.join('$MANTA_EXPORT_PATH$','react.uni'))\n";

const string smoke_export_high = "\n\
print('Exporting grids')\n\
xl_density.save(os.path.join('$MANTA_EXPORT_PATH$','xl_density.uni'))\n\
xl_flags.save(os.path.join('$MANTA_EXPORT_PATH$','xl_flags.uni'))\n\
if using_colors:\n\
  xl_color_r.save(os.path.join('$MANTA_EXPORT_PATH$','xl_color_r.uni'))\n\
  xl_color_g.save(os.path.join('$MANTA_EXPORT_PATH$','xl_color_g.uni'))\n\
  xl_color_b.save(os.path.join('$MANTA_EXPORT_PATH$','xl_color_b.uni'))\n\
if using_fire:\n\
  xl_flame.save(os.path.join('$MANTA_EXPORT_PATH$','xl_flame.uni'))\n\
  xl_fuel.save(os.path.join('$MANTA_EXPORT_PATH$','xl_fuel.uni'))\n\
  xl_react.save(os.path.join('$MANTA_EXPORT_PATH$','xl_react.uni'))\n";

//////////////////////////////////////////////////////////////////////
// IMPORT GRIDS
//////////////////////////////////////////////////////////////////////

const string smoke_import_low = "\n\
def import_grids_low():\n\
  print('Importing grids')\n\
  density.load('$MANTA_EXPORT_PATH$density.uni')\n\
  flags.load('$MANTA_EXPORT_PATH$flags.uni')\n\
  forces.load('$MANTA_EXPORT_PATH$forces.uni')\n\
  inflow_grid.load('$MANTA_EXPORT_PATH$inflow_low.uni')\n\
  fuel_inflow.load('$MANTA_EXPORT_PATH$fuel_inflow.uni')\n\
  \n\
  if using_colors:\n\
    color_r.load('$MANTA_EXPORT_PATH$color_r.uni')\n\
    color_g.load('$MANTA_EXPORT_PATH$color_g.uni')\n\
    color_b.load('$MANTA_EXPORT_PATH$color_b.uni')\n\
  \n\
  if using_heat:\n\
    heat.load('$MANTA_EXPORT_PATH$heat.uni')\n\
  \n\
  if using_fire:\n\
    flame.load('$MANTA_EXPORT_PATH$flame.uni')\n\
    fuel.load('$MANTA_EXPORT_PATH$fuel.uni')\n\
    react.load('$MANTA_EXPORT_PATH$react.uni')\n";

const string smoke_import_high = "\n\
def import_grids_high():\n\
  print('Importing grids')\n\
  vel.load('$MANTA_EXPORT_PATH$vel.uni')\n\
  xl_density.load('$MANTA_EXPORT_PATH$xl_density.uni')\n\
  xl_flags.load('$MANTA_EXPORT_PATH$xl_flags.uni')\n\
  if using_colors:\n\
    xl_color_r.load('$MANTA_EXPORT_PATH$xl_color_r.uni')\n\
    xl_color_g.load('$MANTA_EXPORT_PATH$xl_color_g.uni')\n\
    xl_color_b.load('$MANTA_EXPORT_PATH$xl_color_b.uni')\n\
  if using_fire:\n\
    xl_flame.load('$MANTA_EXPORT_PATH$xl_flame.uni')\n\
    xl_fuel.load('$MANTA_EXPORT_PATH$xl_fuel.uni')\n\
    xl_react.load('$MANTA_EXPORT_PATH$xl_react.uni')\n";

//////////////////////////////////////////////////////////////////////
// INFLOW
//////////////////////////////////////////////////////////////////////

const string smoke_inflow_low = "\n\
def apply_inflow():\n\
  print('Applying inflow')\n\
  #inflow_grid.multConst(0.1)\n\
  #fuel_inflow.multConst(0.1)\n\
  density.add(inflow_grid)\n\
  fuel.add(fuel_inflow)\n";

const string smoke_inflow_high = "\n\
 # TODO\n";
