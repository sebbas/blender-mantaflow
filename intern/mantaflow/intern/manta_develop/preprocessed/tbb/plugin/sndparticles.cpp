

// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).

/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2016 Sebastian Barschkis, Nils Thuerey
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL)
 * http://www.gnu.org/licenses
 *
 * Secondary particle modeling plugin
 *
 ******************************************************************************/

#include "grid.h"
#include "levelset.h"
#include "particle.h"

using namespace std;
namespace Manta {

//! helper to calculate particle radius factor to cover the diagonal of a cell in 2d/3d
inline Real calculateRadiusFactor(Grid<Real> &grid, Real factor)
{
  return (grid.is3D() ? sqrt(3.) : sqrt(2.)) * (factor + .01);
}

//! adjust number of snd particles. optional parameters for life (0 = infinite life allowed)

void adjustSndParts(BasicParticleSystem &parts,
                    FlagGrid &flags,
                    LevelsetGrid &phi,
                    ParticleDataImpl<Vec3> &partVel,
                    ParticleDataImpl<Real> *partLife = NULL,
                    int maxDroplet = 16,
                    int maxBubble = 16,
                    int maxFloater = 16,
                    int maxTracer = 16)
{
  const Real dt = flags.getParent()->getDt();
  Real radiusFactor = 1.;
  const Real DROP_THRESH = calculateRadiusFactor(phi, radiusFactor);  // cell diagonal
  const Real FLOAT_THRESH = calculateRadiusFactor(phi, radiusFactor);

  Grid<int> numDroplet(flags.getParent());
  Grid<int> numBubble(flags.getParent());
  Grid<int> numFloater(flags.getParent());
  Grid<int> numTracer(flags.getParent());

  // Delete invalid particles. Then set new position to those that survived
  for (IndexInt idx = 0; idx < (int)parts.size(); idx++) {
    if (!parts.isActive(idx))
      continue;

    Vec3 p1 = parts.getPos(idx);
    Vec3i p1i = toVec3i(p1);

    // Kill parts at invalid positions
    if (p1i.x < 0 || p1i.y < 0 || p1i.z < 0) {
      parts.kill(idx);
      continue;
    }
    if (p1i.x > flags.getSizeX() - 1 || p1i.y > flags.getSizeY() - 1 ||
        p1i.z > flags.getSizeZ() - 1) {
      parts.kill(idx);
      continue;
    }

    // Try to save float / tracer particle by pushing it into the valid region
    Real phiv = phi.getInterpolated(parts.getPos(idx));
    if ((parts.getStatus(idx) & ParticleBase::PFOAM &&
         (phiv > FLOAT_THRESH || phiv < -FLOAT_THRESH)) ||
        (parts.getStatus(idx) & ParticleBase::PTRACER && phiv > 0.)) {
      Vec3 grad = getGradient(phi, p1i.x, p1i.y, p1i.z);
      if (normalize(grad) > VECTOR_EPSILON) {
        int direction = (phiv > 0.) ? -1 : 1;
        parts.setPos(idx, parts.getPos(idx) + direction * grad);

        // Update values for new position
        p1 = parts.getPos(idx);
        p1i = toVec3i(p1);
        phiv = phi.getInterpolated(parts.getPos(idx));
      }
    }

    // Next particle position (euler step)
    Vec3 p2 = parts.getPos(idx) + partVel[idx] * dt;
    Vec3i p2i = toVec3i(p2);

    // Kill particles depending on type. Especially those that were not converted to other particle
    // type
    if (parts.isSpray(idx) && phiv < -DROP_THRESH) {
      parts.kill(idx);
      continue;
    }
    if (parts.isBubble(idx) && phiv > 0.) {
      parts.kill(idx);
      continue;
    }
    if (parts.isFoam(idx) && (phiv > FLOAT_THRESH || phiv < -FLOAT_THRESH)) {
      parts.kill(idx);
      continue;
    }
    if (parts.isTracer(idx) && phiv > 0.) {
      parts.kill(idx);
      continue;
    }

    // Kill particles depending on current age
    if (partLife && (*partLife)[idx] <= 0.) {
      parts.kill(idx);
      continue;
    }

    // Kill particle if current or next position is invalid, ie outside or in obstacle
    if (!flags.isInBounds(p1i) || flags.isObstacle(p1i) || !flags.isInBounds(p2i) ||
        flags.isObstacle(p2i)) {
      parts.kill(idx);
      continue;
    }

    // Kill excess particles
    if (parts.isSpray(idx) && numDroplet(p1i) > maxDroplet) {
      parts.kill(idx);
      continue;
    }
    else {
      numDroplet(p1i) += 1;
    }
    if (parts.isBubble(idx) && numBubble(p1i) > maxBubble) {
      parts.kill(idx);
      continue;
    }
    else {
      numBubble(p1i) += 1;
    }
    if (parts.isFoam(idx) && numFloater(p1i) > maxFloater) {
      parts.kill(idx);
      continue;
    }
    else {
      numFloater(p1i) += 1;
    }
    if (parts.isTracer(idx) && numTracer(p1i) > maxTracer) {
      parts.kill(idx);
      continue;
    }
    else {
      numTracer(p1i) += 1;
    }
  }
  parts.doCompress();
}
static PyObject *_W_0(PyObject *_self, PyObject *_linargs, PyObject *_kwds)
{
  try {
    PbArgs _args(_linargs, _kwds);
    FluidSolver *parent = _args.obtainParent();
    bool noTiming = _args.getOpt<bool>("notiming", -1, 0);
    pbPreparePlugin(parent, "adjustSndParts", !noTiming);
    PyObject *_retval = 0;
    {
      ArgLocker _lock;
      BasicParticleSystem &parts = *_args.getPtr<BasicParticleSystem>("parts", 0, &_lock);
      FlagGrid &flags = *_args.getPtr<FlagGrid>("flags", 1, &_lock);
      LevelsetGrid &phi = *_args.getPtr<LevelsetGrid>("phi", 2, &_lock);
      ParticleDataImpl<Vec3> &partVel = *_args.getPtr<ParticleDataImpl<Vec3>>(
          "partVel", 3, &_lock);
      ParticleDataImpl<Real> *partLife = _args.getPtrOpt<ParticleDataImpl<Real>>(
          "partLife", 4, NULL, &_lock);
      int maxDroplet = _args.getOpt<int>("maxDroplet", 5, 16, &_lock);
      int maxBubble = _args.getOpt<int>("maxBubble", 6, 16, &_lock);
      int maxFloater = _args.getOpt<int>("maxFloater", 7, 16, &_lock);
      int maxTracer = _args.getOpt<int>("maxTracer", 8, 16, &_lock);
      _retval = getPyNone();
      adjustSndParts(
          parts, flags, phi, partVel, partLife, maxDroplet, maxBubble, maxFloater, maxTracer);
      _args.check();
    }
    pbFinalizePlugin(parent, "adjustSndParts", !noTiming);
    return _retval;
  }
  catch (std::exception &e) {
    pbSetError("adjustSndParts", e.what());
    return 0;
  }
}
static const Pb::Register _RP_adjustSndParts("", "adjustSndParts", _W_0);
extern "C" {
void PbRegister_adjustSndParts()
{
  KEEP_UNUSED(_RP_adjustSndParts);
}
}

//! update velocities. set new particle position. optional: convert between particle types,
//! partLife update

void updateSndParts(LevelsetGrid &phi,
                    FlagGrid &flags,
                    MACGrid &vel,
                    Vec3 gravity,
                    BasicParticleSystem &parts,
                    ParticleDataImpl<Vec3> &partVel,
                    ParticleDataImpl<Real> *partLife = NULL,
                    Real riseBubble = 0.5,
                    Real lifeDroplet = 30.0,
                    Real lifeBubble = 30.0,
                    Real lifeFloater = 30.0,
                    Real lifeTracer = 30.0)
{
  RandomStream mRand(9832);
  const Vec3 grav = gravity * flags.getParent()->getDt() / flags.getParent()->getDx();
  const Real dt = flags.getParent()->getDt();
  const Real framelength = flags.getParent()->getFrameLength();
  Real radiusFactor = 1.;
  const Real DROP_THRESH = 0.5f * calculateRadiusFactor(phi, radiusFactor);  // half cell diagonal

  for (IndexInt idx = 0; idx < (int)parts.size(); idx++) {
    if (!parts.isActive(idx))
      continue;

    // Update all already existing particles
    if ((parts.getStatus(idx) & ParticleBase::PNEW) == 0) {

      // Update particle velocity
      if (parts.isSpray(idx)) {
        partVel[idx] += grav;
      }
      else if (parts.isBubble(idx)) {
        Vec3 buoyancy = (-1) * grav * riseBubble;
        Vec3 randomVel = vel.getInterpolated(parts[idx].pos) * mRand.getFloat(0.25, 0.5);
        partVel[idx] += buoyancy + randomVel;
      }
      else if (parts.isFoam(idx) || parts.isTracer(idx)) {
        partVel[idx] = vel.getInterpolated(parts[idx].pos);
      }

      // Increase particle life
      if (partLife && (*partLife)[idx] > 0.f) {
        (*partLife)[idx] -= dt / framelength;
      }
    }

    // Update all new particles
    if (parts.getStatus(idx) & ParticleBase::PNEW) {

      // Init new particles (any type) with flow velocity
      partVel[idx] = vel.getInterpolated(parts[idx].pos);

      // Init particle life
      if (partLife) {
        if (parts.isSpray(idx))
          (*partLife)[idx] = lifeDroplet;
        if (parts.isBubble(idx))
          (*partLife)[idx] = lifeBubble;
        if (parts.isFoam(idx))
          (*partLife)[idx] = lifeFloater;
        if (parts.isTracer(idx))
          (*partLife)[idx] = lifeTracer;
      }
      // Make sure "new" flag gets removed
      parts.setStatus(idx, parts.getStatus(idx) & ~ParticleBase::PNEW);

      // New particle done here. Dont try to convert to other type
      continue;
    }

    // Set next particle position
    Vec3 pos = parts.getPos(idx) + partVel[idx] * dt;
    parts.setPos(idx, pos);

    // Get phiv for current and next particle position
    pos = parts.getPos(idx);
    Real phiv = phi.getInterpolated(pos);

    // Convert particle type
    if (parts.isSpray(idx) && phiv < -DROP_THRESH) {
      parts.setStatus(idx, ParticleBase::PNEW | ParticleBase::PBUBBLE);
    }
    else if (parts.isBubble(idx) && phiv > 0.) {
      parts.setStatus(idx, ParticleBase::PNEW | ParticleBase::PFOAM);
    }
  }
}
static PyObject *_W_1(PyObject *_self, PyObject *_linargs, PyObject *_kwds)
{
  try {
    PbArgs _args(_linargs, _kwds);
    FluidSolver *parent = _args.obtainParent();
    bool noTiming = _args.getOpt<bool>("notiming", -1, 0);
    pbPreparePlugin(parent, "updateSndParts", !noTiming);
    PyObject *_retval = 0;
    {
      ArgLocker _lock;
      LevelsetGrid &phi = *_args.getPtr<LevelsetGrid>("phi", 0, &_lock);
      FlagGrid &flags = *_args.getPtr<FlagGrid>("flags", 1, &_lock);
      MACGrid &vel = *_args.getPtr<MACGrid>("vel", 2, &_lock);
      Vec3 gravity = _args.get<Vec3>("gravity", 3, &_lock);
      BasicParticleSystem &parts = *_args.getPtr<BasicParticleSystem>("parts", 4, &_lock);
      ParticleDataImpl<Vec3> &partVel = *_args.getPtr<ParticleDataImpl<Vec3>>(
          "partVel", 5, &_lock);
      ParticleDataImpl<Real> *partLife = _args.getPtrOpt<ParticleDataImpl<Real>>(
          "partLife", 6, NULL, &_lock);
      Real riseBubble = _args.getOpt<Real>("riseBubble", 7, 0.5, &_lock);
      Real lifeDroplet = _args.getOpt<Real>("lifeDroplet", 8, 30.0, &_lock);
      Real lifeBubble = _args.getOpt<Real>("lifeBubble", 9, 30.0, &_lock);
      Real lifeFloater = _args.getOpt<Real>("lifeFloater", 10, 30.0, &_lock);
      Real lifeTracer = _args.getOpt<Real>("lifeTracer", 11, 30.0, &_lock);
      _retval = getPyNone();
      updateSndParts(phi,
                     flags,
                     vel,
                     gravity,
                     parts,
                     partVel,
                     partLife,
                     riseBubble,
                     lifeDroplet,
                     lifeBubble,
                     lifeFloater,
                     lifeTracer);
      _args.check();
    }
    pbFinalizePlugin(parent, "updateSndParts", !noTiming);
    return _retval;
  }
  catch (std::exception &e) {
    pbSetError("updateSndParts", e.what());
    return 0;
  }
}
static const Pb::Register _RP_updateSndParts("", "updateSndParts", _W_1);
extern "C" {
void PbRegister_updateSndParts()
{
  KEEP_UNUSED(_RP_updateSndParts);
}
}

//! sample new particles of given type. control amount of particles with amount and threshold
//! fields

void sampleSndParts(LevelsetGrid &phi,
                    LevelsetGrid &phiIn,
                    FlagGrid &flags,
                    MACGrid &vel,
                    BasicParticleSystem &parts,
                    int type,
                    Real amountDroplet,
                    Real amountFloater,
                    Real amountTracer,
                    Real thresholdDroplet)
{
  RandomStream mRand(9832);
  int a;
  Real radiusFactor = 1.;

  const Real DROP_THRESH = 0.5f * calculateRadiusFactor(phi, radiusFactor);  // half cell diagonal
  const Real FLOAT_THRESH = 0.5f * calculateRadiusFactor(phi, radiusFactor);

  // Split amount value into sampling steps (integral part) and probability (fractional part)
  float samplesDroplet, probDroplet, samplesFloater, probFloater, samplesTracer, probTracer;

  // Split amount variables into sample count and sample probability per cell
  float *amount, *samples, *probability;
  for (a = 0; a < 3; a++) {
    if (a == 0) {
      amount = &amountDroplet;
      samples = &samplesDroplet;
      probability = &probDroplet;
    }
    if (a == 1) {
      amount = &amountFloater;
      samples = &samplesFloater;
      probability = &probFloater;
    }
    if (a == 2) {
      amount = &amountTracer;
      samples = &samplesTracer;
      probability = &probTracer;
    }

    // Actual 'amount variable' splitting
    (*probability) = modf((double)(*amount), samples);
    (*probability) =
        ((*probability) == 0) ?
            1.0f :
            (*probability);  // e.g. map amount 1.0 to 100 percent probability (instead of 0.0)
    (*samples) = ceil(
        *amount);  // e.g. 0.1 amount samples once, 1.0 as well, 1.1 samples twice, ...
  }

  FOR_IJK_BND(phi, 0)
  {
    if (flags.isObstacle(i, j, k))
      continue;
    if (!flags.isFluid(i, j, k) && !flags.isEmpty(i, j, k))
      continue;

    const Vec3 pos = Vec3(i, j, k);

    // Droplets sampling
    for (a = 0; a < samplesDroplet; ++a) {

      // Surrounding fluid vel fast enough to generate particle?
      if (fabs(vel(i, j, k).x) < thresholdDroplet && fabs(vel(i, j, k).y) < thresholdDroplet &&
          fabs(vel(i, j, k).z) < thresholdDroplet)
        continue;

      // Only seed if random num exceeds given amount probability
      if (mRand.getFloat(0., 1.) > probDroplet)
        continue;

      if (type & ParticleBase::PSPRAY) {
        // Only generate drop particles at surface
        if (phi(i, j, k) < -DROP_THRESH || phi(i, j, k) > 0.)
          continue;

        // Only generate drops in convex regions
        Vec3 grad = getGradient(phi, i, j, k);
        Vec3 velC = vel.getCentered(i, j, k);
        if (dot(getNormalized(grad), getNormalized(velC)) < 0.75)
          continue;

        parts.addBuffered(pos + mRand.getVec3(), ParticleBase::PSPRAY);
      }
    }

    // Floater sampling
    for (a = 0; a < samplesFloater; ++a) {

      // Only seed if random num exceeds given amount probability
      if (mRand.getFloat(0., 1.) > probFloater)
        continue;

      if (type & ParticleBase::PFOAM) {
        // Only generate float particles at surface
        if (phiIn(i, j, k) < -FLOAT_THRESH || phiIn(i, j, k) > FLOAT_THRESH)
          continue;

        parts.addBuffered(pos + mRand.getVec3(), ParticleBase::PFOAM);
      }
    }

    // Tracer sampling
    for (a = 0; a < samplesTracer; ++a) {

      // Only seed if random num exceeds given amount probability
      if (mRand.getFloat(0., 1.) > probTracer)
        continue;

      if (type & ParticleBase::PTRACER) {
        // Only generate tracer particles inside fluid
        if (phiIn(i, j, k) > 0.)
          continue;

        parts.addBuffered(pos + mRand.getVec3(), ParticleBase::PTRACER);
      }
    }
  }
  // Insert buffered particles into particle system now.
  parts.insertBufferedParticles();
}
static PyObject *_W_2(PyObject *_self, PyObject *_linargs, PyObject *_kwds)
{
  try {
    PbArgs _args(_linargs, _kwds);
    FluidSolver *parent = _args.obtainParent();
    bool noTiming = _args.getOpt<bool>("notiming", -1, 0);
    pbPreparePlugin(parent, "sampleSndParts", !noTiming);
    PyObject *_retval = 0;
    {
      ArgLocker _lock;
      LevelsetGrid &phi = *_args.getPtr<LevelsetGrid>("phi", 0, &_lock);
      LevelsetGrid &phiIn = *_args.getPtr<LevelsetGrid>("phiIn", 1, &_lock);
      FlagGrid &flags = *_args.getPtr<FlagGrid>("flags", 2, &_lock);
      MACGrid &vel = *_args.getPtr<MACGrid>("vel", 3, &_lock);
      BasicParticleSystem &parts = *_args.getPtr<BasicParticleSystem>("parts", 4, &_lock);
      int type = _args.get<int>("type", 5, &_lock);
      Real amountDroplet = _args.get<Real>("amountDroplet", 6, &_lock);
      Real amountFloater = _args.get<Real>("amountFloater", 7, &_lock);
      Real amountTracer = _args.get<Real>("amountTracer", 8, &_lock);
      Real thresholdDroplet = _args.get<Real>("thresholdDroplet", 9, &_lock);
      _retval = getPyNone();
      sampleSndParts(phi,
                     phiIn,
                     flags,
                     vel,
                     parts,
                     type,
                     amountDroplet,
                     amountFloater,
                     amountTracer,
                     thresholdDroplet);
      _args.check();
    }
    pbFinalizePlugin(parent, "sampleSndParts", !noTiming);
    return _retval;
  }
  catch (std::exception &e) {
    pbSetError("sampleSndParts", e.what());
    return 0;
  }
}
static const Pb::Register _RP_sampleSndParts("", "sampleSndParts", _W_2);
extern "C" {
void PbRegister_sampleSndParts()
{
  KEEP_UNUSED(_RP_sampleSndParts);
}
}

}  // namespace Manta
