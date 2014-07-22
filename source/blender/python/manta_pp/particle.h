




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Base class for particle systems
 *
 ******************************************************************************/

#ifndef _PARTICLE_H
#define _PARTICLE_H

#include <vector>
#include "grid.h"
#include "vectorbase.h"
#include "integrator.h"
#include "randomstream.h"
namespace Manta {

// fwd decl
template<class T> class Grid;
class ParticleDataBase;
template<class T> class ParticleDataImpl;

//! Baseclass for particle systems. Does not implement any data
class ParticleBase : public PbClass {public:
	enum SystemType { BASE=0, PARTICLE, VORTEX, FILAMENT, FLIP, TURBULENCE, INDEX };
	
	enum ParticleStatus {
		PNONE         = 0,
		PNEW          = (1<<1),  // particles newly created in this step
		PDELETE       = (1<<10), // mark as deleted, will be deleted in next compress() step
		PINVALID      = (1<<30), // unused
	};

	ParticleBase(FluidSolver* parent); static int _W_0 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = Pb::objFromPy(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "ParticleBase::ParticleBase" ); { ArgLocker _lock; FluidSolver* parent = _args.getPtr<FluidSolver >("parent",0,&_lock);  obj = new ParticleBase(parent); obj->registerObject(_self, &_args); _args.check(); } pbFinalizePlugin(obj->getParent(),"ParticleBase::ParticleBase" ); return 0; } catch(std::exception& e) { pbSetError("ParticleBase::ParticleBase",e.what()); return -1; } }
	virtual ~ParticleBase();

	//! copy all the particle data thats registered with the other particle system to this one
	virtual void cloneParticleData(ParticleBase* nm);

	virtual SystemType getType() const { return BASE; }
	virtual std::string infoString() const; 
	virtual ParticleBase* clone() { assertMsg( false , "Dont use, override..."); return NULL; } 

	// slow virtual function to query size, do not use in kernels! use size() instead
	virtual int getSizeSlow() const { assertMsg( false , "Dont use, override..."); return 0; } 

	//! add a position as potential candidate for new particle (todo, make usable from parallel threads)
	inline void addBuffered(const Vec3& pos);

	//! debug info about pdata
	std::string debugInfoPdata();

	// particle data functions

	//! create a particle data object
	PbClass* create(PbType type, PbTypeVec T=PbTypeVec(), const std::string& name = ""); static PyObject* _W_1 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); ParticleBase* pbo = dynamic_cast<ParticleBase*>(Pb::objFromPy(_self)); pbPreparePlugin(pbo->getParent(), "ParticleBase::create"); PyObject *_retval = 0; { ArgLocker _lock; PbType type = _args.get<PbType >("type",0,&_lock); PbTypeVec T = _args.getOpt<PbTypeVec >("T",1,PbTypeVec(),&_lock); const std::string& name = _args.getOpt<std::string >("name",2,"",&_lock);  pbo->_args.copy(_args);  _retval = toPy(pbo->create(type,T,name));  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"ParticleBase::create"); return _retval; } catch(std::exception& e) { pbSetError("ParticleBase::create",e.what()); return 0; } }
	//! add a particle data field, set its parent particle-system pointer
	void registerPdata(ParticleDataBase* pdata);
	void registerPdataReal(ParticleDataImpl<Real>* pdata);
	void registerPdataVec3(ParticleDataImpl<Vec3>* pdata);
	void registerPdataInt (ParticleDataImpl<int >* pdata);
	//! remove a particle data entry
	void deregister(ParticleDataBase* pdata);
	//! add one zero entry to all data fields
	void addAllPdata();
	// note - deletion of pdata is handled in compress function

	//! how many are there?
	int getNumPdata() const { return mPartData.size(); }
	//! access one of the fields
	ParticleDataBase* getPdata(int i) { return mPartData[i]; }

protected:  
	//! new particle candidates
	std::vector<Vec3> mNewBuffer;

	//! allow automatic compression / resize? disallowed for, eg, flip particle systems
	bool mAllowCompress;

	//! store particle data , each pointer has its own storage vector of a certain type (int, real, vec3)
	std::vector<ParticleDataBase*> mPartData;
	//! lists of different types, for fast operations w/o virtual function calls (all calls necessary per particle)
	std::vector< ParticleDataImpl<Real> *> mPdataReal;
	std::vector< ParticleDataImpl<Vec3> *> mPdataVec3;
	std::vector< ParticleDataImpl<int> *>  mPdataInt; 	//! indicate that pdata of this particle system is copied, and needs to be freed
	bool mFreePdata; public: PbArgs _args;}
#define _C_ParticleBase
;


//! Main class for particle systems
/*! Basetype S must at least contain flag, pos fields */
template<class S> class ParticleSystem : public ParticleBase {public:    
	ParticleSystem(FluidSolver* parent) :ParticleBase(parent),mDeletes(0),mDeleteChunk(0){} static int _W_2 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = Pb::objFromPy(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "ParticleSystem::ParticleSystem" ); { ArgLocker _lock; FluidSolver* parent = _args.getPtr<FluidSolver >("parent",0,&_lock);  obj = new ParticleSystem(parent); obj->registerObject(_self, &_args); _args.check(); } pbFinalizePlugin(obj->getParent(),"ParticleSystem::ParticleSystem" ); return 0; } catch(std::exception& e) { pbSetError("ParticleSystem::ParticleSystem",e.what()); return -1; } }
	virtual ~ParticleSystem() {};
	
	virtual SystemType getType() const { return S::getType(); };
	
	// accessors
	inline S& operator[](int idx)             { DEBUG_ONLY(checkPartIndex(idx)); return mData[idx]; }
	inline const S& operator[](int idx) const { DEBUG_ONLY(checkPartIndex(idx)); return mData[idx]; }
	// return size of container
	inline int size() const { return mData.size(); } static PyObject* _W_3 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); ParticleSystem* pbo = dynamic_cast<ParticleSystem*>(Pb::objFromPy(_self)); pbPreparePlugin(pbo->getParent(), "ParticleSystem::size"); PyObject *_retval = 0; { ArgLocker _lock;  pbo->_args.copy(_args);  _retval = toPy(pbo->size());  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"ParticleSystem::size"); return _retval; } catch(std::exception& e) { pbSetError("ParticleSystem::size",e.what()); return 0; } }
	// slow virtual function of base class, also returns size
	virtual int getSizeSlow() const { return size(); }

	// query status
	inline int  getStatus(int idx) { DEBUG_ONLY(checkPartIndex(idx)); return mData[idx].flag; }
	inline bool isActive(int idx)  { DEBUG_ONLY(checkPartIndex(idx)); return (mData[idx].flag & PDELETE) == 0; }
	
	//! safe accessor for python
	void setPos(int idx, const Vec3& pos) { DEBUG_ONLY(checkPartIndex(idx)); mData[idx].pos = pos; } static PyObject* _W_4 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); ParticleSystem* pbo = dynamic_cast<ParticleSystem*>(Pb::objFromPy(_self)); pbPreparePlugin(pbo->getParent(), "ParticleSystem::setPos"); PyObject *_retval = 0; { ArgLocker _lock; int idx = _args.get<int >("idx",0,&_lock); const Vec3& pos = _args.get<Vec3 >("pos",1,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->setPos(idx,pos);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"ParticleSystem::setPos"); return _retval; } catch(std::exception& e) { pbSetError("ParticleSystem::setPos",e.what()); return 0; } }
	Vec3 getPos(int idx) { DEBUG_ONLY(checkPartIndex(idx)); return mData[idx].pos; } static PyObject* _W_5 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); ParticleSystem* pbo = dynamic_cast<ParticleSystem*>(Pb::objFromPy(_self)); pbPreparePlugin(pbo->getParent(), "ParticleSystem::getPos"); PyObject *_retval = 0; { ArgLocker _lock; int idx = _args.get<int >("idx",0,&_lock);  pbo->_args.copy(_args);  _retval = toPy(pbo->getPos(idx));  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"ParticleSystem::getPos"); return _retval; } catch(std::exception& e) { pbSetError("ParticleSystem::getPos",e.what()); return 0; } }
	//! copy all positions into pdata vec3 field
	void getPosPdata(ParticleDataImpl<Vec3>& target); static PyObject* _W_6 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); ParticleSystem* pbo = dynamic_cast<ParticleSystem*>(Pb::objFromPy(_self)); pbPreparePlugin(pbo->getParent(), "ParticleSystem::getPosPdata"); PyObject *_retval = 0; { ArgLocker _lock; ParticleDataImpl<Vec3>& target = *_args.getPtr<ParticleDataImpl<Vec3> >("target",0,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->getPosPdata(target);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"ParticleSystem::getPosPdata"); return _retval; } catch(std::exception& e) { pbSetError("ParticleSystem::getPosPdata",e.what()); return 0; } }
	void setPosPdata(ParticleDataImpl<Vec3>& source); static PyObject* _W_7 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); ParticleSystem* pbo = dynamic_cast<ParticleSystem*>(Pb::objFromPy(_self)); pbPreparePlugin(pbo->getParent(), "ParticleSystem::setPosPdata"); PyObject *_retval = 0; { ArgLocker _lock; ParticleDataImpl<Vec3>& source = *_args.getPtr<ParticleDataImpl<Vec3> >("source",0,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->setPosPdata(source);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"ParticleSystem::setPosPdata"); return _retval; } catch(std::exception& e) { pbSetError("ParticleSystem::setPosPdata",e.what()); return 0; } }
	//! transform coordinate system from one grid size to another (usually upon load)
	void transformPositions( Vec3i dimOld, Vec3i dimNew );

	//! explicitly trigger compression from outside
	void doCompress() { if ( mDeletes > mDeleteChunk) compress(); }
	//! insert buffered positions as new particles, update additional particle data
	void insertBufferedParticles();
	//! resize data vector, and all pdata fields
	void resizeAll(int newsize);
	
	// adding and deleting 
	inline void kill(int idx);
	int add(const S& data);
	// remove all particles, init 0 length arrays (also pdata)
	void clear(); static PyObject* _W_8 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); ParticleSystem* pbo = dynamic_cast<ParticleSystem*>(Pb::objFromPy(_self)); pbPreparePlugin(pbo->getParent(), "ParticleSystem::clear"); PyObject *_retval = 0; { ArgLocker _lock;  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->clear();  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"ParticleSystem::clear"); return _retval; } catch(std::exception& e) { pbSetError("ParticleSystem::clear",e.what()); return 0; } }
			
	//! Advect particle in grid velocity field
	void advectInGrid(FlagGrid& flags, MACGrid& vel, int integrationMode, bool deleteInObstacle=true ); static PyObject* _W_9 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); ParticleSystem* pbo = dynamic_cast<ParticleSystem*>(Pb::objFromPy(_self)); pbPreparePlugin(pbo->getParent(), "ParticleSystem::advectInGrid"); PyObject *_retval = 0; { ArgLocker _lock; FlagGrid& flags = *_args.getPtr<FlagGrid >("flags",0,&_lock); MACGrid& vel = *_args.getPtr<MACGrid >("vel",1,&_lock); int integrationMode = _args.get<int >("integrationMode",2,&_lock); bool deleteInObstacle = _args.getOpt<bool >("deleteInObstacle",3,true ,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->advectInGrid(flags,vel,integrationMode,deleteInObstacle);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"ParticleSystem::advectInGrid"); return _retval; } catch(std::exception& e) { pbSetError("ParticleSystem::advectInGrid",e.what()); return 0; } }
	
	//! Project particles outside obstacles
	void projectOutside(Grid<Vec3>& gradient); static PyObject* _W_10 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); ParticleSystem* pbo = dynamic_cast<ParticleSystem*>(Pb::objFromPy(_self)); pbPreparePlugin(pbo->getParent(), "ParticleSystem::projectOutside"); PyObject *_retval = 0; { ArgLocker _lock; Grid<Vec3>& gradient = *_args.getPtr<Grid<Vec3> >("gradient",0,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->projectOutside(gradient);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"ParticleSystem::projectOutside"); return _retval; } catch(std::exception& e) { pbSetError("ParticleSystem::projectOutside",e.what()); return 0; } }
	
	virtual ParticleBase* clone();
	virtual std::string infoString() const;

	//! debugging
	inline void checkPartIndex(int idx) const;
	
protected:  
	//! deletion count , and interval for re-compressing 
	int mDeletes, mDeleteChunk;    
	//! the particle data
	std::vector<S> mData;    
 	//! reduce storage , called by doCompress
	virtual void compress();  public: PbArgs _args;}
#define _C_ParticleSystem
;

//******************************************************************************

//! Simplest data class for particle systems
struct BasicParticleData {
public:
	BasicParticleData() : pos(0.), flag(0) {}
	BasicParticleData(const Vec3& p) : pos(p), flag(0) {}
	static ParticleBase::SystemType getType() { return ParticleBase::PARTICLE; }

	//! data
	Vec3 pos;
	int  flag;
};

class BasicParticleSystem : public ParticleSystem<BasicParticleData> {public:
	BasicParticleSystem(FluidSolver* parent); static int _W_11 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = Pb::objFromPy(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "BasicParticleSystem::BasicParticleSystem" ); { ArgLocker _lock; FluidSolver* parent = _args.getPtr<FluidSolver >("parent",0,&_lock);  obj = new BasicParticleSystem(parent); obj->registerObject(_self, &_args); _args.check(); } pbFinalizePlugin(obj->getParent(),"BasicParticleSystem::BasicParticleSystem" ); return 0; } catch(std::exception& e) { pbSetError("BasicParticleSystem::BasicParticleSystem",e.what()); return -1; } }
	
	//! file io
	void save(std::string name); static PyObject* _W_12 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); BasicParticleSystem* pbo = dynamic_cast<BasicParticleSystem*>(Pb::objFromPy(_self)); pbPreparePlugin(pbo->getParent(), "BasicParticleSystem::save"); PyObject *_retval = 0; { ArgLocker _lock; std::string name = _args.get<std::string >("name",0,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->save(name);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"BasicParticleSystem::save"); return _retval; } catch(std::exception& e) { pbSetError("BasicParticleSystem::save",e.what()); return 0; } }
	void load(std::string name); static PyObject* _W_13 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); BasicParticleSystem* pbo = dynamic_cast<BasicParticleSystem*>(Pb::objFromPy(_self)); pbPreparePlugin(pbo->getParent(), "BasicParticleSystem::load"); PyObject *_retval = 0; { ArgLocker _lock; std::string name = _args.get<std::string >("name",0,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->load(name);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"BasicParticleSystem::load"); return _retval; } catch(std::exception& e) { pbSetError("BasicParticleSystem::load",e.what()); return 0; } }

	// save to text file
	void writeParticlesText(std::string name);
	// other output formats
	void writeParticlesRawPositionsGz(std::string name);
	void writeParticlesRawVelocityGz(std::string name);

	// add particles in python
	void addParticle(Vec3 pos) { add(BasicParticleData(pos)); } static PyObject* _W_14 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); BasicParticleSystem* pbo = dynamic_cast<BasicParticleSystem*>(Pb::objFromPy(_self)); pbPreparePlugin(pbo->getParent(), "BasicParticleSystem::addParticle"); PyObject *_retval = 0; { ArgLocker _lock; Vec3 pos = _args.get<Vec3 >("pos",0,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->addParticle(pos);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"BasicParticleSystem::addParticle"); return _retval; } catch(std::exception& e) { pbSetError("BasicParticleSystem::addParticle",e.what()); return 0; } }
 	// dangerous, get low level access - avoid usage, only used in vortex filament advection for now
	std::vector<BasicParticleData>& getData() { return mData; } public: PbArgs _args;}
#define _C_BasicParticleSystem
;


//******************************************************************************

//! Index into other particle system
//  used for grid based neighborhood searches on generic particle systems (stores
//  only active particles, and reduces copied data)
//  note - pos & flag are disabled here, do not use!
struct ParticleIndexData {
public:
	ParticleIndexData() : sourceIndex(0) {}
	static ParticleBase::SystemType getType() { return ParticleBase::INDEX; }

	int  sourceIndex; // index of this particle in the original particle system
	// note - the following two are needed for template instantiation, but not used
	// for the particle index system (use values from original one!)
	static Vec3 pos;  // do not use... 
	static int  flag; // not needed usally 
	//Vec3 pos; // enable for debugging
};

class ParticleIndexSystem : public ParticleSystem<ParticleIndexData> {public:
	ParticleIndexSystem(FluidSolver* parent) :ParticleSystem<ParticleIndexData>(parent){} static int _W_15 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = Pb::objFromPy(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "ParticleIndexSystem::ParticleIndexSystem" ); { ArgLocker _lock; FluidSolver* parent = _args.getPtr<FluidSolver >("parent",0,&_lock);  obj = new ParticleIndexSystem(parent); obj->registerObject(_self, &_args); _args.check(); } pbFinalizePlugin(obj->getParent(),"ParticleIndexSystem::ParticleIndexSystem" ); return 0; } catch(std::exception& e) { pbSetError("ParticleIndexSystem::ParticleIndexSystem",e.what()); return -1; } };
	 	//! we only need a resize function...
	void resize(int size) { mData.resize(size); } public: PbArgs _args;}
#define _C_ParticleIndexSystem
;



//******************************************************************************

//! Particle set with connectivity

template<class DATA, class CON> class ConnectedParticleSystem : public ParticleSystem<DATA> {public:
	ConnectedParticleSystem(FluidSolver* parent) :ParticleSystem<DATA>(parent){} static int _W_16 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = Pb::objFromPy(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "ConnectedParticleSystem::ConnectedParticleSystem" ); { ArgLocker _lock; FluidSolver* parent = _args.getPtr<FluidSolver >("parent",0,&_lock);  obj = new ConnectedParticleSystem(parent); obj->registerObject(_self, &_args); _args.check(); } pbFinalizePlugin(obj->getParent(),"ConnectedParticleSystem::ConnectedParticleSystem" ); return 0; } catch(std::exception& e) { pbSetError("ConnectedParticleSystem::ConnectedParticleSystem",e.what()); return -1; } }
	
	// accessors
	inline bool isSegActive(int i) { return (mSegments[i].flag & ParticleBase::PDELETE) == 0; }    
	inline int segSize() const { return mSegments.size(); }    
	inline CON& seg(int i) { return mSegments[i]; }
	inline const CON& seg(int i) const { return mSegments[i]; }
		
	virtual ParticleBase* clone();
	
protected:
	std::vector<CON> mSegments; 	virtual void compress();     public: PbArgs _args;}
#define _C_ConnectedParticleSystem
;

//******************************************************************************

//! abstract interface for particle data
class ParticleDataBase : public PbClass {public:
	ParticleDataBase(FluidSolver* parent); static int _W_17 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = Pb::objFromPy(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "ParticleDataBase::ParticleDataBase" ); { ArgLocker _lock; FluidSolver* parent = _args.getPtr<FluidSolver >("parent",0,&_lock);  obj = new ParticleDataBase(parent); obj->registerObject(_self, &_args); _args.check(); } pbFinalizePlugin(obj->getParent(),"ParticleDataBase::ParticleDataBase" ); return 0; } catch(std::exception& e) { pbSetError("ParticleDataBase::ParticleDataBase",e.what()); return -1; } }
	virtual ~ParticleDataBase(); 

	enum PdataType { UNKNOWN=0, DATA_INT, DATA_REAL, DATA_VEC3 };

	// interface functions, using assert instead of pure virtual for python compatibility
	virtual int  size() const { assertMsg( false , "Dont use, override..."); return 0; } 
	virtual void add()        { assertMsg( false , "Dont use, override..."); return;   }
	virtual ParticleDataBase* clone() { assertMsg( false , "Dont use, override..."); return NULL; }
	virtual PdataType getType() const { assertMsg( false , "Dont use, override..."); return UNKNOWN; } 
	virtual void resize(int size)     { assertMsg( false , "Dont use, override..."); return;  }
	virtual void copyValueSlow(int from, int to) { assertMsg( false , "Dont use, override..."); return;  }

	//! set base pointer
	void setParticleSys(ParticleBase* set) { mpParticleSys = set; }

	//! debugging
	inline void checkPartIndex(int idx) const;

protected: 	ParticleBase* mpParticleSys; public: PbArgs _args;}
#define _C_ParticleDataBase
;


//! abstract interface for particle data

template<class T> class ParticleDataImpl : public ParticleDataBase {public:
	ParticleDataImpl(FluidSolver* parent); static int _W_18 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = Pb::objFromPy(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "ParticleDataImpl::ParticleDataImpl" ); { ArgLocker _lock; FluidSolver* parent = _args.getPtr<FluidSolver >("parent",0,&_lock);  obj = new ParticleDataImpl(parent); obj->registerObject(_self, &_args); _args.check(); } pbFinalizePlugin(obj->getParent(),"ParticleDataImpl::ParticleDataImpl" ); return 0; } catch(std::exception& e) { pbSetError("ParticleDataImpl::ParticleDataImpl",e.what()); return -1; } }
	ParticleDataImpl(FluidSolver* parent, ParticleDataImpl<T>* other);
	virtual ~ParticleDataImpl();

	//! access data
	inline T& get(int idx)            { DEBUG_ONLY(checkPartIndex(idx)); return mData[idx]; }
	inline const T get(int idx) const { DEBUG_ONLY(checkPartIndex(idx)); return mData[idx]; }
	inline T& operator[](int idx)            { DEBUG_ONLY(checkPartIndex(idx)); return mData[idx]; }
	inline const T operator[](int idx) const { DEBUG_ONLY(checkPartIndex(idx)); return mData[idx]; }

	// set all values to 0, note - different from particleSystem::clear! doesnt modify size of array (has to stay in sync with parent system)
	void clear(); static PyObject* _W_19 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); ParticleDataImpl* pbo = dynamic_cast<ParticleDataImpl*>(Pb::objFromPy(_self)); pbPreparePlugin(pbo->getParent(), "ParticleDataImpl::clear"); PyObject *_retval = 0; { ArgLocker _lock;  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->clear();  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"ParticleDataImpl::clear"); return _retval; } catch(std::exception& e) { pbSetError("ParticleDataImpl::clear",e.what()); return 0; } }

	//! set grid from which to get data...
	void setSource(Grid<T>* grid, bool isMAC=false ); static PyObject* _W_20 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); ParticleDataImpl* pbo = dynamic_cast<ParticleDataImpl*>(Pb::objFromPy(_self)); pbPreparePlugin(pbo->getParent(), "ParticleDataImpl::setSource"); PyObject *_retval = 0; { ArgLocker _lock; Grid<T>* grid = _args.getPtr<Grid<T> >("grid",0,&_lock); bool isMAC = _args.getOpt<bool >("isMAC",1,false ,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->setSource(grid,isMAC);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"ParticleDataImpl::setSource"); return _retval; } catch(std::exception& e) { pbSetError("ParticleDataImpl::setSource",e.what()); return 0; } }

	// particle data base interface
	virtual int  size() const;
	virtual void add();
	virtual ParticleDataBase* clone();
	virtual PdataType getType() const;
	virtual void resize(int s);
	virtual void copyValueSlow(int from, int to);

	// fast inlined functions for per particle operations
	inline void copyValue(int from, int to) { get(to) = get(from); } 
	void initNewValue(int idx, Vec3 pos);
	
	//! file io
	void save(std::string name); static PyObject* _W_21 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); ParticleDataImpl* pbo = dynamic_cast<ParticleDataImpl*>(Pb::objFromPy(_self)); pbPreparePlugin(pbo->getParent(), "ParticleDataImpl::save"); PyObject *_retval = 0; { ArgLocker _lock; std::string name = _args.get<std::string >("name",0,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->save(name);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"ParticleDataImpl::save"); return _retval; } catch(std::exception& e) { pbSetError("ParticleDataImpl::save",e.what()); return 0; } }
	void load(std::string name); static PyObject* _W_22 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); ParticleDataImpl* pbo = dynamic_cast<ParticleDataImpl*>(Pb::objFromPy(_self)); pbPreparePlugin(pbo->getParent(), "ParticleDataImpl::load"); PyObject *_retval = 0; { ArgLocker _lock; std::string name = _args.get<std::string >("name",0,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->load(name);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"ParticleDataImpl::load"); return _retval; } catch(std::exception& e) { pbSetError("ParticleDataImpl::load",e.what()); return 0; } }
protected:
	//! data storage
	std::vector<T> mData; 

	//! optionally , we might have an associated grid from which to grab new data
	Grid<T>* mpGridSource; 	//! unfortunately , we need to distinguish mac vs regular vec3
	bool mGridSourceMAC; public: PbArgs _args;}
#define _C_ParticleDataImpl
;






//******************************************************************************
// Implementation
//******************************************************************************

const int DELETE_PART = 20; // chunk size for compression

void ParticleBase::addBuffered(const Vec3& pos) {
	mNewBuffer.push_back(pos);
}
   
template<class S>
void ParticleSystem<S>::clear() {
	mDeleteChunk = mDeletes = 0;
	this->resizeAll(0); // instead of mData.clear
}

template<class S>
int ParticleSystem<S>::add(const S& data) {
	mData.push_back(data); 
	mDeleteChunk = mData.size() / DELETE_PART;
	this->addAllPdata();
	return mData.size()-1;
}

template<class S>
inline void ParticleSystem<S>::kill(int idx)     { 
	assertMsg(idx>=0 && idx<size(), "Index out of bounds");
	mData[idx].flag |= PDELETE; 
	if ( (++mDeletes > mDeleteChunk) && (mAllowCompress) ) compress(); 
}

template<class S>
void ParticleSystem<S>::getPosPdata(ParticleDataImpl<Vec3>& target) {
	for(int i=0; i<(int)this->size(); ++i) {
		target[i] = this->getPos(i);
	}
}
template<class S>
void ParticleSystem<S>::setPosPdata(ParticleDataImpl<Vec3>& target) {
	for(int i=0; i<(int)this->size(); ++i) {
		this->getPos(i) = target[i];
	}
}

template<class S>
void ParticleSystem<S>::transformPositions( Vec3i dimOld, Vec3i dimNew )
{
	Vec3 factor = calcGridSizeFactor( dimNew, dimOld );
	for(int i=0; i<(int)this->size(); ++i) {
		this->setPos(i, this->getPos(i) * factor );
	}
}

// check for deletion/invalid position, otherwise return velocity



template <class S>  struct GridAdvectKernel : public KernelBase { GridAdvectKernel(std::vector<S>& p, const MACGrid& vel, const FlagGrid& flags, Real dt, bool deleteInObstacle ) :  KernelBase(p.size()) ,p(p),vel(vel),flags(flags),dt(dt),deleteInObstacle(deleteInObstacle) ,u((size))  { run(); }  inline void op(int idx, std::vector<S>& p, const MACGrid& vel, const FlagGrid& flags, Real dt, bool deleteInObstacle  ,std::vector<Vec3> & u)  {
	if (p[idx].flag & ParticleBase::PDELETE) {
		u[idx] =_0;
	} else if (!flags.isInBounds(p[idx].pos,1) || flags.isObstacle(p[idx].pos)) {
		u[idx] = _0;

		// for simple tracer particles, its convenient to delete particles right away
		// for other sim types, eg flip, we can try to fix positions later on
		if(deleteInObstacle) 
			p[idx].flag |= ParticleBase::PDELETE;
	} else {
		u[idx] = vel.getInterpolated(p[idx].pos) * dt;
	}
}   inline operator std::vector<Vec3> () { return u; } inline std::vector<Vec3>  & getRet() { return u; }  inline std::vector<S>& getArg0() { return p; } typedef std::vector<S> type0;inline const MACGrid& getArg1() { return vel; } typedef MACGrid type1;inline const FlagGrid& getArg2() { return flags; } typedef FlagGrid type2;inline Real& getArg3() { return dt; } typedef Real type3;inline bool& getArg4() { return deleteInObstacle; } typedef bool type4; void run() {  const int _sz = size; for (int i=0; i < _sz; i++) op(i, p,vel,flags,dt,deleteInObstacle,u);  } std::vector<S>& p; const MACGrid& vel; const FlagGrid& flags; Real dt; bool deleteInObstacle;  std::vector<Vec3>  u;  };;

// final check after advection to make sure particles haven't escaped
// (similar to particle advection kernel)

template <class S>  struct KnDeleteInObstacle : public KernelBase { KnDeleteInObstacle(std::vector<S>& p, const FlagGrid& flags) :  KernelBase(p.size()) ,p(p),flags(flags)   { run(); }  inline void op(int idx, std::vector<S>& p, const FlagGrid& flags )  {
	if (p[idx].flag & ParticleBase::PDELETE) return;
	if (!flags.isInBounds(p[idx].pos,1) || flags.isObstacle(p[idx].pos)) {
		p[idx].flag |= ParticleBase::PDELETE;
	} 
}   inline std::vector<S>& getArg0() { return p; } typedef std::vector<S> type0;inline const FlagGrid& getArg1() { return flags; } typedef FlagGrid type1; void run() {  const int _sz = size; for (int i=0; i < _sz; i++) op(i, p,flags);  } std::vector<S>& p; const FlagGrid& flags;   };
// at least make sure all particles are inside domain

template <class S>  struct KnClampPositions : public KernelBase { KnClampPositions(std::vector<S>& p, const FlagGrid& flags) :  KernelBase(p.size()) ,p(p),flags(flags)   { run(); }  inline void op(int idx, std::vector<S>& p, const FlagGrid& flags )  {
	if (p[idx].flag & ParticleBase::PDELETE) return;
	if (!flags.isInBounds(p[idx].pos,0) ) {
		p[idx].pos = clamp( p[idx].pos, Vec3(0.), toVec3(flags.getSize())-Vec3(1.) );
	} 
}   inline std::vector<S>& getArg0() { return p; } typedef std::vector<S> type0;inline const FlagGrid& getArg1() { return flags; } typedef FlagGrid type1; void run() {  const int _sz = size; for (int i=0; i < _sz; i++) op(i, p,flags);  } std::vector<S>& p; const FlagGrid& flags;   };

// advection plugin
template<class S>
void ParticleSystem<S>::advectInGrid(FlagGrid& flags, MACGrid& vel, int integrationMode, bool deleteInObstacle ) {
	GridAdvectKernel<S> kernel(mData, vel, flags, getParent()->getDt(), deleteInObstacle );
	integratePointSet(kernel, integrationMode);
	if(deleteInObstacle) KnDeleteInObstacle<S>( mData, flags);
	else                 KnClampPositions<S>  ( mData, flags);
}



template <class S>  struct KnProjectParticles : public KernelBase { KnProjectParticles(ParticleSystem<S>& part, Grid<Vec3>& gradient) :  KernelBase(part.size()) ,part(part),gradient(gradient)   { run(); }  inline void op(int idx, ParticleSystem<S>& part, Grid<Vec3>& gradient )  {
	static RandomStream rand (3123984);
	const double jlen = 0.1;
	
	if (part.isActive(idx)) {
		// project along levelset gradient
		Vec3 p = part[idx].pos;
		if (gradient.isInBounds(p)) {
			Vec3 n = gradient.getInterpolated(p);
			Real dist = normalize(n);
			Vec3 dx = n * (-dist + jlen * (1 + rand.getReal()));
			p += dx;            
		}
		// clamp to outer boundaries (+jitter)
		const double jlen = 0.1;
		Vec3 jitter = jlen * rand.getVec3();
		part[idx].pos = clamp(p, Vec3(1,1,1)+jitter, toVec3(gradient.getSize()-1)-jitter);
	}
}   inline ParticleSystem<S>& getArg0() { return part; } typedef ParticleSystem<S> type0;inline Grid<Vec3>& getArg1() { return gradient; } typedef Grid<Vec3> type1; void run() {  const int _sz = size; for (int i=0; i < _sz; i++) op(i, part,gradient);  } ParticleSystem<S>& part; Grid<Vec3>& gradient;   };

template<class S>
void ParticleSystem<S>::projectOutside(Grid<Vec3>& gradient) {
	KnProjectParticles<S>(*this, gradient);
}

template<class S>
void ParticleSystem<S>::resizeAll(int size) {
	// resize all buffers to target size in 1 go
	mData.resize(size);
	for(int i=0; i<(int)mPartData.size(); ++i)
		mPartData[i]->resize(size);
}

template<class S>
void ParticleSystem<S>::compress() {
	int nextRead = mData.size();
	for (int i=0; i<(int)mData.size(); i++) {
		while ((mData[i].flag & PDELETE) != 0) {
			nextRead--;
			mData[i] = mData[nextRead];
			// ugly, but prevent virtual function calls here:
			for(int pd=0; pd<(int)mPdataReal.size(); ++pd) mPdataReal[pd]->copyValue(nextRead, i);
			for(int pd=0; pd<(int)mPdataVec3.size(); ++pd) mPdataVec3[pd]->copyValue(nextRead, i);
			for(int pd=0; pd<(int)mPdataInt .size(); ++pd) mPdataInt [pd]->copyValue(nextRead, i);
			mData[nextRead].flag = PINVALID;
		}
	}
	if(nextRead<(int)mData.size()) debMsg("Deleted "<<((int)mData.size() - nextRead)<<" particles", 1); // debug info

	resizeAll(nextRead);
	mDeletes = 0;
	mDeleteChunk = mData.size() / DELETE_PART;
}

//! insert buffered positions as new particles, update additional particle data
template<class S>
void ParticleSystem<S>::insertBufferedParticles() {
	if(mNewBuffer.size()==0) return;
	int newCnt = mData.size();
	resizeAll(newCnt + mNewBuffer.size());

	// clear new flag everywhere
	for(int i=0; i<(int)mData.size(); ++i) mData[i].flag &= ~PNEW;

	for(int i=0; i<(int)mNewBuffer.size(); ++i) {
		// note, other fields are not initialized here...
		mData[newCnt].pos  = mNewBuffer[i];
		mData[newCnt].flag = PNEW;
		// now init pdata fields from associated grids...
		for(int pd=0; pd<(int)mPdataReal.size(); ++pd) 
			mPdataReal[pd]->initNewValue(newCnt, mNewBuffer[i] );
		for(int pd=0; pd<(int)mPdataVec3.size(); ++pd) 
			mPdataVec3[pd]->initNewValue(newCnt, mNewBuffer[i] );
		for(int pd=0; pd<(int)mPdataInt.size(); ++pd) 
			mPdataInt[pd]->initNewValue(newCnt, mNewBuffer[i] );
		newCnt++;
	}
	if(mNewBuffer.size()>0) debMsg("Added & initialized "<<(int)mNewBuffer.size()<<" particles", 1); // debug info
	mNewBuffer.clear();
}


template<class DATA, class CON>
void ConnectedParticleSystem<DATA,CON>::compress() {
	const int sz = ParticleSystem<DATA>::size();
	int *renumber_back = new int[sz];
	int *renumber = new int[sz];
	for (int i=0; i<sz; i++)
		renumber[i] = renumber_back[i] = -1;
		
	// reorder elements
	std::vector<DATA>& data = ParticleSystem<DATA>::mData;
	int nextRead = sz;
	for (int i=0; i<nextRead; i++) {
		if ((data[i].flag & ParticleBase::PDELETE) != 0) {
			nextRead--;
			data[i] = data[nextRead];
			data[nextRead].flag = 0;           
			renumber_back[i] = nextRead;
		} else 
			renumber_back[i] = i;
	}
	
	// acceleration structure
	for (int i=0; i<nextRead; i++)
		renumber[renumber_back[i]] = i;
	
	// rename indices in filaments
	for (int i=0; i<(int)mSegments.size(); i++)
		mSegments[i].renumber(renumber);
		
	ParticleSystem<DATA>::mData.resize(nextRead);
	ParticleSystem<DATA>::mDeletes = 0;
	ParticleSystem<DATA>::mDeleteChunk = ParticleSystem<DATA>::size() / DELETE_PART;
	
	delete[] renumber;
	delete[] renumber_back;
}

template<class S>
ParticleBase* ParticleSystem<S>::clone() {
	ParticleSystem<S>* nm = new ParticleSystem<S>(getParent());
	if(this->mAllowCompress) compress();
	
	nm->mData = mData;
	nm->setName(getName());
	this->cloneParticleData(nm);
	return nm;
}

template<class DATA,class CON>
ParticleBase* ConnectedParticleSystem<DATA,CON>::clone() {
	ConnectedParticleSystem<DATA,CON>* nm = new ConnectedParticleSystem<DATA,CON>(this->getParent());
	if(this->mAllowCompress) compress();
	
	nm->mData = this->mData;
	nm->mSegments = mSegments;
	nm->setName(this->getName());
	this->cloneParticleData(nm);
	return nm;
}

template<class S>  
std::string ParticleSystem<S>::infoString() const { 
	std::stringstream s;
	s << "ParticleSys '" << getName() << "' [" << size() << " parts";
	if(this->getNumPdata()>0) s<< " "<< this->getNumPdata()<<" pd";
	s << "]";
	return s.str();
}
	
template<class S>  
inline void ParticleSystem<S>::checkPartIndex(int idx) const {
	int mySize = this->size();
	if (idx<0 || idx > mySize ) {
		errMsg( "ParticleBase " << " size " << mySize << " : index " << idx << " out of bound " );
	}
}
	
inline void ParticleDataBase::checkPartIndex(int idx) const {
	int mySize = this->size();
	if (idx<0 || idx > mySize ) {
		errMsg( "ParticleData " << " size " << mySize << " : index " << idx << " out of bound " );
	}
	if ( mpParticleSys && mpParticleSys->getSizeSlow()!=mySize ) {
		errMsg( "ParticleData " << " size " << mySize << " does not match parent! (" << mpParticleSys->getSizeSlow() << ") " );
	}
}

// set contents to zero, as for a grid
template<class T>
void ParticleDataImpl<T>::clear() {
	for(int i=0; i<(int)mData.size(); ++i) mData[i] = 0.;
}


} // namespace

#endif



