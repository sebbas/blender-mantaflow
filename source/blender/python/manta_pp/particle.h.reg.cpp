




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep link).




#include "particle.h"
namespace Manta {
#ifdef _C_BasicParticleSystem
 static const Pb::Register _R_3 ("BasicParticleSystem","BasicParticleSystem","ParticleSystem<BasicParticleData>"); template<> const char* Namify<BasicParticleSystem >::S = "BasicParticleSystem"; 
 static const Pb::Register _R_4 ("BasicParticleSystem","BasicParticleSystem",BasicParticleSystem::_W_11); 
 static const Pb::Register _R_5 ("BasicParticleSystem","save",BasicParticleSystem::_W_12); 
 static const Pb::Register _R_6 ("BasicParticleSystem","load",BasicParticleSystem::_W_13); 
 static const Pb::Register _R_7 ("BasicParticleSystem","addParticle",BasicParticleSystem::_W_14); 
#endif
#ifdef _C_ParticleBase
 static const Pb::Register _R_8 ("ParticleBase","ParticleBase","PbClass"); template<> const char* Namify<ParticleBase >::S = "ParticleBase"; 
 static const Pb::Register _R_9 ("ParticleBase","ParticleBase",ParticleBase::_W_0); 
 static const Pb::Register _R_10 ("ParticleBase","create",ParticleBase::_W_1); 
#endif
#ifdef _C_ParticleDataBase
 static const Pb::Register _R_11 ("ParticleDataBase","ParticleDataBase","PbClass"); template<> const char* Namify<ParticleDataBase >::S = "ParticleDataBase"; 
 static const Pb::Register _R_12 ("ParticleDataBase","ParticleDataBase",ParticleDataBase::_W_17); 
#endif
#ifdef _C_ParticleDataImpl
 static const Pb::Register _R_13 ("ParticleDataImpl<int>","ParticleDataImpl<int>","ParticleDataBase"); template<> const char* Namify<ParticleDataImpl<int> >::S = "ParticleDataImpl<int>"; 
 static const Pb::Register _R_14 ("ParticleDataImpl<int>","ParticleDataImpl",ParticleDataImpl<int>::_W_18); 
 static const Pb::Register _R_15 ("ParticleDataImpl<int>","clear",ParticleDataImpl<int>::_W_19); 
 static const Pb::Register _R_16 ("ParticleDataImpl<int>","setSource",ParticleDataImpl<int>::_W_20); 
 static const Pb::Register _R_17 ("ParticleDataImpl<int>","save",ParticleDataImpl<int>::_W_21); 
 static const Pb::Register _R_18 ("ParticleDataImpl<int>","load",ParticleDataImpl<int>::_W_22); 
 static const Pb::Register _R_19 ("ParticleDataImpl<Real>","ParticleDataImpl<Real>","ParticleDataBase"); template<> const char* Namify<ParticleDataImpl<Real> >::S = "ParticleDataImpl<Real>"; 
 static const Pb::Register _R_20 ("ParticleDataImpl<Real>","ParticleDataImpl",ParticleDataImpl<Real>::_W_18); 
 static const Pb::Register _R_21 ("ParticleDataImpl<Real>","clear",ParticleDataImpl<Real>::_W_19); 
 static const Pb::Register _R_22 ("ParticleDataImpl<Real>","setSource",ParticleDataImpl<Real>::_W_20); 
 static const Pb::Register _R_23 ("ParticleDataImpl<Real>","save",ParticleDataImpl<Real>::_W_21); 
 static const Pb::Register _R_24 ("ParticleDataImpl<Real>","load",ParticleDataImpl<Real>::_W_22); 
 static const Pb::Register _R_25 ("ParticleDataImpl<Vec3>","ParticleDataImpl<Vec3>","ParticleDataBase"); template<> const char* Namify<ParticleDataImpl<Vec3> >::S = "ParticleDataImpl<Vec3>"; 
 static const Pb::Register _R_26 ("ParticleDataImpl<Vec3>","ParticleDataImpl",ParticleDataImpl<Vec3>::_W_18); 
 static const Pb::Register _R_27 ("ParticleDataImpl<Vec3>","clear",ParticleDataImpl<Vec3>::_W_19); 
 static const Pb::Register _R_28 ("ParticleDataImpl<Vec3>","setSource",ParticleDataImpl<Vec3>::_W_20); 
 static const Pb::Register _R_29 ("ParticleDataImpl<Vec3>","save",ParticleDataImpl<Vec3>::_W_21); 
 static const Pb::Register _R_30 ("ParticleDataImpl<Vec3>","load",ParticleDataImpl<Vec3>::_W_22); 
#endif
#ifdef _C_ParticleIndexSystem
 static const Pb::Register _R_31 ("ParticleIndexSystem","ParticleIndexSystem","ParticleSystem<ParticleIndexData>"); template<> const char* Namify<ParticleIndexSystem >::S = "ParticleIndexSystem"; 
 static const Pb::Register _R_32 ("ParticleIndexSystem","ParticleIndexSystem",ParticleIndexSystem::_W_15); 
#endif
#ifdef _C_ParticleSystem
 static const Pb::Register _R_33 ("ParticleSystem<BasicParticleData>","ParticleSystem<BasicParticleData>","ParticleBase"); template<> const char* Namify<ParticleSystem<BasicParticleData> >::S = "ParticleSystem<BasicParticleData>"; 
 static const Pb::Register _R_34 ("ParticleSystem<BasicParticleData>","ParticleSystem",ParticleSystem<BasicParticleData>::_W_2); 
 static const Pb::Register _R_35 ("ParticleSystem<BasicParticleData>","size",ParticleSystem<BasicParticleData>::_W_3); 
 static const Pb::Register _R_36 ("ParticleSystem<BasicParticleData>","setPos",ParticleSystem<BasicParticleData>::_W_4); 
 static const Pb::Register _R_37 ("ParticleSystem<BasicParticleData>","getPos",ParticleSystem<BasicParticleData>::_W_5); 
 static const Pb::Register _R_38 ("ParticleSystem<BasicParticleData>","getPosPdata",ParticleSystem<BasicParticleData>::_W_6); 
 static const Pb::Register _R_39 ("ParticleSystem<BasicParticleData>","setPosPdata",ParticleSystem<BasicParticleData>::_W_7); 
 static const Pb::Register _R_40 ("ParticleSystem<BasicParticleData>","clear",ParticleSystem<BasicParticleData>::_W_8); 
 static const Pb::Register _R_41 ("ParticleSystem<BasicParticleData>","advectInGrid",ParticleSystem<BasicParticleData>::_W_9); 
 static const Pb::Register _R_42 ("ParticleSystem<BasicParticleData>","projectOutside",ParticleSystem<BasicParticleData>::_W_10); 
 static const Pb::Register _R_43 ("ParticleSystem<ParticleIndexData>","ParticleSystem<ParticleIndexData>","ParticleBase"); template<> const char* Namify<ParticleSystem<ParticleIndexData> >::S = "ParticleSystem<ParticleIndexData>"; 
 static const Pb::Register _R_44 ("ParticleSystem<ParticleIndexData>","ParticleSystem",ParticleSystem<ParticleIndexData>::_W_2); 
 static const Pb::Register _R_45 ("ParticleSystem<ParticleIndexData>","size",ParticleSystem<ParticleIndexData>::_W_3); 
 static const Pb::Register _R_46 ("ParticleSystem<ParticleIndexData>","setPos",ParticleSystem<ParticleIndexData>::_W_4); 
 static const Pb::Register _R_47 ("ParticleSystem<ParticleIndexData>","getPos",ParticleSystem<ParticleIndexData>::_W_5); 
 static const Pb::Register _R_48 ("ParticleSystem<ParticleIndexData>","getPosPdata",ParticleSystem<ParticleIndexData>::_W_6); 
 static const Pb::Register _R_49 ("ParticleSystem<ParticleIndexData>","setPosPdata",ParticleSystem<ParticleIndexData>::_W_7); 
 static const Pb::Register _R_50 ("ParticleSystem<ParticleIndexData>","clear",ParticleSystem<ParticleIndexData>::_W_8); 
 static const Pb::Register _R_51 ("ParticleSystem<ParticleIndexData>","advectInGrid",ParticleSystem<ParticleIndexData>::_W_9); 
 static const Pb::Register _R_52 ("ParticleSystem<ParticleIndexData>","projectOutside",ParticleSystem<ParticleIndexData>::_W_10); 
#endif
static const Pb::Register _R_0 ("ParticleDataImpl<int>","PdataInt","");
static const Pb::Register _R_1 ("ParticleDataImpl<Real>","PdataReal","");
static const Pb::Register _R_2 ("ParticleDataImpl<Vec3>","PdataVec3","");
}