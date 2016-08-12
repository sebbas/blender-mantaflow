




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep link).




#include "particle.h"
namespace Manta {
#ifdef _C_BasicParticleSystem
 static const Pb::Register _R_3 ("BasicParticleSystem","BasicParticleSystem","ParticleSystem<BasicParticleData>"); template<> const char* Namify<BasicParticleSystem >::S = "BasicParticleSystem"; 
 static const Pb::Register _R_4 ("BasicParticleSystem","BasicParticleSystem",BasicParticleSystem::_W_10); 
 static const Pb::Register _R_5 ("BasicParticleSystem","save",BasicParticleSystem::_W_11); 
 static const Pb::Register _R_6 ("BasicParticleSystem","load",BasicParticleSystem::_W_12); 
 static const Pb::Register _R_7 ("BasicParticleSystem","addParticle",BasicParticleSystem::_W_13); 
 static const Pb::Register _R_8 ("BasicParticleSystem","printParts",BasicParticleSystem::_W_14); 
#endif
#ifdef _C_ParticleBase
 static const Pb::Register _R_9 ("ParticleBase","ParticleBase","PbClass"); template<> const char* Namify<ParticleBase >::S = "ParticleBase"; 
 static const Pb::Register _R_10 ("ParticleBase","ParticleBase",ParticleBase::_W_0); 
 static const Pb::Register _R_11 ("ParticleBase","create",ParticleBase::_W_1); 
#endif
#ifdef _C_ParticleDataBase
 static const Pb::Register _R_12 ("ParticleDataBase","ParticleDataBase","PbClass"); template<> const char* Namify<ParticleDataBase >::S = "ParticleDataBase"; 
 static const Pb::Register _R_13 ("ParticleDataBase","ParticleDataBase",ParticleDataBase::_W_17); 
#endif
#ifdef _C_ParticleDataImpl
 static const Pb::Register _R_14 ("ParticleDataImpl<int>","ParticleDataImpl<int>","ParticleDataBase"); template<> const char* Namify<ParticleDataImpl<int> >::S = "ParticleDataImpl<int>"; 
 static const Pb::Register _R_15 ("ParticleDataImpl<int>","ParticleDataImpl",ParticleDataImpl<int>::_W_18); 
 static const Pb::Register _R_16 ("ParticleDataImpl<int>","clear",ParticleDataImpl<int>::_W_19); 
 static const Pb::Register _R_17 ("ParticleDataImpl<int>","setSource",ParticleDataImpl<int>::_W_20); 
 static const Pb::Register _R_18 ("ParticleDataImpl<int>","setConst",ParticleDataImpl<int>::_W_21); 
 static const Pb::Register _R_19 ("ParticleDataImpl<int>","copyFrom",ParticleDataImpl<int>::_W_22); 
 static const Pb::Register _R_20 ("ParticleDataImpl<int>","add",ParticleDataImpl<int>::_W_23); 
 static const Pb::Register _R_21 ("ParticleDataImpl<int>","sub",ParticleDataImpl<int>::_W_24); 
 static const Pb::Register _R_22 ("ParticleDataImpl<int>","addConst",ParticleDataImpl<int>::_W_25); 
 static const Pb::Register _R_23 ("ParticleDataImpl<int>","addScaled",ParticleDataImpl<int>::_W_26); 
 static const Pb::Register _R_24 ("ParticleDataImpl<int>","mult",ParticleDataImpl<int>::_W_27); 
 static const Pb::Register _R_25 ("ParticleDataImpl<int>","multConst",ParticleDataImpl<int>::_W_28); 
 static const Pb::Register _R_26 ("ParticleDataImpl<int>","clamp",ParticleDataImpl<int>::_W_29); 
 static const Pb::Register _R_27 ("ParticleDataImpl<int>","getMaxAbsValue",ParticleDataImpl<int>::_W_30); 
 static const Pb::Register _R_28 ("ParticleDataImpl<int>","getMaxValue",ParticleDataImpl<int>::_W_31); 
 static const Pb::Register _R_29 ("ParticleDataImpl<int>","getMinValue",ParticleDataImpl<int>::_W_32); 
 static const Pb::Register _R_30 ("ParticleDataImpl<int>","printPdata",ParticleDataImpl<int>::_W_33); 
 static const Pb::Register _R_31 ("ParticleDataImpl<int>","save",ParticleDataImpl<int>::_W_34); 
 static const Pb::Register _R_32 ("ParticleDataImpl<int>","load",ParticleDataImpl<int>::_W_35); 
 static const Pb::Register _R_33 ("ParticleDataImpl<Real>","ParticleDataImpl<Real>","ParticleDataBase"); template<> const char* Namify<ParticleDataImpl<Real> >::S = "ParticleDataImpl<Real>"; 
 static const Pb::Register _R_34 ("ParticleDataImpl<Real>","ParticleDataImpl",ParticleDataImpl<Real>::_W_18); 
 static const Pb::Register _R_35 ("ParticleDataImpl<Real>","clear",ParticleDataImpl<Real>::_W_19); 
 static const Pb::Register _R_36 ("ParticleDataImpl<Real>","setSource",ParticleDataImpl<Real>::_W_20); 
 static const Pb::Register _R_37 ("ParticleDataImpl<Real>","setConst",ParticleDataImpl<Real>::_W_21); 
 static const Pb::Register _R_38 ("ParticleDataImpl<Real>","copyFrom",ParticleDataImpl<Real>::_W_22); 
 static const Pb::Register _R_39 ("ParticleDataImpl<Real>","add",ParticleDataImpl<Real>::_W_23); 
 static const Pb::Register _R_40 ("ParticleDataImpl<Real>","sub",ParticleDataImpl<Real>::_W_24); 
 static const Pb::Register _R_41 ("ParticleDataImpl<Real>","addConst",ParticleDataImpl<Real>::_W_25); 
 static const Pb::Register _R_42 ("ParticleDataImpl<Real>","addScaled",ParticleDataImpl<Real>::_W_26); 
 static const Pb::Register _R_43 ("ParticleDataImpl<Real>","mult",ParticleDataImpl<Real>::_W_27); 
 static const Pb::Register _R_44 ("ParticleDataImpl<Real>","multConst",ParticleDataImpl<Real>::_W_28); 
 static const Pb::Register _R_45 ("ParticleDataImpl<Real>","clamp",ParticleDataImpl<Real>::_W_29); 
 static const Pb::Register _R_46 ("ParticleDataImpl<Real>","getMaxAbsValue",ParticleDataImpl<Real>::_W_30); 
 static const Pb::Register _R_47 ("ParticleDataImpl<Real>","getMaxValue",ParticleDataImpl<Real>::_W_31); 
 static const Pb::Register _R_48 ("ParticleDataImpl<Real>","getMinValue",ParticleDataImpl<Real>::_W_32); 
 static const Pb::Register _R_49 ("ParticleDataImpl<Real>","printPdata",ParticleDataImpl<Real>::_W_33); 
 static const Pb::Register _R_50 ("ParticleDataImpl<Real>","save",ParticleDataImpl<Real>::_W_34); 
 static const Pb::Register _R_51 ("ParticleDataImpl<Real>","load",ParticleDataImpl<Real>::_W_35); 
 static const Pb::Register _R_52 ("ParticleDataImpl<Vec3>","ParticleDataImpl<Vec3>","ParticleDataBase"); template<> const char* Namify<ParticleDataImpl<Vec3> >::S = "ParticleDataImpl<Vec3>"; 
 static const Pb::Register _R_53 ("ParticleDataImpl<Vec3>","ParticleDataImpl",ParticleDataImpl<Vec3>::_W_18); 
 static const Pb::Register _R_54 ("ParticleDataImpl<Vec3>","clear",ParticleDataImpl<Vec3>::_W_19); 
 static const Pb::Register _R_55 ("ParticleDataImpl<Vec3>","setSource",ParticleDataImpl<Vec3>::_W_20); 
 static const Pb::Register _R_56 ("ParticleDataImpl<Vec3>","setConst",ParticleDataImpl<Vec3>::_W_21); 
 static const Pb::Register _R_57 ("ParticleDataImpl<Vec3>","copyFrom",ParticleDataImpl<Vec3>::_W_22); 
 static const Pb::Register _R_58 ("ParticleDataImpl<Vec3>","add",ParticleDataImpl<Vec3>::_W_23); 
 static const Pb::Register _R_59 ("ParticleDataImpl<Vec3>","sub",ParticleDataImpl<Vec3>::_W_24); 
 static const Pb::Register _R_60 ("ParticleDataImpl<Vec3>","addConst",ParticleDataImpl<Vec3>::_W_25); 
 static const Pb::Register _R_61 ("ParticleDataImpl<Vec3>","addScaled",ParticleDataImpl<Vec3>::_W_26); 
 static const Pb::Register _R_62 ("ParticleDataImpl<Vec3>","mult",ParticleDataImpl<Vec3>::_W_27); 
 static const Pb::Register _R_63 ("ParticleDataImpl<Vec3>","multConst",ParticleDataImpl<Vec3>::_W_28); 
 static const Pb::Register _R_64 ("ParticleDataImpl<Vec3>","clamp",ParticleDataImpl<Vec3>::_W_29); 
 static const Pb::Register _R_65 ("ParticleDataImpl<Vec3>","getMaxAbsValue",ParticleDataImpl<Vec3>::_W_30); 
 static const Pb::Register _R_66 ("ParticleDataImpl<Vec3>","getMaxValue",ParticleDataImpl<Vec3>::_W_31); 
 static const Pb::Register _R_67 ("ParticleDataImpl<Vec3>","getMinValue",ParticleDataImpl<Vec3>::_W_32); 
 static const Pb::Register _R_68 ("ParticleDataImpl<Vec3>","printPdata",ParticleDataImpl<Vec3>::_W_33); 
 static const Pb::Register _R_69 ("ParticleDataImpl<Vec3>","save",ParticleDataImpl<Vec3>::_W_34); 
 static const Pb::Register _R_70 ("ParticleDataImpl<Vec3>","load",ParticleDataImpl<Vec3>::_W_35); 
#endif
#ifdef _C_ParticleIndexSystem
 static const Pb::Register _R_71 ("ParticleIndexSystem","ParticleIndexSystem","ParticleSystem<ParticleIndexData>"); template<> const char* Namify<ParticleIndexSystem >::S = "ParticleIndexSystem"; 
 static const Pb::Register _R_72 ("ParticleIndexSystem","ParticleIndexSystem",ParticleIndexSystem::_W_15); 
#endif
#ifdef _C_ParticleSystem
 static const Pb::Register _R_73 ("ParticleSystem<BasicParticleData>","ParticleSystem<BasicParticleData>","ParticleBase"); template<> const char* Namify<ParticleSystem<BasicParticleData> >::S = "ParticleSystem<BasicParticleData>"; 
 static const Pb::Register _R_74 ("ParticleSystem<BasicParticleData>","ParticleSystem",ParticleSystem<BasicParticleData>::_W_2); 
 static const Pb::Register _R_75 ("ParticleSystem<BasicParticleData>","setPos",ParticleSystem<BasicParticleData>::_W_3); 
 static const Pb::Register _R_76 ("ParticleSystem<BasicParticleData>","getPos",ParticleSystem<BasicParticleData>::_W_4); 
 static const Pb::Register _R_77 ("ParticleSystem<BasicParticleData>","getPosPdata",ParticleSystem<BasicParticleData>::_W_5); 
 static const Pb::Register _R_78 ("ParticleSystem<BasicParticleData>","setPosPdata",ParticleSystem<BasicParticleData>::_W_6); 
 static const Pb::Register _R_79 ("ParticleSystem<BasicParticleData>","clear",ParticleSystem<BasicParticleData>::_W_7); 
 static const Pb::Register _R_80 ("ParticleSystem<BasicParticleData>","advectInGrid",ParticleSystem<BasicParticleData>::_W_8); 
 static const Pb::Register _R_81 ("ParticleSystem<BasicParticleData>","projectOutside",ParticleSystem<BasicParticleData>::_W_9); 
 static const Pb::Register _R_82 ("ParticleSystem<ParticleIndexData>","ParticleSystem<ParticleIndexData>","ParticleBase"); template<> const char* Namify<ParticleSystem<ParticleIndexData> >::S = "ParticleSystem<ParticleIndexData>"; 
 static const Pb::Register _R_83 ("ParticleSystem<ParticleIndexData>","ParticleSystem",ParticleSystem<ParticleIndexData>::_W_2); 
 static const Pb::Register _R_84 ("ParticleSystem<ParticleIndexData>","setPos",ParticleSystem<ParticleIndexData>::_W_3); 
 static const Pb::Register _R_85 ("ParticleSystem<ParticleIndexData>","getPos",ParticleSystem<ParticleIndexData>::_W_4); 
 static const Pb::Register _R_86 ("ParticleSystem<ParticleIndexData>","getPosPdata",ParticleSystem<ParticleIndexData>::_W_5); 
 static const Pb::Register _R_87 ("ParticleSystem<ParticleIndexData>","setPosPdata",ParticleSystem<ParticleIndexData>::_W_6); 
 static const Pb::Register _R_88 ("ParticleSystem<ParticleIndexData>","clear",ParticleSystem<ParticleIndexData>::_W_7); 
 static const Pb::Register _R_89 ("ParticleSystem<ParticleIndexData>","advectInGrid",ParticleSystem<ParticleIndexData>::_W_8); 
 static const Pb::Register _R_90 ("ParticleSystem<ParticleIndexData>","projectOutside",ParticleSystem<ParticleIndexData>::_W_9); 
#endif
static const Pb::Register _R_0 ("ParticleDataImpl<int>","PdataInt","");
static const Pb::Register _R_1 ("ParticleDataImpl<Real>","PdataReal","");
static const Pb::Register _R_2 ("ParticleDataImpl<Vec3>","PdataVec3","");
}