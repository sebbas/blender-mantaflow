




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep link).




#include "grid.h"
namespace Manta {
#ifdef _C_FlagGrid
 static const Pb::Register _R_3 ("FlagGrid","FlagGrid","Grid<int>"); template<> const char* Namify<FlagGrid >::S = "FlagGrid"; 
 static const Pb::Register _R_4 ("FlagGrid","FlagGrid",FlagGrid::_W_25); 
 static const Pb::Register _R_5 ("FlagGrid","initDomain",FlagGrid::_W_26); 
 static const Pb::Register _R_6 ("FlagGrid","updateFromLevelset",FlagGrid::_W_27); 
 static const Pb::Register _R_7 ("FlagGrid","fillGrid",FlagGrid::_W_28); 
#endif
#ifdef _C_Grid
 static const Pb::Register _R_8 ("Grid<int>","Grid<int>","GridBase"); template<> const char* Namify<Grid<int> >::S = "Grid<int>"; 
 static const Pb::Register _R_9 ("Grid<int>","Grid",Grid<int>::_W_1); 
 static const Pb::Register _R_10 ("Grid<int>","save",Grid<int>::_W_2); 
 static const Pb::Register _R_11 ("Grid<int>","load",Grid<int>::_W_3); 
 static const Pb::Register _R_12 ("Grid<int>","clear",Grid<int>::_W_4); 
 static const Pb::Register _R_13 ("Grid<int>","copyFrom",Grid<int>::_W_5); 
 static const Pb::Register _R_14 ("Grid<int>","add",Grid<int>::_W_6); 
 static const Pb::Register _R_15 ("Grid<int>","sub",Grid<int>::_W_7); 
 static const Pb::Register _R_16 ("Grid<int>","setConst",Grid<int>::_W_8); 
 static const Pb::Register _R_17 ("Grid<int>","addConst",Grid<int>::_W_9); 
 static const Pb::Register _R_18 ("Grid<int>","addScaled",Grid<int>::_W_10); 
 static const Pb::Register _R_19 ("Grid<int>","mult",Grid<int>::_W_11); 
 static const Pb::Register _R_20 ("Grid<int>","multConst",Grid<int>::_W_12); 
 static const Pb::Register _R_21 ("Grid<int>","clamp",Grid<int>::_W_13); 
 static const Pb::Register _R_22 ("Grid<int>","getMaxAbs",Grid<int>::_W_14); 
 static const Pb::Register _R_23 ("Grid<int>","getMax",Grid<int>::_W_15); 
 static const Pb::Register _R_24 ("Grid<int>","getMin",Grid<int>::_W_16); 
 static const Pb::Register _R_25 ("Grid<int>","setBound",Grid<int>::_W_17); 
 static const Pb::Register _R_26 ("Grid<int>","setBoundNeumann",Grid<int>::_W_18); 
 static const Pb::Register _R_27 ("Grid<int>","getMaxAbsValue",Grid<int>::_W_19); 
 static const Pb::Register _R_28 ("Grid<int>","getMaxValue",Grid<int>::_W_20); 
 static const Pb::Register _R_29 ("Grid<int>","getMinValue",Grid<int>::_W_21); 
 static const Pb::Register _R_30 ("Grid<int>","getDataPointer",Grid<int>::_W_22); 
 static const Pb::Register _R_31 ("Grid<int>","printGrid",Grid<int>::_W_23); 
 static const Pb::Register _R_32 ("Grid<Real>","Grid<Real>","GridBase"); template<> const char* Namify<Grid<Real> >::S = "Grid<Real>"; 
 static const Pb::Register _R_33 ("Grid<Real>","Grid",Grid<Real>::_W_1); 
 static const Pb::Register _R_34 ("Grid<Real>","save",Grid<Real>::_W_2); 
 static const Pb::Register _R_35 ("Grid<Real>","load",Grid<Real>::_W_3); 
 static const Pb::Register _R_36 ("Grid<Real>","clear",Grid<Real>::_W_4); 
 static const Pb::Register _R_37 ("Grid<Real>","copyFrom",Grid<Real>::_W_5); 
 static const Pb::Register _R_38 ("Grid<Real>","add",Grid<Real>::_W_6); 
 static const Pb::Register _R_39 ("Grid<Real>","sub",Grid<Real>::_W_7); 
 static const Pb::Register _R_40 ("Grid<Real>","setConst",Grid<Real>::_W_8); 
 static const Pb::Register _R_41 ("Grid<Real>","addConst",Grid<Real>::_W_9); 
 static const Pb::Register _R_42 ("Grid<Real>","addScaled",Grid<Real>::_W_10); 
 static const Pb::Register _R_43 ("Grid<Real>","mult",Grid<Real>::_W_11); 
 static const Pb::Register _R_44 ("Grid<Real>","multConst",Grid<Real>::_W_12); 
 static const Pb::Register _R_45 ("Grid<Real>","clamp",Grid<Real>::_W_13); 
 static const Pb::Register _R_46 ("Grid<Real>","getMaxAbs",Grid<Real>::_W_14); 
 static const Pb::Register _R_47 ("Grid<Real>","getMax",Grid<Real>::_W_15); 
 static const Pb::Register _R_48 ("Grid<Real>","getMin",Grid<Real>::_W_16); 
 static const Pb::Register _R_49 ("Grid<Real>","setBound",Grid<Real>::_W_17); 
 static const Pb::Register _R_50 ("Grid<Real>","setBoundNeumann",Grid<Real>::_W_18); 
 static const Pb::Register _R_51 ("Grid<Real>","getMaxAbsValue",Grid<Real>::_W_19); 
 static const Pb::Register _R_52 ("Grid<Real>","getMaxValue",Grid<Real>::_W_20); 
 static const Pb::Register _R_53 ("Grid<Real>","getMinValue",Grid<Real>::_W_21); 
 static const Pb::Register _R_54 ("Grid<Real>","getDataPointer",Grid<Real>::_W_22); 
 static const Pb::Register _R_55 ("Grid<Real>","printGrid",Grid<Real>::_W_23); 
 static const Pb::Register _R_56 ("Grid<Vec3>","Grid<Vec3>","GridBase"); template<> const char* Namify<Grid<Vec3> >::S = "Grid<Vec3>"; 
 static const Pb::Register _R_57 ("Grid<Vec3>","Grid",Grid<Vec3>::_W_1); 
 static const Pb::Register _R_58 ("Grid<Vec3>","save",Grid<Vec3>::_W_2); 
 static const Pb::Register _R_59 ("Grid<Vec3>","load",Grid<Vec3>::_W_3); 
 static const Pb::Register _R_60 ("Grid<Vec3>","clear",Grid<Vec3>::_W_4); 
 static const Pb::Register _R_61 ("Grid<Vec3>","copyFrom",Grid<Vec3>::_W_5); 
 static const Pb::Register _R_62 ("Grid<Vec3>","add",Grid<Vec3>::_W_6); 
 static const Pb::Register _R_63 ("Grid<Vec3>","sub",Grid<Vec3>::_W_7); 
 static const Pb::Register _R_64 ("Grid<Vec3>","setConst",Grid<Vec3>::_W_8); 
 static const Pb::Register _R_65 ("Grid<Vec3>","addConst",Grid<Vec3>::_W_9); 
 static const Pb::Register _R_66 ("Grid<Vec3>","addScaled",Grid<Vec3>::_W_10); 
 static const Pb::Register _R_67 ("Grid<Vec3>","mult",Grid<Vec3>::_W_11); 
 static const Pb::Register _R_68 ("Grid<Vec3>","multConst",Grid<Vec3>::_W_12); 
 static const Pb::Register _R_69 ("Grid<Vec3>","clamp",Grid<Vec3>::_W_13); 
 static const Pb::Register _R_70 ("Grid<Vec3>","getMaxAbs",Grid<Vec3>::_W_14); 
 static const Pb::Register _R_71 ("Grid<Vec3>","getMax",Grid<Vec3>::_W_15); 
 static const Pb::Register _R_72 ("Grid<Vec3>","getMin",Grid<Vec3>::_W_16); 
 static const Pb::Register _R_73 ("Grid<Vec3>","setBound",Grid<Vec3>::_W_17); 
 static const Pb::Register _R_74 ("Grid<Vec3>","setBoundNeumann",Grid<Vec3>::_W_18); 
 static const Pb::Register _R_75 ("Grid<Vec3>","getMaxAbsValue",Grid<Vec3>::_W_19); 
 static const Pb::Register _R_76 ("Grid<Vec3>","getMaxValue",Grid<Vec3>::_W_20); 
 static const Pb::Register _R_77 ("Grid<Vec3>","getMinValue",Grid<Vec3>::_W_21); 
 static const Pb::Register _R_78 ("Grid<Vec3>","getDataPointer",Grid<Vec3>::_W_22); 
 static const Pb::Register _R_79 ("Grid<Vec3>","printGrid",Grid<Vec3>::_W_23); 
#endif
#ifdef _C_GridBase
 static const Pb::Register _R_80 ("GridBase","GridBase","PbClass"); template<> const char* Namify<GridBase >::S = "GridBase"; 
 static const Pb::Register _R_81 ("GridBase","GridBase",GridBase::_W_0); 
#endif
#ifdef _C_MACGrid
 static const Pb::Register _R_82 ("MACGrid","MACGrid","Grid<Vec3>"); template<> const char* Namify<MACGrid >::S = "MACGrid"; 
 static const Pb::Register _R_83 ("MACGrid","MACGrid",MACGrid::_W_24); 
#endif
static const Pb::Register _R_0 ("Grid<int>","IntGrid","");
static const Pb::Register _R_1 ("Grid<Real>","RealGrid","");
static const Pb::Register _R_2 ("Grid<Vec3>","VecGrid","");
}