




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep link).




#include "grid.h"
namespace Manta {
#ifdef _C_FlagGrid
 static const Pb::Register _R_3 ("FlagGrid","FlagGrid","Grid<int>"); template<> const char* Namify<FlagGrid >::S = "FlagGrid"; 
 static const Pb::Register _R_4 ("FlagGrid","FlagGrid",FlagGrid::_W_20); 
 static const Pb::Register _R_5 ("FlagGrid","initDomain",FlagGrid::_W_21); 
 static const Pb::Register _R_6 ("FlagGrid","initBoundaries",FlagGrid::_W_22); 
 static const Pb::Register _R_7 ("FlagGrid","updateFromLevelset",FlagGrid::_W_23); 
 static const Pb::Register _R_8 ("FlagGrid","fillGrid",FlagGrid::_W_24); 
#endif
#ifdef _C_Grid
 static const Pb::Register _R_9 ("Grid<int>","Grid<int>","GridBase"); template<> const char* Namify<Grid<int> >::S = "Grid<int>"; 
 static const Pb::Register _R_10 ("Grid<int>","Grid",Grid<int>::_W_1); 
 static const Pb::Register _R_11 ("Grid<int>","save",Grid<int>::_W_2); 
 static const Pb::Register _R_12 ("Grid<int>","load",Grid<int>::_W_3); 
 static const Pb::Register _R_13 ("Grid<int>","clear",Grid<int>::_W_4); 
 static const Pb::Register _R_14 ("Grid<int>","add",Grid<int>::_W_5); 
 static const Pb::Register _R_15 ("Grid<int>","sub",Grid<int>::_W_6); 
 static const Pb::Register _R_16 ("Grid<int>","setAdd",Grid<int>::_W_7); 
 static const Pb::Register _R_17 ("Grid<int>","setSub",Grid<int>::_W_8); 
 static const Pb::Register _R_18 ("Grid<int>","addConstReal",Grid<int>::_W_9); 
 static const Pb::Register _R_19 ("Grid<int>","multiply",Grid<int>::_W_10); 
 static const Pb::Register _R_20 ("Grid<int>","multiplyConstReal",Grid<int>::_W_11); 
 static const Pb::Register _R_21 ("Grid<int>","addScaledReal",Grid<int>::_W_12); 
 static const Pb::Register _R_22 ("Grid<int>","copyFrom",Grid<int>::_W_13); 
 static const Pb::Register _R_23 ("Grid<int>","clamp",Grid<int>::_W_14); 
 static const Pb::Register _R_24 ("Grid<int>","getMaxAbsValue",Grid<int>::_W_15); 
 static const Pb::Register _R_25 ("Grid<int>","getMaxValue",Grid<int>::_W_16); 
 static const Pb::Register _R_26 ("Grid<int>","getMinValue",Grid<int>::_W_17); 
 static const Pb::Register _R_27 ("Grid<int>","printGrid",Grid<int>::_W_18); 
 static const Pb::Register _R_28 ("Grid<Real>","Grid<Real>","GridBase"); template<> const char* Namify<Grid<Real> >::S = "Grid<Real>"; 
 static const Pb::Register _R_29 ("Grid<Real>","Grid",Grid<Real>::_W_1); 
 static const Pb::Register _R_30 ("Grid<Real>","save",Grid<Real>::_W_2); 
 static const Pb::Register _R_31 ("Grid<Real>","load",Grid<Real>::_W_3); 
 static const Pb::Register _R_32 ("Grid<Real>","clear",Grid<Real>::_W_4); 
 static const Pb::Register _R_33 ("Grid<Real>","add",Grid<Real>::_W_5); 
 static const Pb::Register _R_34 ("Grid<Real>","sub",Grid<Real>::_W_6); 
 static const Pb::Register _R_35 ("Grid<Real>","setAdd",Grid<Real>::_W_7); 
 static const Pb::Register _R_36 ("Grid<Real>","setSub",Grid<Real>::_W_8); 
 static const Pb::Register _R_37 ("Grid<Real>","addConstReal",Grid<Real>::_W_9); 
 static const Pb::Register _R_38 ("Grid<Real>","multiply",Grid<Real>::_W_10); 
 static const Pb::Register _R_39 ("Grid<Real>","multiplyConstReal",Grid<Real>::_W_11); 
 static const Pb::Register _R_40 ("Grid<Real>","addScaledReal",Grid<Real>::_W_12); 
 static const Pb::Register _R_41 ("Grid<Real>","copyFrom",Grid<Real>::_W_13); 
 static const Pb::Register _R_42 ("Grid<Real>","clamp",Grid<Real>::_W_14); 
 static const Pb::Register _R_43 ("Grid<Real>","getMaxAbsValue",Grid<Real>::_W_15); 
 static const Pb::Register _R_44 ("Grid<Real>","getMaxValue",Grid<Real>::_W_16); 
 static const Pb::Register _R_45 ("Grid<Real>","getMinValue",Grid<Real>::_W_17); 
 static const Pb::Register _R_46 ("Grid<Real>","printGrid",Grid<Real>::_W_18); 
 static const Pb::Register _R_47 ("Grid<Vec3>","Grid<Vec3>","GridBase"); template<> const char* Namify<Grid<Vec3> >::S = "Grid<Vec3>"; 
 static const Pb::Register _R_48 ("Grid<Vec3>","Grid",Grid<Vec3>::_W_1); 
 static const Pb::Register _R_49 ("Grid<Vec3>","save",Grid<Vec3>::_W_2); 
 static const Pb::Register _R_50 ("Grid<Vec3>","load",Grid<Vec3>::_W_3); 
 static const Pb::Register _R_51 ("Grid<Vec3>","clear",Grid<Vec3>::_W_4); 
 static const Pb::Register _R_52 ("Grid<Vec3>","add",Grid<Vec3>::_W_5); 
 static const Pb::Register _R_53 ("Grid<Vec3>","sub",Grid<Vec3>::_W_6); 
 static const Pb::Register _R_54 ("Grid<Vec3>","setAdd",Grid<Vec3>::_W_7); 
 static const Pb::Register _R_55 ("Grid<Vec3>","setSub",Grid<Vec3>::_W_8); 
 static const Pb::Register _R_56 ("Grid<Vec3>","addConstReal",Grid<Vec3>::_W_9); 
 static const Pb::Register _R_57 ("Grid<Vec3>","multiply",Grid<Vec3>::_W_10); 
 static const Pb::Register _R_58 ("Grid<Vec3>","multiplyConstReal",Grid<Vec3>::_W_11); 
 static const Pb::Register _R_59 ("Grid<Vec3>","addScaledReal",Grid<Vec3>::_W_12); 
 static const Pb::Register _R_60 ("Grid<Vec3>","copyFrom",Grid<Vec3>::_W_13); 
 static const Pb::Register _R_61 ("Grid<Vec3>","clamp",Grid<Vec3>::_W_14); 
 static const Pb::Register _R_62 ("Grid<Vec3>","getMaxAbsValue",Grid<Vec3>::_W_15); 
 static const Pb::Register _R_63 ("Grid<Vec3>","getMaxValue",Grid<Vec3>::_W_16); 
 static const Pb::Register _R_64 ("Grid<Vec3>","getMinValue",Grid<Vec3>::_W_17); 
 static const Pb::Register _R_65 ("Grid<Vec3>","printGrid",Grid<Vec3>::_W_18); 
#endif
#ifdef _C_GridBase
 static const Pb::Register _R_66 ("GridBase","GridBase","PbClass"); template<> const char* Namify<GridBase >::S = "GridBase"; 
 static const Pb::Register _R_67 ("GridBase","GridBase",GridBase::_W_0); 
#endif
#ifdef _C_MACGrid
 static const Pb::Register _R_68 ("MACGrid","MACGrid","Grid<Vec3>"); template<> const char* Namify<MACGrid >::S = "MACGrid"; 
 static const Pb::Register _R_69 ("MACGrid","MACGrid",MACGrid::_W_19); 
#endif
static const Pb::Register _R_0 ("Grid<int>","IntGrid","");
static const Pb::Register _R_1 ("Grid<Real>","RealGrid","");
static const Pb::Register _R_2 ("Grid<Vec3>","VecGrid","");
}