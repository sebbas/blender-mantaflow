




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep link).




#include "shapes.h"
namespace Manta {
#ifdef _C_Box
 static const Pb::Register _R_0 ("Box","Box","Shape"); template<> const char* Namify<Box >::S = "Box"; 
 static const Pb::Register _R_1 ("Box","Box",Box::_W_9); 
#endif
#ifdef _C_Cylinder
 static const Pb::Register _R_2 ("Cylinder","Cylinder","Shape"); template<> const char* Namify<Cylinder >::S = "Cylinder"; 
 static const Pb::Register _R_3 ("Cylinder","Cylinder",Cylinder::_W_11); 
 static const Pb::Register _R_4 ("Cylinder","setRadius",Cylinder::_W_12); 
 static const Pb::Register _R_5 ("Cylinder","setZ",Cylinder::_W_13); 
#endif
#ifdef _C_NullShape
 static const Pb::Register _R_6 ("NullShape","NullShape","Shape"); template<> const char* Namify<NullShape >::S = "NullShape"; 
 static const Pb::Register _R_7 ("NullShape","NullShape",NullShape::_W_8); 
#endif
#ifdef _C_Shape
 static const Pb::Register _R_8 ("Shape","Shape","PbClass"); template<> const char* Namify<Shape >::S = "Shape"; 
 static const Pb::Register _R_9 ("Shape","Shape",Shape::_W_0); 
 static const Pb::Register _R_10 ("Shape","applyToGrid",Shape::_W_1); 
 static const Pb::Register _R_11 ("Shape","applyToGridSmooth",Shape::_W_2); 
 static const Pb::Register _R_12 ("Shape","computeLevelset",Shape::_W_3); 
 static const Pb::Register _R_13 ("Shape","collideMesh",Shape::_W_4); 
 static const Pb::Register _R_14 ("Shape","getCenter",Shape::_W_5); 
 static const Pb::Register _R_15 ("Shape","setCenter",Shape::_W_6); 
 static const Pb::Register _R_16 ("Shape","getExtent",Shape::_W_7); 
#endif
#ifdef _C_Slope
 static const Pb::Register _R_17 ("Slope","Slope","Shape"); template<> const char* Namify<Slope >::S = "Slope"; 
 static const Pb::Register _R_18 ("Slope","Slope",Slope::_W_14); 
#endif
#ifdef _C_Sphere
 static const Pb::Register _R_19 ("Sphere","Sphere","Shape"); template<> const char* Namify<Sphere >::S = "Sphere"; 
 static const Pb::Register _R_20 ("Sphere","Sphere",Sphere::_W_10); 
#endif
}