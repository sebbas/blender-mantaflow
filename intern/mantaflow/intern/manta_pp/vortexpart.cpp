




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/Users/sbarschkis/Developer/Mantaflow/blenderIntegration/mantaflowgit/source/vortexpart.cpp"
/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Vortex particles
 * (warning, the vortex methods are currently experimental, and not fully supported!)
 *
 ******************************************************************************/

#include "vortexpart.h"
#include "integrator.h"
#include "mesh.h"

using namespace std;
namespace Manta {

// vortex particle effect: (cyl coord around wp)
// u = -|wp|*rho*exp( (-rho^2-z^2)/(2sigma^2) ) e_phi
inline Vec3 VortexKernel(const Vec3& p, const vector<VortexParticleData>& vp, Real scale) {
	Vec3 u(0.0);
	for (size_t i=0; i<vp.size(); i++) {
		if (vp[i].flag & ParticleBase::PDELETE) continue;
		
		// cutoff radius
		const Vec3 r = p - vp[i].pos;
		const Real rlen2 = normSquare(r);   
		const Real sigma2 = square(vp[i].sigma);
		if (rlen2 > 6.0 * sigma2 || rlen2 < 1e-8) continue;
		
		// split vortex strength
		Vec3 vortNorm = vp[i].vorticity;
		Real strength = normalize(vortNorm) * scale;
	
		// transform in cylinder coordinate system
		const Real rlen = sqrt(rlen2);
		const Real z = dot(r, vortNorm);
		const Vec3 ePhi = cross(r, vortNorm) / rlen;
		const Real rho2 = rlen2 - z*z;
	
		Real vortex = 0;
		if (rho2 > 1e-10) {
			// evaluate Kernel      
			vortex = strength * sqrt(rho2) * exp (rlen2 * -0.5/sigma2);  
		}
		u += vortex * ePhi;
	}
	return u;
}


 struct KnVpAdvectMesh : public KernelBase { KnVpAdvectMesh(vector<Node>& nodes, const vector<VortexParticleData>& vp, Real scale) :  KernelBase(nodes.size()) ,nodes(nodes),vp(vp),scale(scale) ,u((size))  { runMessage(); run(); }   inline void op(IndexInt idx, vector<Node>& nodes, const vector<VortexParticleData>& vp, Real scale ,vector<Vec3> & u)  {
	if (nodes[idx].flags & Mesh::NfFixed)
		u[idx] = 0.0;
	else
		u[idx] = VortexKernel(nodes[idx].pos, vp, scale);
}    inline operator vector<Vec3> () { return u; } inline vector<Vec3>  & getRet() { return u; }  inline vector<Node>& getArg0() { return nodes; } typedef vector<Node> type0;inline const vector<VortexParticleData>& getArg1() { return vp; } typedef vector<VortexParticleData> type1;inline Real& getArg2() { return scale; } typedef Real type2; void runMessage() { debMsg("Executing kernel KnVpAdvectMesh ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (IndexInt i = 0; i < _sz; i++) op(i,nodes,vp,scale,u);  }   } vector<Node>& nodes; const vector<VortexParticleData>& vp; Real scale;  vector<Vec3>  u;  };
#line 56 "vortexpart.cpp"




 struct KnVpAdvectSelf : public KernelBase { KnVpAdvectSelf(vector<VortexParticleData>& vp, Real scale) :  KernelBase(vp.size()) ,vp(vp),scale(scale) ,u((size))  { runMessage(); run(); }   inline void op(IndexInt idx, vector<VortexParticleData>& vp, Real scale ,vector<Vec3> & u)  {
	if (vp[idx].flag & ParticleBase::PDELETE) 
		u[idx] = 0.0;
	else
		u[idx] = VortexKernel(vp[idx].pos, vp, scale);
}    inline operator vector<Vec3> () { return u; } inline vector<Vec3>  & getRet() { return u; }  inline vector<VortexParticleData>& getArg0() { return vp; } typedef vector<VortexParticleData> type0;inline Real& getArg1() { return scale; } typedef Real type1; void runMessage() { debMsg("Executing kernel KnVpAdvectSelf ", 2); debMsg("Kernel range" << " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 3); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for 
  for (IndexInt i = 0; i < _sz; i++) op(i,vp,scale,u);  }   } vector<VortexParticleData>& vp; Real scale;  vector<Vec3>  u;  };
#line 64 "vortexpart.cpp"


	
VortexParticleSystem::VortexParticleSystem(FluidSolver* parent) :
	ParticleSystem<VortexParticleData>(parent)
{ 
}

void VortexParticleSystem::advectSelf(Real scale, int integrationMode) {
	KnVpAdvectSelf kernel(mData, scale* getParent()->getDt());
	integratePointSet( kernel, integrationMode);    
}

void VortexParticleSystem::applyToMesh(Mesh& mesh, Real scale, int integrationMode) {
	KnVpAdvectMesh kernel(mesh.getNodeData(), mData, scale* getParent()->getDt());
	integratePointSet( kernel, integrationMode);    
}

ParticleBase* VortexParticleSystem::clone() {
	VortexParticleSystem* nm = new VortexParticleSystem(getParent());
	compress();
	
	nm->mData = mData;
	nm->setName(getName());
	return nm;
}

	

} // namespace


