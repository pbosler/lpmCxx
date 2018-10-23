#include "LpmAosEdge.hpp"
#include "LpmAosParticleFactory.hpp"
#include "LpmGll.hpp"
#include <sstream>
#include <iostream>
#include <iomanip>

namespace Lpm {
namespace Aos {

template <int ndim> KidEdgeArrays<ndim> Edge<ndim>::divide(ParticleSet<ndim>& particles, const scalar_type radius, const GeometryType geom) const {
	KidEdgeArrays<ndim> result;
	// find edge midpoints in physical and material space
	Vec<ndim> midpt;
	Vec<ndim> lagmidpt;
	if (geom == SPHERICAL_SURFACE_GEOMETRY) {
		midpt = this->sphMidpoint(particles, radius);
		lagmidpt = this->sphLagMidpoint(particles, radius);
	}
	else if (geom == PLANAR_GEOMETRY or geom == CARTESIAN_3D_GEOMETRY) {
		midpt = this->midpoint(particles);
		lagmidpt = this->lagMidpoint(particles);
	}
	const index_type particle_insert_point = particles.n();
	result.newOrigs[0] = this->_orig;
	result.newDests[0] = particle_insert_point;
	result.newLefts[0] = this->_left;
	result.newRights[0] = this->_right;
	
	result.newOrigs[1] = particle_insert_point;
	result.newDests[1] = this->_dest;
	result.newLefts[1] = this->_left;
	result.newRights[1] = this->_right;
	
	particles.insert(midpt, lagmidpt, 0.0, true);
	return result;
}

template <int ndim> std::string KidEdgeArrays<ndim>::infoString() const {
	std::ostringstream ss;
	std::cout << std::setw(10) << "origs" << std::setw(10) << "dest" << std::setw(10) << "lefts" 
		<< std::setw(10) << "rights" << std::setw(10) << "mid0" << std::setw(10) << "mid1" << std::endl;
	for (int i=0; i<newOrigs.size(); ++i) {
		std::cout << std::setw(10) << newOrigs[i] << std::setw(10) << newDests[i] << std::setw(10) << newLefts[i] 
			<< std::setw(10) << newRights[i] << std::setw(10) << newMids[i][0] << std::setw(10) << newMids[i][1] << std::endl;
	}
	return ss.str();
}

template <int ndim> KidEdgeArrays<ndim> QuadraticEdge<ndim>::divide(ParticleSet<ndim>& particles, const scalar_type radius, const GeometryType geom) const {
	KidEdgeArrays<ndim> result;
	const Vec<ndim> origcrd = this->origCrd(particles);
	const Vec<ndim> origlagcrd = this->origLagCrd(particles);
	const Vec<ndim> midcrd = this->mid0Crd(particles);
	const Vec<ndim> midlagcrd = this->mid0LagCrd(particles);	
	const Vec<ndim> destcrd = this->destCrd(particles);
	const Vec<ndim> destlagcrd = this->destLagCrd(particles);

	Vec<ndim> kid0midpt;
	Vec<ndim> kid1midpt;
	Vec<ndim> kid0lagmidpt;
	Vec<ndim> kid1lagmidpt;
	if (geom == SPHERICAL_SURFACE_GEOMETRY) {
		kid0midpt = origcrd.sphereMidpoint(midcrd, radius);
		kid0lagmidpt = origlagcrd.sphereMidpoint(midlagcrd, radius);
		kid1midpt = midcrd.sphereMidpoint(destcrd, radius);
		kid1lagmidpt = midlagcrd.sphereMidpoint(destlagcrd, radius);
	}
	else if (geom == PLANAR_GEOMETRY or geom == CARTESIAN_3D_GEOMETRY) {
		kid0midpt = origcrd.midpoint(midcrd);
		kid0lagmidpt = origlagcrd.midpoint(midlagcrd);
		kid1midpt = midcrd.midpoint(destcrd);
		kid1lagmidpt = midlagcrd.midpoint(destlagcrd);
	}
	const index_type particle_insert_point = particles.n();
	result.newOrigs[0] = this->_orig;
	result.newDests[0] = this->_midpts[0];
	result.newLefts[0] = this->_left;
	result.newRights[0] = this->_right;
	result.newMids[0][0] = particle_insert_point;
	result.newOrigs[1] = this->_midpts[0];
	result.newDests[1] = this->_dest;
	result.newLefts[1] = this->_left;
	result.newRights[1] = this->_right;
	result.newMids[1][0] = particle_insert_point + 1;
	
	particles.insert(kid0midpt, kid0lagmidpt, 0.0, true);
	particles.insert(kid1midpt, kid1lagmidpt, 0.0, true);
		
// 	result.particles.push_back(particles.getFactory()->createParticle(kid0midpt, kid0lagmidpt));
// 	result.particles.push_back(particles.getFactory()->createParticle(kid1midpt, kid1lagmidpt));
// 	result.newPhysCrds[0] = kid0midpt;
// 	result.newLagCrds[0] = kid0lagmidpt;
// 	result.newPhysCrds[1] = kid1midpt;
// 	result.newLagCrds[1] = kid1lagmidpt;
	return result;	
}

template <int ndim> KidEdgeArrays<ndim> CubicEdge<ndim>::divide(ParticleSet<ndim>& particles, const scalar_type radius, const GeometryType geom) const {
	KidEdgeArrays<ndim> result;
	const Vec<ndim> origcrd = this->origCrd(particles);
	const Vec<ndim> origlagcrd = this->origLagCrd(particles);
// 	const Vec<ndim> mid0crd = this->mid0Crd(particles);
// 	const Vec<ndim> mid0lagcrd = this->mid0LagCrd(particles);	
// 	const Vec<ndim> mid1crd = this->mid1Crd(particles);
// 	const Vec<ndim> mid1lagcrd = this->mid1LagCrd(particles);		
	const Vec<ndim> destcrd = this->destCrd(particles);
	const Vec<ndim> destlagcrd = this->destLagCrd(particles);
	
	Vec<ndim> parentmidpt;
	Vec<ndim> parentlagmidpt;
	Vec<ndim> kid0mid0crd;
	Vec<ndim> kid0mid0lagcrd;
	Vec<ndim> kid0mid1crd;
	Vec<ndim> kid0mid1lagcrd;
	Vec<ndim> kid1mid0crd;
	Vec<ndim> kid1mid0lagcrd;
	Vec<ndim> kid1mid1crd;
	Vec<ndim> kid1mid1lagcrd;
	if (geom == SPHERICAL_SURFACE_GEOMETRY) {
		parentmidpt = this->sphMidpoint(particles, radius);
		parentlagmidpt = this->sphLagMidpoint(particles, radius);
		kid0mid0crd = pointAlongCircle(origcrd, parentmidpt, -CubicGLL<ndim>::sqrt5);
		kid0mid0lagcrd = pointAlongCircle(origlagcrd, parentlagmidpt, -CubicGLL<ndim>::sqrt5);
		kid0mid1crd = pointAlongCircle(origcrd, parentmidpt, CubicGLL<ndim>::sqrt5);
		kid0mid1lagcrd = pointAlongCircle(origlagcrd, parentlagmidpt, CubicGLL<ndim>::sqrt5);
		kid1mid0crd = pointAlongCircle(parentmidpt, destcrd, -CubicGLL<ndim>::sqrt5);
		kid1mid0lagcrd = pointAlongCircle(parentlagmidpt, destlagcrd, -CubicGLL<ndim>::sqrt5);
		kid1mid1crd = pointAlongCircle(parentmidpt, destcrd, CubicGLL<ndim>::sqrt5);
		kid1mid1lagcrd = pointAlongCircle(parentlagmidpt, destlagcrd, CubicGLL<ndim>::sqrt5);		
	}
	else if (geom == PLANAR_GEOMETRY or geom == CARTESIAN_3D_GEOMETRY) {
		parentmidpt = this->midpoint(particles);
		parentlagmidpt = this->lagMidpoint(particles);
		kid0mid0crd = pointAlongChord(origcrd, parentmidpt, -CubicGLL<ndim>::sqrt5);
		kid0mid0lagcrd = pointAlongChord(origlagcrd, parentlagmidpt, -CubicGLL<ndim>::sqrt5);
		kid0mid1crd = pointAlongChord(origcrd, parentmidpt, CubicGLL<ndim>::sqrt5);
		kid0mid1lagcrd = pointAlongChord(origlagcrd, parentlagmidpt, CubicGLL<ndim>::sqrt5);
		kid1mid0crd = pointAlongChord(parentmidpt, destcrd, -CubicGLL<ndim>::sqrt5);
		kid1mid0lagcrd = pointAlongChord(parentlagmidpt, destlagcrd, -CubicGLL<ndim>::sqrt5);
		kid1mid1crd = pointAlongChord(parentmidpt, destcrd, CubicGLL<ndim>::sqrt5);
		kid1mid1lagcrd = pointAlongChord(parentlagmidpt, destlagcrd, CubicGLL<ndim>::sqrt5);
	}
// 	result.newPhysCrds[0] = kid0mid0crd;
// 	result.newLagCrds[0] = kid0mid0lagcrd;
// 	result.newPhysCrds[1] = kid0mid1crd;
// 	result.newLagCrds[1] = kid0mid1lagcrd;
// 	result.newPhysCrds[2] = parentmidpt;
// 	result.newLagCrds[2] = parentlagmidpt;
// 	result.newPhysCrds[3] = kid1mid0crd;
// 	result.newLagCrds[3] = kid1mid0lagcrd;
// 	result.newPhysCrds[4] = kid1mid1crd;
// 	result.newLagCrds[4] = kid1mid1lagcrd;
	const index_type particle_insert_point = particles.n();
	particles.move(this->_midpts[0], kid0mid0crd, kid0mid0lagcrd);
	particles.insert(kid0mid1crd, kid0mid1lagcrd, 0.0, true);
	particles.insert(parentmidpt, parentlagmidpt, 0.0, true);
	particles.insert(kid1mid0crd, kid1mid0lagcrd, 0.0, true);
	particles.move(this->_midpts[1], kid1mid1crd, kid1mid1lagcrd);
	
	result.newOrigs[0] = this->_orig;
	result.newMids[0][0] = this->_midpts[0];
	result.newMids[0][1] = particle_insert_point;
	result.newDests[0] = particle_insert_point+1;
	result.newLefts[0] = this->_left;
	result.newRights[0] = this->_right;
	result.newOrigs[1] = particle_insert_point+1;
	result.newMids[1][0] = particle_insert_point+2;
	result.newMids[1][1] = this->_midpts[1];
	result.newDests[1] = this->_dest;
	result.newLefts[1] = this->_left;
	result.newRights[1] = this->_right;
	return result;
}

template <int ndim> Vec<ndim> Edge<ndim>::midpoint(const ParticleSet<ndim>& particles) const {
    Particle<ndim>* optr = particles.getPtr(_orig);
    Particle<ndim>* dptr = particles.getPtr(_dest);
    return optr->physCrd().midpoint(dptr->physCrd());
}

template <int ndim> Vec<ndim> Edge<ndim>::lagMidpoint(const ParticleSet<ndim>& particles) const {
    Particle<ndim>* optr = particles.getPtr(_orig);
    Particle<ndim>* dptr = particles.getPtr(_dest);
    return optr->lagCrd().midpoint(dptr->lagCrd());
}

template <int ndim> Vec<ndim> Edge<ndim>::sphMidpoint(const ParticleSet<ndim>& particles, const scalar_type radius) const {
    Particle<ndim>* optr = particles.getPtr(_orig);
    Particle<ndim>* dptr = particles.getPtr(_dest);
    return optr->physCrd().sphereMidpoint(dptr->physCrd(), radius);
}

template <int ndim> Vec<ndim> Edge<ndim>::sphLagMidpoint(const ParticleSet<ndim>& particles, const scalar_type radius) const {
    Particle<ndim>* optr = particles.getPtr(_orig);
    Particle<ndim>* dptr = particles.getPtr(_dest);
    return optr->lagCrd().sphereMidpoint(dptr->lagCrd(), radius);
}

template <int ndim> Vec<ndim> Edge<ndim>::edgeVector(const ParticleSet<ndim>& particles) const {
    Particle<ndim>* optr = particles.getPtr(_orig);
    Particle<ndim>* dptr = particles.getPtr(_dest);
    return dptr->physCrd() - optr->physCrd();
}

template <int ndim> scalar_type Edge<ndim>::eucLength(const ParticleSet<ndim>& particles) const {
    Particle<ndim>* optr = particles.getPtr(_orig);
    Particle<ndim>* dptr = particles.getPtr(_dest);
    return optr->physCrd().dist(dptr->physCrd());
}

template <int ndim> scalar_type Edge<ndim>::sphLength(const ParticleSet<ndim>& particles, const scalar_type radius) const {
    Particle<ndim>* optr = particles.getPtr(_orig);
    Particle<ndim>* dptr = particles.getPtr(_dest);
    return optr->physCrd().sphereDist(dptr->physCrd(), radius);
}

template <int ndim> std::string Edge<ndim>::infoString() const {
    std::ostringstream oss;
    oss << "edge info:" << std::endl;
    oss << "\torig = " << _orig << std::endl;
    oss << "\tdest = " << _dest << std::endl;
    oss << "\tleft = " << _left << std::endl;
    oss << "\tright= " << _right << std::endl;
    oss << "\thasKids = " << (this->hasKids() ? "true" : "false") << std::endl;
    oss << "\tkids = " << _kids[0] << " " << _kids[1] << std::endl;
    oss << "\tparent = " << _parent << std::endl;
    oss << "\tinterior_indices = " << _midpts[0] << " " << _midpts[1] << std::endl;
    return oss.str();
}


// template <int ndim> Edge<ndim>::unitTangent(const ParticleSet<ndim>& particles) const {
// 	return edgeVector(particles).normalize();
// }
// 
// template <int ndim> Edge<ndim>::unitNormal(const ParticleSet<ndim>& particles, const GeometryType geom) const {
// 	Vec<ndim> result;
// 	switch (geom) : {
// 		case (PLANAR_GEOMETRY) : {
// 			break;
// 		}
// 		case (SPHERICAL_SURFACE_GEOMETRY) : {
// 			break;
// 		}
// 	}
// 	return result;
// }

template struct KidEdgeArrays<2>;
template struct KidEdgeArrays<3>;
template class Edge<2>;
template class Edge<3>;
template class QuadraticEdge<2>;
template class QuadraticEdge<3>;
template class CubicEdge<2>;
template class CubicEdge<3>;
}
}