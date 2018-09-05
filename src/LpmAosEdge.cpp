#include "LpmAosEdge.hpp"
#include "LpmAosParticleFactory.hpp"
#include <sstream>
#include <iostream>
#include <iomanip>

namespace Lpm {
namespace Aos {

template <int ndim> KidEdgeArrays<ndim> Edge<ndim>::divide(const ParticleSet<ndim>& particles, const scalar_type radius, const GeometryType geom) const {
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
	
	result.particles.push_back(particles.getFactory()->createParticle(midpt, lagmidpt));
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
	for (int i=0; i<particles.size(); ++i) {
		std::cout << particles[i]->infoString();
	}
	return ss.str();
}

template <int ndim> KidEdgeArrays<ndim> QuadraticEdge<ndim>::divide(const ParticleSet<ndim>& particles, const scalar_type radius, const GeometryType geom) const {
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
		kid0midpt = sphereMidpoint(origcrd, midcrd, radius);
		kid0lagmidpt = sphereMidpoint(origlagcrd, midlagcrd, radius);
		kid1midpt = sphereMidpoint(midcrd, destcrd, radius);
		kid1lagmidpt = sphereMidpoint(midlagcrd, destlagcrd, radius);
	}
	else if (geom == PLANAR_GEOMETRY or geom == CARTESIAN_3D_GEOMETRY) {
		kid0midpt = midpoint(origcrd, midcrd);
		kid0lagmidpt = midpoint(origlagcrd, midlagcrd);
		kid1midpt = midpoint(midcrd, destcrd);
		kid1lagmidpt= midpoint(midlagcrd, destlagcrd);
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
	
	result.particles.push_back(particles.getFactory()->createParticle(kid0midpt, kid0lagmidpt));
	result.particles.push_back(particles.getFactory()->createParticle(kid1midpt, kid1lagmidpt));
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
}
}