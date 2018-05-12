#include "LpmAosEdge.hpp"
#include "LpmAosTypes.hpp"
#include "LpmAosParticle.hpp"
#include <sstream>
#include <iostream>

namespace Lpm {

template <int ndim> Vec<ndim> Edge::midpoint(const ParticleSet<ndim>& particles) const {
    Particle<ndim>* optr = particles.getPtr(_orig);
    Particle<ndim>* dptr = particles.getPtr(_dest);
    return optr->physCrd().midpoint(dptr->physCrd());
}

template <int ndim> Vec<ndim> Edge::lagMidpoint(const ParticleSet<ndim>& particles) const {
    Particle<ndim>* optr = particles.getPtr(_orig);
    Particle<ndim>* dptr = particles.getPtr(_dest);
    return optr->lagCrd().midpoint(dptr->lagCrd());
}

template <int ndim> Vec<ndim> Edge::sphMidpoint(const ParticleSet<ndim>& particles, const scalar_type radius) const {
    Particle<ndim>* optr = particles.getPtr(_orig);
    Particle<ndim>* dptr = particles.getPtr(_dest);
    return optr->physCrd().sphMidpoint(dptr->physCrd(), radius);
}

template <int ndim> Vec<ndim> Edge::sphLagMidpoint(const ParticleSet<ndim>& particles, const scalar_type radius) const {
    Particle<ndim>* optr = particles.getPtr(_orig);
    Particle<ndim>* dptr = particles.getPtr(_dest);
    return optr->lagCrd().sphMidpoint(dptr->lagCrd(), radius);
}

template <int ndim> Vec<ndim> Edge::edgeVector(const ParticleSet<ndim>& particles) const {
    Particle<ndim>* optr = particles.getPtr(_orig);
    Particle<ndim>* dptr = particles.getPtr(_dest);
    return dptr->physCrd() - optr->physCrd();
}

template <int ndim> scalar_type Edge::eucLength(const ParticleSet<ndim>& particles) const {
    Particle<ndim>* optr = particles.getPtr(_orig);
    Particle<ndim>* dptr = particles.getPtr(_dest);
    return optr->physCrd().dist(dptr->physCrd());
}

template <int ndim> scalar_type Edge::sphLength(const ParticleSet<ndim>& particles, const scalar_type radius) const {
    Particle<ndim>* optr = particles.getPtr(_orig);
    Particle<ndim>* dptr = particles.getPtr(_dest);
    return optr->physCrd().sphereDist(dptr->physCrd(), radius);
}

std::string Edge::infoString() const {
    std::ostringstream oss;
    oss << "edge info:" << std::endl;
    oss << "\torig = " << _orig << std::endl;
    oss << "\tdest = " << _dest << std::endl;
    oss << "\tleft = " << _left << std::endl;
    oss << "\tright= " << _right << std::endl;
    oss << "\thasKids = " << (this->hasKids() ? "true" : "false") << std::endl;
    oss << "\tkids = " << _kids[0] << " " << _kids[1] << std::endl;
    oss << "\tparent = " << _parent << std::endl;
    return oss.str();
}

}
