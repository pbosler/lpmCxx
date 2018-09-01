#include "LpmAosEdge.hpp"
#include "LpmAosTypes.hpp"
#include "LpmAosParticle.hpp"
#include <sstream>
#include <iostream>

namespace Lpm {
namespace Aos {

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


template class Edge<2>;
template class Edge<3>;
}
}