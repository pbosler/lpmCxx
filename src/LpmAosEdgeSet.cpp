#include "LpmAosEdgeSet.hpp"
#include <limits>

namespace Lpm {

template <int ndim> scalar_type EdgeSet::maxEucLength(const ParticleSet<ndim>& particles) const {
    scalar_type result(0.0);
    for (index_type i=0; i<_edges.size(); ++i) {
        if (_edges[i]->isLeaf() && _edges[i]->eucLength(particles) > result)
            result = _edges[i]->eucLength(particles);
    }
    return result;
}

template <int ndim> scalar_type EdgeSet::minEucLength(const ParticleSet<ndim>& particles) const {
    scalar_type result = std::numeric_limits<scalar_type>::max();
    for (index_type i=0; i<_edges.size(); ++i) {
        if (_edges[i]->isLeaf() && _edges[i]->eucLength(particles) < result)
            result = _edges[i]->eucLength(particles);
    }
    return result;
}

template <int ndim> scalar_type EdgeSet::maxSphLength(const ParticleSet<ndim>& particles, const scalar_type radius) const {
    scalar_type result(0.0);
    for (index_type i=0; i<_edges.size(); ++i) {
        if (_edges[i]->isLeaf() && _edges[i]->sphLength(particles) > result)
            result = _edges[i]->sphLength(particles, radius);
    }
    return result;
}

template <int ndim> scalar_type EdgeSet::minSphLength(const ParticleSet<ndim>& particles, const scalar_type radius) const {
    scalar_type result = std::numeric_limits<scalar_type>::max();
    for (index_type i=0; i<_edges.size(); ++i) {
        if (_edges[i]->isLeaf() && _edges[i]->sphLength(particles) < result)
            result = _edges[i]->sphLength(particles, radius);
    }
    return result;
}

void EdgeSet::insert(const index_type origID, const index_type destID, const index_type leftID, const index_type rightID) {
    _edges.push_back(_factory->createEdge(const index_type origID, const index_type destID,
        const index_type leftID, const index_type rightID));
    _nActive +=1;
}

template <int ndim> void EdgeSet::divide(const index_type ind, ParticleSet<ndim>& particles, const scalar_type radius) {
    if (dynamic_cast<QuadraticEdge*>(_edges[0].get())) {
        //TODO
    }
    else if (dynamic_cast<CubicEdge*>(_edges[0].get())) {
        const Vec<ndim> physOrig = particles.getPtr(_edges[ind]->orig())->physCrd();
        const Vec<ndim> lagOrig = particles.getPtr(_edges[ind]->orig())->lagCrd();
        const Vec<ndim> physDest = particles.getPtr(_edges[ind]->dest())->physCrd();
        const Vec<ndim> lagDest = particles.getPtr(_edges[ind]->dest())->lagCrd();

        // TODO
    }
    else {
        Vec<ndim> midpt;
        Vec<ndim> lagMidpt;
        const index_type particle_insert = particles.n();
        const index_type edge_insert = _edges.size();
        if (_geom == SPHERICAL_SURFACE_GEOMETRY) {
            midpt = _edges[ind]->sphMidpoint(particles, radius);
            lagMidpt = _edges[ind]->sphLagMidpoint(particles, radius);
        }
        else if (_geom == PLANAR_GEOMETRY || _geom == CARTESIAN_3D_GEOMETRY) {
            midpt = _edges[ind]->midpoint(particles);
            lagMidpt = _edges[ind]->lagMidpoint(particles);
        }
        const index_type lface = _edges[ind]->left();
        const index_type rface = _edges[ind]->right();
        particles->insert(midpt, lagMidpt);
        this->insert(_edges[ind]->orig(), particle_insert, lface, rface);
        this->insert(particle_insert, _edges[ind]->dest(), lface, rface);
        _edges[ind]->setKids(edge_insert, edge_insert+1);
        _edges[edge_insert]->setParent(ind);
        _edges[edge_insert+1]->setParent(ind);
    }
    _nActive += 1;
}

}

