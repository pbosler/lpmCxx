#include "LpmAosEdgeSet.hpp"
#include <limits>
#include <exception>
#include <sstream>
#include <iostream>

namespace Lpm {
namespace Aos {

template <int ndim> scalar_type EdgeSet<ndim>::maxEucLength(const ParticleSet<ndim>& particles) const {
    scalar_type result(0.0);
    for (index_type i=0; i<_edges.size(); ++i) {
        if (_edges[i]->isLeaf() && _edges[i]->eucLength(particles) > result)
            result = _edges[i]->eucLength(particles);
    }
    return result;
}

template <int ndim> scalar_type EdgeSet<ndim>::minEucLength(const ParticleSet<ndim>& particles) const {
    scalar_type result = std::numeric_limits<scalar_type>::max();
    for (index_type i=0; i<_edges.size(); ++i) {
        if (_edges[i]->isLeaf() && _edges[i]->eucLength(particles) < result)
            result = _edges[i]->eucLength(particles);
    }
    return result;
}

template <int ndim> scalar_type EdgeSet<ndim>::maxSphLength(const ParticleSet<ndim>& particles, const scalar_type radius) const {
    scalar_type result(0.0);
    for (index_type i=0; i<_edges.size(); ++i) {
        if (_edges[i]->isLeaf() && _edges[i]->sphLength(particles) > result)
            result = _edges[i]->sphLength(particles, radius);
    }
    return result;
}

template <int ndim> scalar_type EdgeSet<ndim>::minSphLength(const ParticleSet<ndim>& particles, const scalar_type radius) const {
    scalar_type result = std::numeric_limits<scalar_type>::max();
    for (index_type i=0; i<_edges.size(); ++i) {
        if (_edges[i]->isLeaf() && _edges[i]->sphLength(particles) < result)
            result = _edges[i]->sphLength(particles, radius);
    }
    return result;
}

template <int ndim> void EdgeSet<ndim>::insert(const index_type origID, const index_type destID, const index_type leftID, const index_type rightID, 
    const std::array<index_type,2>& mids) {
    if (_edges.size() + 1 >= _nMax) {
        throw std::out_of_range("EdgeSet::insert _nMax exceeded.");
    }
    _edges.push_back(_factory->createEdge(origID, destID, leftID, rightID, mids));
    _nActive +=1;
}

template <int ndim> void EdgeSet<ndim>::divide(const index_type ind, ParticleSet<ndim>& particles, const scalar_type radius) {
    if (_edges.size() + 2 >= _nMax) {
        throw std::out_of_range("EdgeSet::divide _nMax exceeded.");
    }
    
    Vec<ndim> midpt;
    Vec<ndim> lagMidpt;
    if (_geom == SPHERICAL_SURFACE_GEOMETRY) {
        midpt = _edges[ind]->sphMidpoint(particles, radius);
        lagMidpt = _edges[ind]->sphLagMidpoint(particles, radius);
    }
    else if (_geom == PLANAR_GEOMETRY || _geom == CARTESIAN_3D_GEOMETRY) {
        midpt = _edges[ind]->midpoint(particles);
        lagMidpt = _edges[ind]->lagMidpoint(particles);
    }
    
    std::array<Vec<ndim>, 4> new_crds;
    std::array<Vec<ndim>, 4> new_lag_crds;
    if (dynamic_cast<QuadraticEdge<ndim>*>(_edges[0].get())) {
        //TODO
        throw std::runtime_error("EdgeSet::divide ERROR: quadratic edges not implemented.");
    }
    else if (dynamic_cast<CubicEdge<ndim>*>(_edges[0].get())) {
        //TODO
        throw std::runtime_error("EdgeSet::divide ERROR: cubic edges not implemented.");
    
        const Vec<ndim> physOrig = particles.getPtr(_edges[ind]->orig())->physCrd();
        const Vec<ndim> lagOrig = particles.getPtr(_edges[ind]->orig())->lagCrd();
        const Vec<ndim> physDest = particles.getPtr(_edges[ind]->dest())->physCrd();
        const Vec<ndim> lagDest = particles.getPtr(_edges[ind]->dest())->lagCrd();
        
        if (_geom == SPHERICAL_SURFACE_GEOMETRY) {
            
        }
        else if (_geom == PLANAR_GEOMETRY || _geom == CARTESIAN_3D_GEOMETRY) {
        }
        
    }
    else {
        Vec<ndim> midpt;
        Vec<ndim> lagMidpt;
        const index_type particle_insert = particles.n();
        const index_type edge_insert = _edges.size();
        
        const index_type lface = _edges[ind]->left();
        const index_type rface = _edges[ind]->right();
        particles.insert(midpt, lagMidpt);
        this->insert(_edges[ind]->orig(), particle_insert, lface, rface);
        this->insert(particle_insert, _edges[ind]->dest(), lface, rface);
        _edges[ind]->setKids(edge_insert, edge_insert+1);
        _edges[edge_insert]->setParent(ind);
        _edges[edge_insert+1]->setParent(ind);
    }
    _nActive -= 1; // parent edge is no longer active.
}

// template <int ndim> void EdgeSet<ndim>::initFromVtkPolydataFile(const std::string& fname) {
//     std::string fullfname(LPM_PARTICLE_SET_DIR);
//     fullfname += "/" + fname;
//     std::ifstream file(fullfname);
//     if (!file.is_open()) {
//         throw std::ios_base::failure("EdgeSet::initFromParticleSetFile ERROR: file " + fullfname + " not found.");
//     }
//     
//     index_type lineNumber = 0;
//     std::string line;
//     
// }

template <int ndim> std::string EdgeSet<ndim>::infoString(const bool printAll) const {
    std::ostringstream ss;
    ss << "EdgeSet info:" << std::endl;
    ss << "\tgeom = " << geometryString(_geom) << std::endl;
    ss << "\tnMax = " << _nMax << std::endl;
    ss << "\tnActive = " << _nActive << std::endl;
    if (printAll) {
        for (index_type i=0; i<_edges.size(); ++i) {
            ss << _edges[i]->infoString();
        }
    }
    return ss.str();
}

#ifdef HAVE_VTK
template <int ndim> vtkSmartPointer<vtkCellArray> EdgeSet<ndim>::toVtkCellArray() const {
	vtkSmartPointer<vtkCellArray> result = vtkSmartPointer<vtkCellArray>::New();
	const index_type ptsPerEdge = _edges[0]->ptsPerEdge();
	for (index_type i=0; i<_edges.size(); ++i) {
		if (_edges[i]->isLeaf()) {
			result->InsertNextCell(ptsPerEdge);
			result->InsertCellPoint(_edges[i]->orig());
			if (_edges[i]->_midpts[0] >= 0) result->InsertCellPoint(_edges[i]->_midpts[0]);
			if (_edges[i]->_midpts[1] >= 0) result->InsertCellPoint(_edges[i]->_midpts[1]);
			result->InsertCellPoint(_edges[i]->dest());
		}
	}
	
	return result;
}
#endif

template class EdgeSet<2>;
template class EdgeSet<3>;

}
}

