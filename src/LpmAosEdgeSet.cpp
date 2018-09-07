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
    const std::vector<index_type>& mids) {
    if (_edges.size() + 1 > _nMax) {
        throw std::out_of_range("EdgeSet::insert _nMax exceeded.");
    }
    _edges.push_back(_factory->createEdge(origID, destID, leftID, rightID, mids));
    _nActive +=1;
}

template <int ndim> void EdgeSet<ndim>::divide(const index_type ind, ParticleSet<ndim>& particles, const scalar_type radius) {
    if (_edges.size() + 2 > _nMax) {
        throw std::out_of_range("EdgeSet::divide _nMax exceeded.");
    }
	const KidEdgeArrays<ndim> kids = _edges[ind]->divide(particles, radius, _geom);
	const index_type edge_insert_pt = _edges.size();
	if (dynamic_cast<QuadraticEdge<ndim>*>(_edges[0].get())) {
		for (short i=0; i<2; ++i) {
			particles.insert(kids.newPhysCrds[i], kids.newLagCrds[i]);
		}
    }
    else if (dynamic_cast<CubicEdge<ndim>*>(_edges[0].get())) {
    	particles.move(_edges[ind]->midpt0(), kids.newPhysCrds[0], kids.newLagCrds[0]);
    	for (short i=1; i<=3; ++i) {
	    	particles.insert(kids.newPhysCrds[i], kids.newLagCrds[i]);
	    }
    	particles.move(_edges[ind]->midpt1(), kids.newPhysCrds[4], kids.newLagCrds[4]);     
    }
    else {
    	particles.insert(kids.newPhysCrds[0], kids.newLagCrds[0]);
    }
    for (short i=0; i<2; ++i) {
		this->insert(kids.newOrigs[i], kids.newDests[i], kids.newLefts[i], kids.newRights[i], kids.newMids[i]);
		_edges[edge_insert_pt+i]->setParent(ind);
	}
    _edges[ind]->setKids(edge_insert_pt, edge_insert_pt+1);
    _nActive -= 1; // parent edge is no longer active.
}

template <int ndim> std::string EdgeSet<ndim>::infoString(const bool printAll) const {
    std::ostringstream ss;
    ss << "EdgeSet info:" << std::endl;
    ss << "\tgeom = " << geometryString(_geom) << std::endl;
    ss << "\tnMax = " << _nMax << std::endl;
    ss << "\tnActive = " << _nActive << std::endl;
    if (printAll) {
        for (index_type i=0; i<_edges.size(); ++i) {
            ss << i << ": " << _edges[i]->infoString();
        }
    }
    return ss.str();
}

#ifdef HAVE_VTK
template <int ndim> vtkSmartPointer<vtkCellArray> EdgeSet<ndim>::toVtkCellArray() const {
	vtkSmartPointer<vtkCellArray> result = vtkSmartPointer<vtkCellArray>::New();
	const index_type ptsPerEdge = _edges[0]->ptsPerEdge();
	std::array<index_type,2> midpts;
	for (index_type i=0; i<_edges.size(); ++i) {
		if (_edges[i]->isLeaf()) {
			result->InsertNextCell(ptsPerEdge);
			result->InsertCellPoint(_edges[i]->orig());
			midpts = _edges[i]->midpts();
			if (midpts[0] >= 0) result->InsertCellPoint(midpts[0]);
			if (midpts[1] >= 0) result->InsertCellPoint(midpts[1]);
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

