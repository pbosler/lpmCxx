#include <iostream>

#include <sstream>
#include "LpmAosFace.hpp"
#include "LpmAosEdge.hpp"
#include "LpmGll.hpp"
#include "LpmUtilities.h"

namespace Lpm {
namespace Aos {

template <int ndim> std::string KidFaceArrays<ndim>::infoString() const {
	std::ostringstream ss;
	ss << "new face data:"<< std::endl;
	for (short i=0; i<4; ++i) {
		ss << "Face " << i << ":" << std::endl;
		ss << "vertices/edge particles: ";
		for (short j=0; j<this->newFaceVerts[0].size(); ++j)
			ss << newFaceVerts[i][j] << " ";
		ss << std::endl << "edges: ";
		for (short j=0; j<this->newFaceEdges[0].size(); ++j) 
			ss << newFaceEdges[i][j] << " ";
		ss << std::endl << "interior particles: ";
		for (short j=0; j<this->newFaceInteriors[0].size(); ++j)
			ss << newFaceInteriors[i][j] << " ";
		ss << std::endl << "child face area: ";
		ss << kidsFaceArea[i];
		ss << std::endl;
	}
	return ss.str();
}

template <int ndim> void Face<ndim>::setArea(const scalar_type ar, ParticleSet<ndim>& particles) {
	_area = ar;
	particles.setWeight(_interiorInds[0], ar);
}

template <int ndim> void Face<ndim>::setArea(ParticleSet<ndim>& particles, const GeometryType geom, const scalar_type radius) {
	_area = this->computeAreaFromCorners(particles, geom, radius);
	particles.setWeight(_interiorInds[0], _area);
}

template <int ndim>  std::vector<Vec<ndim>> Face<ndim>::vertPhysCrds(const ParticleSet<ndim>& particles) const {
	std::vector<Vec<ndim>> result(_vertInds.size());
	for (index_type j=0; j<_vertInds.size(); ++j) {
		result[j] = particles.physCrd(_vertInds[j]);
	}
	return result;
}

template <int ndim> Vec<ndim> Face<ndim>::physBarycenter(const ParticleSet<ndim>& particles) const {
    Vec<ndim> result;
    for (int i=0; i<_vertInds.size(); ++i) {
        result += particles.physCrd(_vertInds[i]);
    }
    result.scaleInPlace(1.0/_vertInds.size());
    return result;
}

template <int ndim> Vec<ndim> Face<ndim>::lagBarycenter(const ParticleSet<ndim>& particles) const {
    Vec<ndim> result;
    for (int i=0; i<_vertInds.size(); ++i) {
        result += particles.lagCrd(_vertInds[i]);
    }
    result.scaleInPlace(1.0/_vertInds.size());
    return result;
}

template <int ndim> Vec<3> Face<ndim>::physSphBarycenter(const ParticleSet<ndim>& particles, 
    const scalar_type radius) const {
    Vec<3> result;
    for (int i=0; i<_vertInds.size(); ++i) {
        result += particles.physCrd(_vertInds[i]);
    }
    result.scaleInPlace(1.0/_vertInds.size());
    result.scaleInPlace(radius);
    return result;
}

template <int ndim> Vec<3> Face<ndim>::lagSphBarycenter(const ParticleSet<ndim>& particles, 
    const scalar_type radius) const {
    Vec<3> result;
    for (int i=0; i<_vertInds.size(); ++i) {    
        result += particles.lagCrd(_vertInds[i]);
    }
    result.scaleInPlace(1.0/_vertInds.size());
    result.scaleInPlace(radius);
    return result;
}

template <int ndim> scalar_type Face<ndim>::scalarIntegral(const std::string field_name, const ParticleSet<ndim>& particles) const {
	scalar_type result = 0.0;
	for (int i=0; i<_vertInds.size(); ++i) {
		result += particles.scalarVal(_vertInds[i], field_name)*particles.weight(_vertInds[i]);
	}
	for (int i=0; i<_interiorInds.size(); ++i) {
		result += particles.scalarVal(_interiorInds[i], field_name)*particles.weight(_interiorInds[i]);
	}
	return result;
}

template <int ndim> std::string Face<ndim>::infoString() const {
    std::ostringstream ss;
    ss << "Face info:" << std::endl;
    ss << "\tvertInds = ";
    for (int i=0; i<_vertInds.size(); ++i) 
        ss << _vertInds[i] << " ";
    ss << std::endl;
    ss << "\tedgeInds = ";
    for (int i=0; i<_edgeInds.size(); ++i)
        ss << _edgeInds[i] << " ";
    ss << std::endl;
    ss << "\tinterior indices = " ;
    for (int i=0; i<_interiorInds.size(); ++i) 
        ss << _interiorInds[i] << " ";
    ss << std::endl;
    ss << "\tparent face = " << _parent << std::endl;
    ss << "\tkids = ";
    for (int i=0; i<4; ++i)
        ss << _kids[i] << " ";
    ss << std::endl;
    ss << "\tarea = " << _area << std::endl;
    return ss.str();
}



template class Face<2>;
template class Face<3>;
template struct KidFaceArrays<2>;
template struct KidFaceArrays<3>;
}
}
