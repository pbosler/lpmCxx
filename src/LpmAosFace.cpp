#include <iostream>
#include <sstream>
#include "LpmAosFace.hpp"
#include "LpmAosEdge.hpp"

namespace Lpm {

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
}

