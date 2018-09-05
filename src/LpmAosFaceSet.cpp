#include "LpmAosFaceSet.hpp"
#include <exception>
#include <limits>
#include <iostream>
#include <iomanip>
#include <sstream>

namespace Lpm {
namespace Aos {

template <int ndim> void FaceSet<ndim>::insert(const ind_vec& intrs, const ind_vec& verts, const ind_vec& edges, 
            const index_type pt, const scalar_type ar) {
    if (_faces.size() + 1 > _nMax) {
        throw std::out_of_range("FaceSet::insert _nMax exceeded.");
    }
    _faces.push_back(_factory->createFace(intrs, verts, edges, pt, ar));
    _nActive += 1;
}

template <int ndim> void FaceSet<ndim>::insert(Face<ndim>* fptr) {
    if (_faces.size() + 1 > _nMax) {
        throw std::out_of_range("FaceSet::insert (ptr) _nMax exceeded");
    }
  // _faces.push_back(std::unique_ptr<Face<ndim>>&)
}

template <int ndim> void FaceSet<ndim>::divide(const index_type ind, ParticleSet<ndim>& particles, EdgeSet<ndim>& edges) {
    if (_faces.size() + 4 > _nMax) {
        throw std::out_of_range("FaceSet::divide _nMax exceeded");
    }
    //this->insert(_faces[ind]->divide(particles, edges));
}

template <int ndim> scalar_type FaceSet<ndim>::minArea() const {
    scalar_type result = std::numeric_limits<scalar_type>::max();
    for (index_type i=0; i<_faces.size(); ++i) {
        if (_faces[i]->area() < result)
            result = _faces[i]->area();
    }
    return result;
}

template <int ndim> scalar_type FaceSet<ndim>::maxArea() const {
    scalar_type result = 0.0;
    for (index_type i=0; i<_faces.size(); ++i) {
        if (_faces[i]->area() > result)
            result = _faces[i]->area();
    }
    return result;
}

template <int ndim> scalar_type FaceSet<ndim>::maxLeafArea() const {
    scalar_type result = 0.0;
    for (index_type i=0; i<_faces.size(); ++i) {
        if (_faces[i]->isLeaf()) {
            if (_faces[i]->area() > result) 
                result = _faces[i]->area();
        }
    }
    return result;
}

template <int ndim> scalar_type FaceSet<ndim>::totalArea() const {
    scalar_type result = 0.0;
    for (index_type i=0; i<_faces.size(); ++i) {
        result += _faces[i]->_area;
    }
    return result;
}

template <int ndim> std::string FaceSet<ndim>::infoString(const bool printAll) const {
    std::ostringstream ss;
    ss << "FaceSet info:" << std::endl;
    ss << "\tgeom = " << geometryString(_geom) << std::endl;
    ss << "\tnMax = " << _nMax << std::endl;
    ss << "\tsize = " << _faces.size() << std::endl;
    ss << "\tnActive = " << _nActive << std::endl;
    ss << "\ttotalArea = " << totalArea() << std::endl;
    ss << "\tmaxFaceArea = " << maxArea() << std::endl;
    ss << "\tminFaceArea = " << minArea() << std::endl;
    if (printAll) {
        for (index_type i=0; i<_faces.size(); ++i) {
            ss << _faces[i]->infoString();
        }
    }
    return ss.str();
}





#ifdef HAVE_VTK
template<int ndim> vtkSmartPointer<vtkCellArray> FaceSet<ndim>::toVtkCellArray() const {
    vtkSmartPointer<vtkCellArray> result = vtkSmartPointer<vtkCellArray>::New();
    const index_type ptsPerFace = _faces[0]->nVerts();
    for(index_type i=0; i<_faces.size(); ++i) {
        if (_faces[i]->isLeaf()) {
            const std::vector<index_type> verts = _faces[i]->vertices();
            result->InsertNextCell(ptsPerFace);
            for(int j=0; j<ptsPerFace; ++j) result->InsertCellPoint(verts[j]);
        }
    }
    return result;
}
#endif

template class FaceSet<2>;
template class FaceSet<3>;
}
}
