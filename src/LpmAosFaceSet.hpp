#ifndef _LPM_AOS_FACE_SET_HPP
#define _LPM_AOS_FACE_SET_HPP

#include <vector>
#include <memory>
#include <string>
#include <exception>
#include <algorithm>
#include "LpmConfig.h"
#include "LpmAosTypes.hpp"
#include "LpmAosParticleSet.hpp"
#include "LpmAosEdgeSet.hpp"
#include "LpmAosFace.hpp"
#include "LpmAosFaceFactory.hpp"

#ifdef HAVE_VTK
#include "vtkSmartPointer.h"
#include "vtkCellArray.h"
#endif

namespace Lpm {
namespace Aos {

template <int ndim> class PolyMesh2d;

template <int ndim> class FaceSet {
    public:
        FaceSet() : _factory(nullptr), _nMax(0), _nActive(0), _geom(PLANAR_GEOMETRY), _radius(1.0) {}
    
        FaceSet(const std::shared_ptr<FaceFactory<ndim>> fac, const index_type nMax, const GeometryType geom, 
            const scalar_type r=1.0) : _factory(fac), _nMax(nMax), _nActive(0), _geom(geom), _radius(r) 
        {
            _faces.reserve(nMax);
        }
    
        void insert(const ind_vec& intrs, const ind_vec& verts, const ind_vec& edges, 
            const index_type parent, const scalar_type ar=0.0);
            
        void divide(const index_type ind, ParticleSet<ndim>& particles, EdgeSet<ndim>& edges);
        
        inline void enrich(const index_type ind, ParticleSet<ndim>& particles, EdgeSet<ndim>& edges) {
            _faces[ind]->enrich(particles, edges);
        }
        
        inline void setArea(const ParticleSet<ndim>& particles, const scalar_type radius=1.0) {
            for (index_type i=0; i<_faces.size(); ++i) {
                if (_faces[i]->isLeaf()) {
                    _faces[i]->setArea(_geom, particles, radius);
                }
                else {
                    _faces[i]->setArea(0.0);
                }
            }
        }
        
        inline index_type n() const {return _faces.size();}
        inline index_type nMax() const {return _nMax;}
        inline index_type nActive() const {return _nActive;}
        inline index_type nLeaves() const {return _nActive;}
        inline index_type nDivided() const {return _faces.size() - _nActive;}

        scalar_type minArea() const;
        scalar_type maxArea() const;
        scalar_type maxLeafArea() const;
        scalar_type totalArea() const;
        
        inline scalar_type area(const index_type ind) const {return _faces[ind]->_area;}
        inline bool isLeaf(const index_type ind) const {return _faces[ind]->isLeaf();}
        
        std::string infoString(const bool printAll = false) const;

#ifdef HAVE_VTK
		vtkSmartPointer<vtkCellArray> toVtkCellArray() const;
#endif
    
    	friend class PolyMesh2d<ndim>;
    protected:
        scalar_type _radius;
        GeometryType _geom;
        index_type _nMax;
        index_type _nActive;
        std::shared_ptr<FaceFactory<ndim>> _factory;
        std::vector<std::unique_ptr<Face<ndim>>> _faces;
};

}
}
#endif
