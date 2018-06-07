#ifndef _LPM_AOS_FACE_SET_HPP
#define _LPM_AOS_FACE_SET_HPP

#include <vector>
#include <memory>
#include <string>
#include <exception>
#include "LpmConfig.h"
#include "LpmAosTypes.hpp"
#include "LpmAosParticleSet.hpp"
#include "LpmAosEdgeSet.hpp"
#include "LpmAosFace.hpp"
#include "LpmAosFaceFactory.hpp"

namespace Lpm {

template <int ndim> class FaceSet {
    public:
        FaceSet(const std::shared_ptr<FaceFactory<ndim>> fac, const index_type nMax, const GeometryType geom) :
             _factory(fac), _nMax(nMax), _geom(geom) {}
    
        void insert(const ind_vec& intrs, const ind_vec& verts, const ind_vec& edges, 
            const index_type pt, const scalar_type ar=0.0);
            
        void divide(const index_type ind, ParticleSet<ndim>& particles, EdgeSet<ndim>& edges);
        
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
        
        std::string infoString() const;
    
    protected:
        GeometryType _geom;
        index_type _nMax;
        index_type _nActive;
        std::shared_ptr<FaceFactory<ndim>> _factory;
        std::vector<std::unique_ptr<Face<ndim>>> _faces;
};

}

#endif
