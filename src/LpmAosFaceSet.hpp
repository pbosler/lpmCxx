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
        FaceSet(const std::shared_ptr<FaceFactory<ndim>> fac, const GeometryType geom) :
             _factory(factory), _geom(geom) {}
    
        void insert(const ind_vec& intrs, const ind_vec& verts, const ind_vec& edges, 
            const index_type pt, const scalar_type ar=0.0)
    
    protected:
        GeometryType _geom;
        index_type _nMax;
        index_type _nActive;
        std::shared_ptr<FaceFactory<ndim>> _factory;
        std::vector<std::unique_ptr<Face<ndim>>> _faces;
};

}

#endif
