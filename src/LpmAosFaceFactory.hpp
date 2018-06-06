#ifndef _LPM_AOS_FACE_FACTORY_HPP
#define _LPM_AOS_FACE_FACTORY_HPP

#include <memory>
#include <vector>
#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmAosFace.hpp"

namespace Lpm {

typedef std::vector<index_type> ind_vec;

template <int ndim> class FaceFactory {
    
    public:
        virtual std::unique_ptr<Face<ndim>> createFace(const ind_vec& intrs, const ind_vec& verts, 
            const ind_vec& edges, const index_type pt, const scalar_type ar=0.0) const = 0;
};

template <int ndim> class TriFaceFactory : public FaceFactory<ndim> {
    public: 
        std::unique_ptr<Face<ndim>> createFace(const ind_vec& intrs, const ind_vec& verts, 
            const ind_vec& edges, const index_type pt, const scalar_type ar=0.0) const {
            return std::unique_ptr<Face<ndim>>(new TriFace<ndim>(intrs, verts, edges, pt, ar));
        }
};

template <int ndim> class QuadFaceFactory : public FaceFactory<ndim> {
    public:
        std::unique_ptr<Face<ndim>> createFace(const ind_vec& intrs, const ind_vec& verts, 
            const ind_vec& edges, const index_type pt, const scalar_type ar=0.0) const {
            return std::unique_ptr<Face<ndim>>(new QuadFace<ndim>(intrs, verts, edges, pt, ar));
        }
};

template <int ndim> class QuadCubicFaceFactory : public FaceFactory<ndim> {
    public:
        std::unique_ptr<Face<ndim>> createFace(const ind_vec& intrs, const ind_vec& verts, 
            const ind_vec& edges, const index_type pt, const scalar_type ar=0.0) const {
            return std::unique_ptr<Face<ndim>>(new QuadCubicFace<ndim>(intrs, verts, edges, pt, ar));
        }
};

}

#endif