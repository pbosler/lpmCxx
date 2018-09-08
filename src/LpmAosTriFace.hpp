#ifndef LPM_AOS_TRI_FACE_HPP
#define LPM_AOS_TRI_FACE_HPP

#include "LpmAosFace.hpp"

namespace Lpm {
namespace Aos {

template <int ndim> class TriFace : public Face<ndim> {
    public:
        TriFace(const ind_vec& intrs, const ind_vec& verts, const ind_vec& edges, 
            const index_type prnt, const scalar_type ar = 0.0) : 
                Face<ndim>(intrs, verts, edges, prnt, ar) {}
       
        KidFaceArrays<ndim> divide(ParticleSet<ndim>& particles, EdgeSet<ndim>& edges,
        	const index_type myIndex, const index_type faceN, 
        	const scalar_type radius=1.0, const GeometryType geom=PLANAR_GEOMETRY);
        
        std::vector<Vec<ndim>> getCorners(const ParticleSet<ndim>& particles) const;
        
        scalar_type computeAreaFromCorners(const ParticleSet<ndim>& particles, const GeometryType geom, 
        	const scalar_type radius) const;
        
        std::vector<Vec<ndim>> getLagCorners(const ParticleSet<ndim>& particles) const;
};

}
}
#endif 