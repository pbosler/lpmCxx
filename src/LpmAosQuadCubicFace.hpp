#ifndef LPM_QUAD_CUBIC_FACE_HPP
#define LPM_QUAD_CUBIC_FACE_HPP

#include "LpmAosFace.hpp"
#include "LpmGll.hpp"

namespace Lpm {
namespace Aos {

template <int ndim> class QuadCubicFace : public Face<ndim> {
    public: 
        QuadCubicFace(const ind_vec& intInds, const ind_vec& vertInds, const ind_vec& edgeInds, 
            const index_type pt, const scalar_type ar=0.0) : Face<ndim>(intInds, vertInds, edgeInds, pt, ar) {}
        
        KidFaceArrays<ndim> divide(ParticleSet<ndim>& particles, EdgeSet<ndim>& edges,
        	const index_type myIndex, const index_type faceInsertPoint, 
        	const scalar_type radius=1.0, const GeometryType geom=PLANAR_GEOMETRY);

		std::vector<Vec<ndim>> getCorners(const ParticleSet<ndim>& particles) const;
		std::vector<Vec<ndim>> getLagCorners(const ParticleSet<ndim>& particles) const;
		        	
        scalar_type computeAreaFromCorners(const ParticleSet<ndim>& particles, const GeometryType geom, 
        	const scalar_type radius) const;
        
        void setArea(const scalar_type ar, ParticleSet<ndim>& particles) override;
        void setArea(ParticleSet<ndim>& particles, const GeometryType geom, const scalar_type radius=1.0) override;
    
    protected:
        static CubicGLL<ndim> gll;
        std::vector<Vec<ndim>> makeInteriors(const std::vector<Vec<ndim>>& corners, const GeometryType geom, const scalar_type radius=1.0) const;
};

template <int ndim> CubicGLL<ndim> QuadCubicFace<ndim>::gll = CubicGLL<ndim>();

}
}
#endif