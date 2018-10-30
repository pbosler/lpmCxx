#ifndef _LPM_AOS_FACE_HPP
#define _LPM_AOS_FACE_HPP

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmAosParticleSet.hpp"
#include "LpmAosEdgeSet.hpp"
#include <array>
#include <memory>
#include <vector>
#include <map>

namespace Lpm {
namespace Aos {

typedef std::vector<index_type> ind_vec;
template <int ndim> class FaceSet;

template <int ndim> struct KidFaceArrays {
	std::vector<ind_vec> newFaceVerts;
	std::vector<ind_vec> newFaceEdges;
	std::vector<ind_vec> newFaceInteriors;
	std::vector<scalar_type> kidsFaceArea;
	
	KidFaceArrays() : newFaceVerts(4,ind_vec(12, -1)), newFaceEdges(4,ind_vec(4,-1)), 
	newFaceInteriors(4,ind_vec(4,-1)), kidsFaceArea(std::vector<scalar_type>(4,0.0)) {} 
	
	KidFaceArrays(const int nverts, const int nedges, const int ninteriors) : newFaceVerts(4, ind_vec(nverts,-1)), 
		newFaceEdges(4, ind_vec(nedges,-1)), newFaceInteriors(4, ind_vec(ninteriors,-1)), 
		kidsFaceArea(std::vector<scalar_type>(4,0.0)){}

	std::string infoString() const;
};	

template <int ndim> class Face {
    typedef std::vector<scalar_type> vfield_type;

    public:
        virtual ~Face() {}
        
        virtual void enrich(ParticleSet<ndim>& particles, EdgeSet<ndim>& edges) {};
        
        virtual KidFaceArrays<ndim> divide(ParticleSet<ndim>& particles, EdgeSet<ndim>& edges,
        	const index_type myIndex, const index_type faceInsertPoint, 
        	const scalar_type radius=1.0, const GeometryType geom=PLANAR_GEOMETRY) = 0;
        	
        virtual std::vector<Vec<ndim>> getCorners(const ParticleSet<ndim>& particles) const = 0;
        virtual std::vector<Vec<ndim>> getLagCorners(const ParticleSet<ndim>& particles) const = 0;

        virtual scalar_type computeAreaFromCorners(const ParticleSet<ndim>& particles, const GeometryType geom, 
        	const scalar_type radius) const = 0;
        
        inline bool hasKids() const {return _kids[0] >= 0;}
        inline bool isDivided() const {return hasKids();}
        inline bool isLeaf() const {return !hasKids();}
        inline bool isActive() const {return isLeaf();}
        
        Vec<ndim> physBarycenter(const ParticleSet<ndim>& particles) const;
        Vec<ndim> lagBarycenter(const ParticleSet<ndim>& particles) const;
        
        Vec<3> physSphBarycenter(const ParticleSet<ndim>& particles, const scalar_type radius = 1.0) const;
        Vec<3> lagSphBarycenter(const ParticleSet<ndim>& particles, const scalar_type radius = 1.0) const;
        
        std::string infoString() const;
        
        inline ind_vec interiors() const {return _interiorInds;}
        inline ind_vec vertices() const {return _vertInds;}
        inline ind_vec edges() const {return _edgeInds;}

		std::vector<Vec<ndim>> vertPhysCrds(const ParticleSet<ndim>& particles) const;
	
        inline index_type nIntrs() const {return _interiorInds.size();}
        inline index_type nVerts() const {return _vertInds.size();}
        
        inline void setInteriorInds(const ind_vec& inds) {_interiorInds = inds;}
        inline void setVerts(const ind_vec& inds) {_vertInds = inds;}
        inline void setEdges(const ind_vec& inds) {_edgeInds = inds;}
        
        inline std::array<index_type, 4> kids() const {return _kids;}
        inline void setKid(const int child, const index_type ind) {_kids[child] = ind;}
        inline void setKids(const std::array<index_type,4>& inds) {_kids = inds;}
        inline void setKids(const std::vector<index_type>& inds) { for (short i=0; i<4; ++i) _kids[i] = inds[i];}
        
        inline index_type parent() const {return _parent;}
        inline void setParent(const index_type ind) {_parent = ind;}

        inline scalar_type area() const {return _area;}
        inline void setArea(const scalar_type ar) {_area = ar;}
        virtual void setArea(const scalar_type ar, ParticleSet<ndim>& particles);
        virtual void setArea(ParticleSet<ndim>& particles, const GeometryType geom, const scalar_type radius);

        scalar_type scalarIntegral(const std::string field_name, const ParticleSet<ndim>& particles) const ;
    
            

    protected:
    	Face(const ind_vec& intInds, const ind_vec& vertInds, const ind_vec& edgeInds, const index_type pt, const scalar_type ar=0.0) : 
                _interiorInds(intInds), _vertInds(vertInds), _edgeInds(edgeInds), _parent(pt), _area(ar) {
             for (int i=0; i<4; ++i) {
                _kids[i] = -1;
            }
        }
    
        Face() {}
    
        ind_vec _interiorInds;
        ind_vec _vertInds;
        ind_vec _edgeInds;
        index_type _parent;
        std::array<index_type, 4> _kids;
        scalar_type _area;
        Vec<ndim> _normal;
        
        std::map<std::string, scalar_type> _sfields;
        std::map<std::string, std::vector<scalar_type>> _vfields;
};
}
}
#endif
