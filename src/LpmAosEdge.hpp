#ifndef _LPM_AOS_EDGES_HPP
#define _LPM_AOS_EDGES_HPP

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmAosTypes.hpp"
#include "LpmAosParticle.hpp"
#include "LpmAosParticleSet.hpp"
#include <array>
#include <vector>
#include <memory>

namespace Lpm {

namespace Aos {

template <int ndim> class EdgeSet;

template <int ndim> struct KidEdgeArrays {
	std::vector<index_type> newOrigs;
	std::vector<index_type> newDests;
	std::vector<index_type> newLefts;
	std::vector<index_type> newRights;
	std::vector<std::vector<index_type>> newMids;
	std::vector<std::unique_ptr<Particle<ndim>>> particles;
	
	KidEdgeArrays(): newOrigs(std::vector<index_type>(2,-1)), newDests(std::vector<index_type>(2,-1)), 
		newLefts(std::vector<index_type>(2,-1)), newRights(std::vector<index_type>(2,-1)),
		newMids(std::vector<std::vector<index_type>>(2,std::vector<index_type>(2,-1))) {particles.reserve(4);}
	
	std::string infoString() const;
};

template <int ndim> class Edge {
    public:
        Edge(const index_type origID, const index_type destID, const index_type leftID, const index_type rightID,
            const std::array<index_type,2>& midpts = std::array<index_type,2>()) :
            _orig(origID), _dest(destID), _left(leftID), _right(rightID),
            _parent(-1) {
            _kids[0] = -1;
            _kids[1] = -1;
            _midpts[0] = midpts[0];
            _midpts[1] = midpts[1];
        }

        virtual ~Edge() {}
        
        virtual KidEdgeArrays<ndim> divide(const ParticleSet<ndim>& particles, const scalar_type radius=1.0, const GeometryType geom=PLANAR_GEOMETRY) const;
        
        /// Enrich a low-order edge to a higher order edge
        virtual void enrich(ParticleSet<ndim>& particles) {};
		
		virtual inline index_type ptsPerEdge() const {return 2;}

		/// True if Edge has been divided
        inline bool hasKids() const {return _kids[0] > 0;}
        /// True if Edge has not been divided
        inline bool isLeaf() const {return !hasKids();}
        /// True if Edge has been divided
        inline bool isDivided() const {return hasKids();}

		/// Returns the index of the Edge's origin particle in a corresponding ParticleSet<ndim> instance
        inline index_type orig() const {return _orig;}
        /// Returns the index of the Edge's destination particle in a corresponding ParticleSet<ndim> instance
        inline index_type dest() const {return _dest;}
        /// Returns the index of the Edge's left Face in a corresponding FaceSet<ndim> instance
        inline index_type left() const {return _left;}
        /// Returns the index of the Edge's right Face in a corresponding FaceSet<ndim> instance
        inline index_type right() const {return _right;}
        /// Returns the index of the Edge's interior particle in a corresponding ParticleSet<ndim> instance
        inline index_type midpt() const {return _midpts[0];}
        /// Returns the index of the Edge's interior particle in a corresponding ParticleSet<ndim> instance
        inline index_type midpt0() const {return _midpts[0];}
        /// Returns the index of the Edge's interior particle in a corresponding ParticleSet<ndim> instance
        inline index_type midpt1() const {return _midpts[1];}
        /// Returns the indices of the Edge's interior particles in a corresponding ParticleSet<ndim> instance
        inline std::array<index_type,2> midpts() const {return _midpts;}
        
        /// Set midpoint indices
        inline void setMidptIdds(const index_type ind0, const index_type ind1) {
        	_midpts[0] = ind0;
        	_midpts[1] = ind1;
        }

		/// Set Left Face index (to a Face in a corresponding FaceSet)
        inline void setLeft(const index_type ind) {_left = ind;}
        /// Set Right Face index (to a Face in a corresponding FaceSet)
        inline void setRight(const index_type ind){_right= ind;}
		/// Return the indices of this edge's children (in a corresponding EdgeSet)
        inline std::array<index_type, 2> kids() const {return _kids;}
        /// Set the indices of this Edge's children
        inline void setKids(const index_type k0, const index_type k1) {
            _kids[0] = k0;
            _kids[1] = k1;
        }
        /// Return the parent edge index (in a corresponding EdgeSet)
        inline index_type parent() const {return _parent;}
        /// Set the parent edge index
        inline void setParent(const index_type ind) {_parent = ind;}

		/// Return the physical coordinate of the Edge origin
        inline Vec<ndim> origCrd(const ParticleSet<ndim>& particles) const {
            return particles.physCrd(_orig);
        }
        
        /// Return the physical cooridnate of the Edge destination
        inline Vec<ndim> destCrd(const ParticleSet<ndim>& particles) const {
            return particles.physCrd(_dest);
        }
        
        /// Return the Lagrangian coordinate of the Edge origin
        inline Vec<ndim> origLagCrd(const ParticleSet<ndim>& particles) const {
        	return particles.lagCrd(_orig);
        }
        
        /// Return the Lagrangian coordiante of the Edge destination
        inline Vec<ndim> destLagCrd(const ParticleSet<ndim>& particles) const {
        	return particles.lagCrd(_dest);
        }
        
        /// Return the physical coordinate of the Edge's first interior point.
        inline Vec<ndim> mid0Crd(const ParticleSet<ndim>& particles) const {
        	return particles.physCrd(_midpts[0]);
        }
        /// Return the Lagrangian coordinate of the Edge 's first interior point.
        inline Vec<ndim> mid0LagCrd(const ParticleSet<ndim>& particles) const {
        	return particles.lagCrd(_midpts[0]);
        }
        /// Return the physical coordinate of the Edge's second interior point.
		inline Vec<ndim> mid1Crd(const ParticleSet<ndim>& particles) const {
        	return particles.physCrd(_midpts[1]);
        }
        /// Return the Lagrangian coordinate of the Edge's second interior point.
        inline Vec<ndim> mid1LagCrd(const ParticleSet<ndim>& particles) const {
        	return particles.lagCrd(_midpts[1]);
        }

		/// Compute the edge midpoint in physical space
        Vec<ndim> midpoint(const ParticleSet<ndim>& particles) const;
        /// Compute the edge midpoint in Lagrangian space
        Vec<ndim> lagMidpoint(const ParticleSet<ndim>& particles) const;
        /// Compute the edge midpoint in physical space
        Vec<ndim> sphMidpoint(const ParticleSet<ndim>& particles, const scalar_type radius=1.0) const;
        /// Compute the edge midpoint in Lagrangian space
        Vec<ndim> sphLagMidpoint(const ParticleSet<ndim>& particles, const scalar_type radius=1.0) const;
        /// Return the Cartesian vector pointing from origin to destination
        Vec<ndim> edgeVector(const ParticleSet<ndim>& particles) const;
        /// Return the Euclidean length of the Edge
        scalar_type eucLength(const ParticleSet<ndim>& particles) const;
        /// Return the great-cricle length of the Edge
        scalar_type sphLength(const ParticleSet<ndim>& particles, const scalar_type radius=1.0) const;

		/// Returns True if the Edge lies on or defines the boundary of the computation domain
        inline bool onBoundary() const {return _left < 0 || _right < 0;}
        
        /// Returns the Edge's unit tangent vector (at the edge midpoint)
        Vec<ndim> unitTangent(const ParticleSet<ndim>& particles) const;
        /// Returns the Eges' unit normal vector (at the edge midpoint)
        Vec<ndim> unitNormal(const ParticleSet<ndim>& particles, const GeometryType geom) const;

		/// Output Edge data to a string (usually then piped to std::cout)
        virtual std::string infoString() const;

		// friend class EdgeSet<ndim>;
    protected:
        index_type _orig;
        index_type _dest;
        index_type _left;
        index_type _right;
        std::array<index_type, 2> _kids;
        index_type _parent;
        std::array<index_type, 2> _midpts;
        
        scalar_type _length;
        Vec<ndim> _normal;
        Vec<ndim> _tangent;
};

template <int ndim> class QuadraticEdge : public Edge<ndim> {
    public:
        QuadraticEdge(const index_type origID, const index_type destID, const index_type leftID, const index_type rightID,
            const std::array<index_type,2>& midpts) : Edge<ndim>(origID, destID, leftID, rightID, midpts) {}

        
        inline index_type ptsPerEdge() const override {return 3;}
        
        KidEdgeArrays<ndim> divide(const ParticleSet<ndim>& particles, const scalar_type radius=1.0, const GeometryType geom=PLANAR_GEOMETRY) const override;
};

template <int ndim> class CubicEdge : public Edge<ndim> {
    public:
        CubicEdge(const index_type origID, const index_type destID, const index_type leftID, const index_type rightID,
            const std::array<index_type,2>& midpts) : Edge<ndim>(origID, destID, leftID, rightID, midpts) {}
            
        inline index_type ptsPerEdge() const override {return 4;}
        
        KidEdgeArrays<ndim> divide(const ParticleSet<ndim>& particles, const scalar_type radius=1.0, const GeometryType geom=PLANAR_GEOMETRY) const override;
        
       // kidsArray divide(const std::shared_ptr<ParticleFactory<ndim>> pfac) const override;
};
}
}
#endif
