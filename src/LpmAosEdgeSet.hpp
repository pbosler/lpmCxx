#ifndef _LPM_AOS_EDGE_SET_HPP
#define _LPM_AOS_EDGE_SET_HPP

#include <vector>
#include <memory>
#include <string>
#include <exception>
#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmAosEdge.hpp"
#include "LpmAosEdgeFactory.hpp"
#include "LpmAosParticleSet.hpp"
#ifdef HAVE_VTK
#include "vtkSmartPointer.h"
#include "vtkCellArray.h"
#endif

namespace Lpm {

namespace Aos {

// template <int ndim> class PolyMesh2d;

/// Collect of Edge Pointers arranged in a binary tree.
template <int ndim> class EdgeSet {
    public:
    	/// Constructor. Initializes an empty set.
        EdgeSet() : _factory(nullptr), _geom(PLANAR_GEOMETRY), _nMax(0), _nActive(0) {}
    
    	/// Constructor. Initializes an empty set.
        EdgeSet(const std::shared_ptr<EdgeFactory<ndim>> factory, const GeometryType geom, const index_type nMax) : 
            _factory(factory), _geom(geom), _nMax(nMax), _nActive(0) {
            _edges.reserve(nMax);
            }

        virtual ~EdgeSet() {}

		/// Return the maximum number of edges allowed in memory.
        inline index_type nMax() const {return _nMax;}
        /// Returns the current number of edges (leaves and interior nodes of the tree) in memory.
        inline index_type n() const {return _edges.size();}
        /// Returns the number of leaves in the edgeset tree
        inline index_type nActive() const {return _nActive;}
        /// Returns the number of leaves in the edgeset tree
        inline index_type nLeaves() const {return _nActive;}
        /// Returns the number of divided edges (interior tree nodes) in the EdgeSet tree
        inline index_type nDivided() const {return _edges.size() - _nActive;}

		/// Returns the maximum Euclidean distance spanned by leaf edges
        scalar_type maxEucLength(const ParticleSet<ndim>& particles) const;
        /// Returns the maximum great-circle distance spanned by leaf edges
        scalar_type maxSphLength(const ParticleSet<ndim>& particles,
            const scalar_type radius=1.0) const;
        /// Returns the minimum Euclidean distance spanned by leaf edges
        scalar_type minEucLength(const ParticleSet<ndim>& particles) const;
        /// Returns the minimum great-circle distance spanned by leaf edges
        scalar_type minSphLength(const ParticleSet<ndim>& particles,
            const scalar_type radius=1.0) const;
        
        /// Returns true if the edge at index ind has been divided
        inline bool isDivided(const index_type ind) const {return _edges[ind]->isDivided();}
        /// Returns the indices (to edges within this EdgeSet) of an edge's children
        inline std::array<index_type, 2> kids(const index_type ind) const {return _edges[ind]->kids();}
        /// Set the left face of an edge
        inline void setLeftFace(const index_type eind, const index_type find) {_edges[eind]->setLeft(find);}
        /// Set the right face of an edge
        inline void setRightFace(const index_type eind, const index_type find) {_edges[eind]->setRight(find);}
        /// Return the origin index of an edge (the index is associated with a ParticleSet instance)
        inline index_type orig(const index_type ind) const {return _edges[ind]->orig();}
        /// Return the destination index of an edge (the index is associated with a ParticleSet instance)
        inline index_type dest(const index_type ind) const {return _edges[ind]->dest();}

		/// Collect info about an EdgeSet, output to string.
        std::string infoString(const bool printAll = false) const;
        
        /// \todo Convert a low-order edge to a higher order egde
        inline void enrich(const index_type ind, ParticleSet<ndim> particles) {_edges[ind]->enrich(particles);}
        
        /// \todo Initialize an EdgeSet instance from a data file.
        void initFromParticleSetFile(const std::string& fname);
        
        /// Returns true if an edge is positively oriented relative to a face.
        /**
        	\warning this does not check to make sure that the edge is connected to the requested face.
        */
        inline bool positiveOrientation(const index_type edgeInd, const index_type faceInd) const {
            return (_edges[edgeInd]->left() == faceInd);
        }
        
        /// Returns the length of an edge
        inline scalar_type length(const index_type ind, const ParticleSet<ndim>& particles) const {
            return (_geom == SPHERICAL_SURFACE_GEOMETRY ? _edges[ind]->sphLength(particles) : _edges[ind]->eucLength(particles));}

		/// Returns true if the edge is on the boundary of the computational domain.
        inline bool onBoundary(const index_type ind) const {return _edges[ind]->onBoundary();}
	
		/// Returns the raw pointer to an Edge within an EdgeSet instance.
        inline Edge<ndim>* getPtr(const index_type ind) const {return _edges[ind].get();}

		/// Inserts a new edge to the EdgeSet
		/**
			New edges are inserted at the end of the linear edges array in memory; they maintain pointers (indices) to 
			their parents and kids to define the tree.
		*/
        void insert(const index_type origID, const index_type destID, const index_type leftID, const index_type rightID, 
            const std::vector<index_type>& interiorIDs = std::vector<index_type>());
        
        /// Divides a particular edge.  Its children are inserted (along with any new particles) to the EdgeSet (and ParticleSet).
        void divide(const index_type ind, ParticleSet<ndim>& particles, const scalar_type radius=1.0);

#ifdef HAVE_VTK
		/// Convert an EdgeSet to a VTK cell array for use with VTK.
		vtkSmartPointer<vtkCellArray> toVtkCellArray() const;
#endif

/* 		friend class PolyMesh2d<ndim>; */

    protected:
        GeometryType _geom;
        index_type _nMax;
        index_type _nActive;
        std::shared_ptr<EdgeFactory<ndim>> _factory;
        std::vector<std::unique_ptr<Edge<ndim>>> _edges;
};

}
}
#endif
