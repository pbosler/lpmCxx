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

template <int ndim> class PolyMesh2d;

template <int ndim> class EdgeSet {
    public:
        EdgeSet() : _factory(nullptr), _geom(PLANAR_GEOMETRY), _nMax(0), _nActive(0) {}
    
        EdgeSet(const std::shared_ptr<EdgeFactory<ndim>> factory, const GeometryType geom, const index_type nMax) : 
            _factory(factory), _geom(geom), _nMax(nMax), _nActive(0) {
            _edges.reserve(nMax);
            }

        virtual ~EdgeSet() {}

        inline index_type nMax() const {return _nMax;}
        inline index_type n() const {return _edges.size();}
        inline index_type nActive() const {return _nActive;}
        inline index_type nLeaves() const {return _nActive;}
        inline index_type nDivided() const {return _edges.size() - _nActive;}

        scalar_type maxEucLength(const ParticleSet<ndim>& particles) const;
        scalar_type maxSphLength(const ParticleSet<ndim>& particles,
            const scalar_type radius=1.0) const;
        scalar_type minEucLength(const ParticleSet<ndim>& particles) const;
        scalar_type minSphLength(const ParticleSet<ndim>& particles,
            const scalar_type radius=1.0) const;
            
        inline bool isDivided(const index_type ind) const {return _edges[ind]->isDivided();}
        inline std::array<index_type, 2> kids(const index_type ind) const {return _edges[ind]->kids();}
        inline void setLeftFace(const index_type eind, const index_type find) {_edges[eind]->setLeft(find);}
        inline void setRightFace(const index_type eind, const index_type find) {_edges[eind]->setRight(find);}
        inline index_type orig(const index_type ind) const {return _edges[ind]->orig();}
        inline index_type dest(const index_type ind) const {return _edges[ind]->dest();}

        std::string infoString(const bool printAll = false) const;
        
        inline void enrich(const index_type ind, ParticleSet<ndim> particles) {_edges[ind]->enrich(particles);}
        
        void initFromParticleSetFile(const std::string& fname);
        
        inline bool positiveOrientation(const index_type edgeInd, const index_type faceInd) const {
            return (_edges[edgeInd]->left() == faceInd);
        }
        
        inline scalar_type length(const index_type ind, const ParticleSet<ndim>& particles) const {
            return (_geom == SPHERICAL_SURFACE_GEOMETRY ? _edges[ind]->sphLength(particles) : _edges[ind]->eucLength(particles));}

        inline bool onBoundary(const index_type ind) const {return _edges[ind]->onBoundary();}

        inline Edge<ndim>* getPtr(const index_type ind) const {return _edges[ind].get();}

        void insert(const index_type origID, const index_type destID, const index_type leftID, const index_type rightID, 
            const std::array<index_type, 2>& interiorIDs = std::array<index_type,2>());
        void divide(const index_type ind, ParticleSet<ndim>& particles, const scalar_type radius=1.0);

#ifdef HAVE_VTK
		vtkSmartPointer<vtkCellArray> toVtkCellArray() const;
#endif

		friend class PolyMesh2d<ndim>;

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
