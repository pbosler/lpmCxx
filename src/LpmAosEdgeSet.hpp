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

namespace Lpm {

template <int ndim> class EdgeSet {
    public:
        EdgeSet(const std::shared_ptr<EdgeFactory<ndim>> factory, const GeometryType geom) : _factory(factory), _geom(geom) {}

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

        std::string infoString() const;
        
        void initFromParticleSetFile(const std::string& fname);
        
        inline scalar_type length(const index_type ind, const ParticleSet<ndim>& particles) const {
            return (_geom == SPHERICAL_SURFACE_GEOMETRY ? _edges[ind]->sphLength(particles) : _edges[ind]->eucLength(particles));}

        inline bool onBoundary(const index_type ind) const {return _edges[ind]->onBoundary();}

        inline Edge<ndim>* getPtr(const index_type ind) const {return _edges[ind].get();}

        void insert(const index_type origID, const index_type destID, const index_type leftID, const index_type rightID, 
            const std::array<index_type, 2>& interiorIDs = std::array<index_type,2>());
        void divide(const index_type ind, ParticleSet<ndim>& particles, const scalar_type radius=1.0);

    protected:
        GeometryType _geom;
        index_type _nMax;
        index_type _nActive;
        std::shared_ptr<EdgeFactory<ndim>> _factory;
        std::vector<std::unique_ptr<Edge<ndim>>> _edges;
};

}
#endif
