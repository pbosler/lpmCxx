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

class EdgeSet {
    public:
        EdgeSet(const std::shared_ptr<EdgeFactory> factory, const GeometryType geom) : _factory(factory), _geom(geom) {}

        virtual ~EdgeSet() {}

        inline index_type nMax() const {return _nMax;}
        inline index_type n() const {return _edges.size();}
        inline index_type nActive() const {return _nActive;}
        inline index_type nLeaves() const {return _nActive;}
        inline index_type nDivided() const {return _edges.size() - _nActive;}

        template <int ndim> scalar_type maxEucLength(const ParticleSet<ndim>& particles) const;
        template <int ndim> scalar_type maxSphLength(const ParticleSet<ndim>& particles,
            const scalar_type radius=1.0) const;
        template <int ndim> scalar_type minEucLength(const ParticleSet<ndim>& particles) const;
        template <int ndim> scalar_type minSphLength(const ParticleSet<ndim>& particles,
            const scalar_type radius=1.0) const;

        virtual std::string infoString() const;

        inline Edge* getPtr(const index_type ind) const {return _edges[i].get();}

        void insert(const index_type origID, const index_type destID, const index_type leftID, const index_type rightID);
        template <int ndim> void divide(const index_type ind, ParticleSet<ndim>& particles, const scalar_type radius=1.0);

    protected:
        GeometryType _geom;
        index_type _nMax;
        index_type _nActive;
        std::shared_ptr<EdgeFactory> _factory;
        std::vector<std::unique_ptr<Edge>> _edges;
};

}
#endif
