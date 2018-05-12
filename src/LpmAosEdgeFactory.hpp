#ifndef _LPM_AOS_EDGE_FACTORY_HPP
#define _LPM_AOS_EDGE_FACTORY_HPP

#include <memory>
#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmAosEdge.hpp"

namespace Lpm {

class EdgeFactory {
    public:
        virtual std::unique_ptr<Edge> createEdge(const index_type origID, const index_type destID,
            const index_type leftID, const index_type rightID) const = 0;
    protected:
};

class LinearEdgeFactory : public EdgeFactory {
    public:
        std::unique_ptr<Edge> createEgde(const index_type origID, const index_type destID,
            const index_type leftID, const index_type rightID) const {
            return std::unique_ptr<Edge>(new Edge(const index_type origID, const index_type destID,
            const index_type leftID, const index_type rightID));
        }
    protected:
};

class QuadraticEdgeFactory : public EdgeFactory {
    public:
        std::unique_ptr<Edge> createEdge(const index_type origID, const index_type destID,
            const index_type leftID, const index_type rightID) const {
            return std::unique_ptr<Edge>(new QuadraticEdge(const index_type origID, const index_type destID,
            const index_type leftID, const index_type rightID));
        }
    protected:
};

class CubicEdgeFactory : public EdgeFactory {
    public:
        std::unique_ptr<Edge> createEdge(const index_type origID, const index_type destID,
            const index_type leftID, const index_type rightID) const {
            return std::unique_ptr<Edge>(new CubicEdge(const index_type origID, const index_type destID,
            const index_type leftID, const index_type rightID));
        }
    protected:
};

}

#endif
