#ifndef _LPM_AOS_EDGE_FACTORY_HPP
#define _LPM_AOS_EDGE_FACTORY_HPP

#include <memory>
#include <vector>
#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmAosEdge.hpp"

namespace Lpm {
namespace Aos {

template <int ndim> class EdgeFactory {
    public:
        virtual std::unique_ptr<Edge<ndim>> createEdge(const index_type origID, const index_type destID,
            const index_type leftID, const index_type rightID, 
            const std::vector<index_type>& mipts=std::vector<index_type>()) const = 0;
    protected:
};

template <int ndim> class LinearEdgeFactory : public EdgeFactory<ndim> {
    public:
        std::unique_ptr<Edge<ndim>> createEdge(const index_type origID, const index_type destID,
            const index_type leftID, const index_type rightID, 
            const std::vector<index_type>& mipts=std::vector<index_type>()) const 
        {
            return std::unique_ptr<Edge<ndim>>(new Edge<ndim>(origID, destID, leftID, rightID));
        }
    protected:
};

template <int ndim> class QuadraticEdgeFactory : public EdgeFactory<ndim> {
    public:
        std::unique_ptr<Edge<ndim>> createEdge(const index_type origID, const index_type destID,
            const index_type leftID, const index_type rightID, 
            const std::vector<index_type>& midpts) const {
            return std::unique_ptr<Edge<ndim>>(new QuadraticEdge<ndim>(origID, destID, leftID, rightID, midpts));
        }
    protected:
};

template <int ndim> class CubicEdgeFactory : public EdgeFactory<ndim> {
    public:
        std::unique_ptr<Edge<ndim>> createEdge(const index_type origID, const index_type destID,
            const index_type leftID, const index_type rightID, 
            const std::vector<index_type>& midpts) const {
            return std::unique_ptr<Edge<ndim>>(new CubicEdge<ndim>(origID, destID, leftID, rightID, midpts));
        }
    protected:
};

}
}
#endif
