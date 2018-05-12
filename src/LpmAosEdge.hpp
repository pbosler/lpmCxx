#ifndef _LPM_AOS_EDGES_HPP
#define _LPM_AOS_EDGES_HPP

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmAosParticle.hpp"
#include "LpmAosParticleSet.hpp"
#include <array>
#include <memory>

namespace Lpm {

class Edge {
    public:
        Edge(const index_type origID, const index_type destID, const index_type leftID, const index_type rightID) :
            _orig(origID), _dest(destID), _left(leftID), _right(rightID),
            _parent(-1) {
            _kids[0] = -1;
            _kids[1] = -1;
        }

        virtual ~Edge() {}

        inline bool hasKids() const {return _kids[0] > 0;}
        inline bool isLeaf() const {return !hasKids();}
        inline bool isDivided() const {return hasKids();}

        inline index_type orig() const {return _orig;}
        inline index_type dest() const {return _dest;}
        inline index_type left() const {return _left;}
        inline index_type right() const {return _right;}

        inline void setLeft(const index_type ind) {_left = ind;}
        inline void setRight(const index_type ind){_right= ind;}

        inline std::array<index_type, 2> kids() const {return _kids;}

        template <int ndim> Vec<ndim> midpoint(const ParticleSet<ndim>& particles) const;
        template <int ndim> Vec<ndim> lagMidpoint(const ParticleSet<ndim>& particles) const;
        template <int ndim> Vec<ndim> edgeVector(const ParticleSet<ndim>& particles) const;
        template <int ndim> scalar_type eucLength(const ParticleSet<ndim>& particles) const;
        template <int ndim> scalar_type sphLength(const ParticleSet<ndim>& particles, const scalar_type radius=1.0) const;

        inline bool onBoundary() const {return _left < 0 || _right < 0;}

        virtual std::string infoString() const;

    protected:
        index_type _orig;
        index_type _dest;
        index_type _left;
        index_type _right;
        std::array<index_type, 2> _kids;
        index_type _parent;
};

class QuadraticEdge : public Edge {
    public:
        QuadraticEdge(const index_type origID, const index_type destID, const index_type leftID, const index_type rightID,
            const index_type midID) : Edge(origID, destID, leftID, rightID), _midpt(midID) {}

        inline index_type midpt() const {return _midpt;}
    protected:
        index_type _midpt;
};

class CubicEdge : public Edge {
    public:
        CubicEdge(const index_type origID, const index_type destID, const index_type leftID, const index_type rightID,
            const index_type midID0, const index_type midID1) : Edge(origID, destID, leftID, rightID) {
            _midpts[0] = midID0;
            _midpts[1] = midID1;
        }

        inline std::array<index_type, 2> midpts() const {return _midpts;}
    protected:
        std::array<index_type, 2> _midpts;
};

}
#endif
