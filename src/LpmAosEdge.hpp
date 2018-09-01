#ifndef _LPM_AOS_EDGES_HPP
#define _LPM_AOS_EDGES_HPP

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmAosParticle.hpp"
#include "LpmAosParticleSet.hpp"
#include <array>
#include <memory>

namespace Lpm {

namespace Aos {

template <int ndim> class EdgeSet;

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
        
        virtual void enrich(ParticleSet<ndim>& particles) {};
		
		virtual inline index_type ptsPerEdge() const {return 2;}

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
        inline void setKids(const index_type k0, const index_type k1) {
            _kids[0] = k0;
            _kids[1] = k1;
        }
        inline index_type parent() const {return _parent;}
        inline void setParent(const index_type ind) {_parent = ind;}

        Vec<ndim> origCrd(const ParticleSet<ndim>& particles) const {
            return particles.getPtr(_orig)->physCrd();
        }
        Vec<ndim> destCrd(const ParticleSet<ndim>& particles) const {
            return particles.getPtr(_dest)->physCrd();
        }

        Vec<ndim> midpoint(const ParticleSet<ndim>& particles) const;
        Vec<ndim> lagMidpoint(const ParticleSet<ndim>& particles) const;
        Vec<ndim> sphMidpoint(const ParticleSet<ndim>& particles, const scalar_type radius=1.0) const;
        Vec<ndim> sphLagMidpoint(const ParticleSet<ndim>& particles, const scalar_type radius=1.0) const;
        Vec<ndim> edgeVector(const ParticleSet<ndim>& particles) const;
        scalar_type eucLength(const ParticleSet<ndim>& particles) const;
        scalar_type sphLength(const ParticleSet<ndim>& particles, const scalar_type radius=1.0) const;

        inline bool onBoundary() const {return _left < 0 || _right < 0;}

        virtual std::string infoString() const;

        virtual void setMidpt(const index_type ind) {}
        virtual void setMidpts(const index_type ind0, const index_type ind1) {}
        virtual index_type midpt() const {return -1;}
        virtual std::array<index_type,2> midpts() const {
            std::array<index_type,2> result;
            result[0] = -1;
            result[1] = -1;
            return result;
        }

		friend class EdgeSet<ndim>;
    protected:
        index_type _orig;
        index_type _dest;
        index_type _left;
        index_type _right;
        std::array<index_type, 2> _kids;
        index_type _parent;
        std::array<index_type, 2> _midpts;
};

template <int ndim> class QuadraticEdge : public Edge<ndim> {
    public:
        QuadraticEdge(const index_type origID, const index_type destID, const index_type leftID, const index_type rightID,
            const std::array<index_type,2>& midpts) : Edge<ndim>(origID, destID, leftID, rightID, midpts) {}
    
        inline index_type midpt() const override {return this->_midpts[0];}
        inline void setMidpt(const index_type ind) override {this->_midpts[0] = ind;}
        
        inline index_type ptsPerEdge() const override {return 3;}
};

template <int ndim> class CubicEdge : public Edge<ndim> {
    public:
        CubicEdge(const index_type origID, const index_type destID, const index_type leftID, const index_type rightID,
            const std::array<index_type,2>& midpts) : Edge<ndim>(origID, destID, leftID, rightID, midpts) {}
        inline std::array<index_type, 2> midpts() const override {return this->_midpts;}
        inline void setMidpts(const index_type ind0, const index_type ind1) override {
            this->_midpts[0] = ind0;
            this->_midpts[1] = ind1;
        }
        inline index_type ptsPerEdge() const override {return 4;}
};
}
}
#endif
