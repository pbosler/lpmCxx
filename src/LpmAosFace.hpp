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

template <int ndim> class Face {
    typedef std::vector<scalar_type> vfield_type;

    public:
        Face(const ind_vec& intInds, const ind_vec& vertInds, const ind_vec& edgeInds, const index_type pt, const scalar_type ar=0.0) : 
                _interiorInds(intInds), _vertInds(vertInds), _edgeInds(edgeInds), _parent(pt), _area(ar) {
             for (int i=0; i<4; ++i) {
                _kids[i] = -1;
            }
        }
        
        virtual ~Face() {}
        
        virtual void enrich(ParticleSet<ndim>& particles, EdgeSet<ndim>& edges) {};
        
        inline bool hasKids() const {return _kids[0] >= 0;}
        inline bool isDivided() const {return hasKids();}
        inline bool isLeaf() const {return !hasKids();}
        
        Vec<ndim> physBarycenter(const ParticleSet<ndim>& particles) const;
        Vec<ndim> lagBarycenter(const ParticleSet<ndim>& particles) const;
        
        Vec<3> physSphBarycenter(const ParticleSet<ndim>& particles, const scalar_type radius = 1.0) const;
        Vec<3> lagSphBarycenter(const ParticleSet<ndim>& particles, const scalar_type radius = 1.0) const;
        
        std::string infoString() const;
        
        inline ind_vec interiors() const {return _interiorInds;}
        inline ind_vec vertices() const {return _vertInds;}
        inline ind_vec edges() const {return _edgeInds;}
        
        inline index_type nVerts() const {return _vertInds.size();}
        
        std::vector<Vec<ndim>> corners() const;
        
        inline void setInteriorInds(const ind_vec& inds) {_interiorInds = inds;}
        inline void setVerts(const ind_vec& inds) {_vertInds = inds;}
        inline void setEdges(const ind_vec& inds) {_edgeInds = inds;}
        
        inline std::array<index_type, 4> kids() const {return _kids;}
        inline void setKid(const int child, const index_type ind) {_kids[child] = ind;}
        inline void setKids(const std::array<index_type,4>& inds) {_kids = inds;}
        
        inline index_type parent() const {return _parent;}
        inline void setParent(const index_type ind) {_parent = ind;}
        
        scalar_type area() const {return _area;}
        inline void setArea(const scalar_type ar) {_area = ar;}
        scalar_type computeArea(const ParticleSet<ndim>& vertexParticles, const ParticleSet<ndim>& faceParticles);
        void setArea(const GeometryType geom, const ParticleSet<ndim>& particles, const scalar_type radius=1.0);    
            
        friend class FaceSet<ndim>;
        
        inline void registerScalarField(const std::string& field_name) {
            this->_sfields.emplace(field_name, 0.0);
        }
        
        inline void registerVectorField(const std::string& field_name) {
            this->_vfields.emplace(field_name, vfield_type(ndim,0.0));
        }
        
        inline std::vector<std::string> scalarFieldNames() const {
            std::vector<std::string> result;
            for (auto& sf : _sfields) {
                result.push_back(sf.first);
            }
            return result;
        }
        
        inline std::vector<std::string> vectorFieldNames() const {
            std::vector<std::string> result;
            for (auto& vf : _vfields) {
                result.push_back(vf.first);
            }
            return result;
        }
        
        inline std::vector<std::string> fieldNames() const {
            std::vector<std::string> result;
            for (auto& sf : _sfields) {
                result.push_back(sf.first);
            }
            for (auto& vf : _vfields) {
                result.push_back(vf.first);
            }
            return result;
        }
        
        inline void setScalar(const std::string& fname, const scalar_type val) {
            _sfields.at(fname) = val;
        }
        
        inline void setVector(const std::string& fname, const vfield_type& val) {
            _vfields.at(fname) = val;
        }
        
        inline void setVector(const std::string& fname, const Vec<ndim>& val) {
            _vfields.at(fname) = val.toStdVec();
        }
        
        inline scalar_type getScalar(const std::string& fname) const {return _sfields.at(fname);}
        
        inline vfield_type getVector(const std::string& fname) const {return _vfields.at(fname);}
        
        
        
        
    protected:

        Face() {}
    
        ind_vec _interiorInds;
        ind_vec _vertInds;
        ind_vec _edgeInds;
        index_type _parent;
        std::array<index_type, 4> _kids;
        scalar_type _area;
        
        std::map<std::string, scalar_type> _sfields;
        std::map<std::string, std::vector<scalar_type>> _vfields;
};

template <int ndim> class QuadFace : public Face<ndim> {
    public:
        QuadFace(const ind_vec intrs, const ind_vec& verts, const ind_vec& edges, 
            const index_type prnt, const scalar_type ar = 0.0) : 
                Face<ndim>(intrs, verts, edges, prnt, ar) {}
};

template <int ndim> class TriFace : public Face<ndim> {
    public:
        TriFace(const ind_vec& intrs, const ind_vec& verts, const ind_vec& edges, 
            const index_type prnt, const scalar_type ar = 0.0) : 
                Face<ndim>(intrs, verts, edges, prnt, ar) {}
};

template <int ndim> class QuadCubicFace : public Face<ndim> {
    public: 
        QuadCubicFace(const ind_vec& intInds, const ind_vec& vertInds, const ind_vec& edgeInds, 
            const index_type pt, const scalar_type ar=0.0) : Face<ndim>(intInds, vertInds, edgeInds, pt, ar) {}
};

}
}
#endif
