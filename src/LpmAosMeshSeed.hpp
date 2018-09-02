#ifndef LPM_AOS_MESH_SEED_HPP
#define LPM_AOS_MESH_SEED_HPP

#include <string>
#include <vector>
#include "LpmTypeDefs.h"
#include "LpmConfig.h"
#include "LpmAosTypes.hpp"
#include "LpmAosParticleFactory.hpp"
#include "LpmAosEdgeFactory.hpp"
#include "LpmAosFaceFactory.hpp"
#include <iostream>
#include <memory>

namespace Lpm {

namespace Aos {

class MeshSeed {
    public:
        typedef std::vector<Vec<2>> crd2_vec_type;
        typedef std::vector<Vec<3>> crd3_vec_type;
        typedef std::vector<index_type> ind_vec_type;
        
        crd2_vec_type r2Crds;
        crd3_vec_type r3Crds;
        ind_vec_type edgeOrigs;
        ind_vec_type edgeDests;
        ind_vec_type edgeLefts;
        ind_vec_type edgeRights;
        std::vector<ind_vec_type> faceVerts;
        std::vector<ind_vec_type> faceEdges;
        
        virtual ~MeshSeed() {};
        
        virtual index_type nFacesAfterUniformRefinement(const index_type maxRecursion) const = 0;
        virtual index_type nEdgesAfterUniformRefinement(const index_type nverts, const index_type nfaces) const = 0;
        virtual index_type nVerticesAfterUniformRefinement(const index_type maxRecursion) const = 0;
        virtual index_type nRootFaces() const = 0;
        virtual std::string idString() const = 0;
        virtual index_type ndim() const = 0;
        virtual GeometryType geometryType() const = 0;
        virtual FaceType faceType() const = 0;
        virtual bool faceCrdsIncluded() const = 0;
        
        inline index_type nVerts() const {return _nVerts;}
        inline index_type nCrds() const {return _nCrds;}
        inline index_type nEdges() const {return _nEdges;}
        inline index_type nFaces() const {return _nFaces;}
        inline short nEdgesPerFace() const {return _nEdgesPerFace;}
        
        std::string infoString() const;
        
        void initFromFile();
        
        void determineMaxAllocations(index_type& nv, index_type& nf, index_type& ne, const index_type maxRec) const;
    
    protected:
        MeshSeed(const std::string fname, const int ndim, const int nverts, const int ncoords, const int nedges, const int nfaces, 
            const int nedgesperface) : _fname(fname), _ndim(ndim), _nVerts(nverts), _nCrds(ncoords), _nEdges(nedges), 
            _nFaces(nfaces), _nEdgesPerFace(nedgesperface) {
        }
        
        int _nCrds;
        int _nVerts;
        int _nEdges;
        int _nFaces;
        std::string _fname;
        short _ndim;
        short _nEdgesPerFace;
};

class TriHexSeed : public MeshSeed {
    public:
        TriHexSeed() : MeshSeed("triHexSeed.dat", 2, 7, 13, 12, 6, 3) {}
        index_type nFacesAfterUniformRefinement(const index_type maxRecursion) const;
        index_type nEdgesAfterUniformRefinement(const index_type nverts, const index_type nfaces) const;
        index_type nVerticesAfterUniformRefinement(const index_type maxRecursion) const;
        inline index_type nRootFaces() const {return 6;}
        inline std::string idString() const {return "triHexPlane";}
        inline index_type ndim() const {return 2;}
        inline GeometryType geometryType() const {return PLANAR_GEOMETRY;}
        inline FaceType faceType() const {return TRI;}
        inline bool faceCrdsIncluded() const {return true;}
};

class QuadRectSeed : public MeshSeed {
    public:
        QuadRectSeed() : MeshSeed("quadRectSeed.dat", 2, 9, 13, 12, 4, 4) {}
        index_type nFacesAfterUniformRefinement(const index_type maxRecursion) const;
        index_type nEdgesAfterUniformRefinement(const index_type nverts, const index_type nfaces) const;
        index_type nVerticesAfterUniformRefinement(const index_type maxRecursion) const;
        inline index_type nRootFaces() const {return 4;}
        inline std::string idString() const {return "quadRectPlane";}
        inline index_type ndim() const {return 2;}
        inline GeometryType geometryType() const {return PLANAR_GEOMETRY;}
        inline FaceType faceType() const {return QUAD;}
        inline bool faceCrdsIncluded() const {return true;}
};

class IcosTriSphereSeed : public MeshSeed {
    public:
        IcosTriSphereSeed() : MeshSeed("icosTriSphereSeed.dat", 3, 12, 12, 30, 20, 3) {}
        index_type nFacesAfterUniformRefinement(const index_type maxRecursion) const;
        index_type nEdgesAfterUniformRefinement(const index_type nverts, const index_type nfaces) const;
        index_type nVerticesAfterUniformRefinement(const index_type maxRecursion) const;
        inline index_type nRootFaces() const {return 20;}
        inline std::string idString() const {return "icosTriSphere";}
        inline index_type ndim() const {return 3;}
        inline GeometryType geometryType() const {return SPHERICAL_SURFACE_GEOMETRY;}
        inline FaceType faceType() const {return TRI;}
        inline bool faceCrdsIncluded() const {return false;}
};

class CubedSphereSeed : public MeshSeed {
    public:
        CubedSphereSeed() : MeshSeed("cubedSphereSeed.dat", 3, 8, 14, 12, 6, 4) {}
        index_type nFacesAfterUniformRefinement(const index_type maxRecursion) const;
        index_type nEdgesAfterUniformRefinement(const index_type nverts, const index_type nfaces) const;
        index_type nVerticesAfterUniformRefinement(const index_type maxRecursion) const;
        inline index_type nRootFaces() const {return 6;}
        inline std::string idString() const {return "cubedSphere";}
        inline index_type ndim() const {return 3;}
        inline GeometryType geometryType() const {return SPHERICAL_SURFACE_GEOMETRY;}
        inline FaceType faceType() const {return QUAD;}
        inline bool faceCrdsIncluded() const {return true;}
};

}
}
#endif
