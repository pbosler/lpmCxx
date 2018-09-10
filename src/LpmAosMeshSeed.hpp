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

template <int ndim> class MeshSeed {
    public:
        typedef std::vector<Vec<ndim>> crd_vec_type;
        typedef std::vector<index_type> ind_vec_type;
        
        crd_vec_type crds;
        ind_vec_type edgeOrigs;
        ind_vec_type edgeDests;
        ind_vec_type edgeLefts;
        ind_vec_type edgeRights;
        std::vector<ind_vec_type> edgeInteriors;
        std::vector<ind_vec_type> faceVerts;
        std::vector<ind_vec_type> faceEdges;
        std::vector<ind_vec_type> faceInteriors;
        
        virtual ~MeshSeed() {};
        
        virtual index_type nFacesAfterUniformRefinement(const index_type maxRecursion) const = 0;
        virtual index_type nEdgesAfterUniformRefinement(const index_type nverts, const index_type nfaces, const index_type recursionLevel) const = 0;
        virtual index_type nVerticesAfterUniformRefinement(const index_type maxRecursion) const = 0;
        virtual index_type nRootFaces() const = 0;
        virtual std::string idString() const = 0;
//         virtual index_type ndim() const = 0;
        virtual GeometryType geometryType() const = 0;
        virtual FaceType faceType() const = 0;
        virtual bool faceCrdsIncluded() const = 0;
        
        inline index_type nVerts() const {return _nVerts;}
        inline index_type nCrds() const {return _nCrds;}
        inline index_type nEdges() const {return _nEdges;}
        inline index_type nFaces() const {return _nFaces;}
        inline short nEdgesPerFace() const {return _nEdgesPerFace;}
        
        std::string infoString() const;
        
        
        
        void determineMaxAllocations(index_type& nv, index_type& nf, index_type& ne, const index_type maxRec) const;
    
    protected:
        MeshSeed(const std::string fname, const int nverts, const int ncoords, const int nedges, const int nfaces, 
            const int nedgesperface) : _fname(fname), _ndim(ndim), _nVerts(nverts), _nCrds(ncoords), _nEdges(nedges), 
            _nFaces(nfaces), _nEdgesPerFace(nedgesperface) {
        }
        
        void initFromFile();
        
        int _nCrds;
        int _nVerts;
        int _nEdges;
        int _nFaces;
        std::string _fname;
        short _ndim;
        short _nEdgesPerFace;
};

class TriHexSeed : public MeshSeed<2> {
    public:
        TriHexSeed() : MeshSeed("triHexSeed.dat", 7, 13, 12, 6, 3) {initFromFile();}
        index_type nFacesAfterUniformRefinement(const index_type maxRecursion) const;
        index_type nEdgesAfterUniformRefinement(const index_type nverts, const index_type nfaces, const index_type recursionLevel) const;
        index_type nVerticesAfterUniformRefinement(const index_type maxRecursion) const;
        inline index_type nRootFaces() const {return 6;}
        inline std::string idString() const {return "triHexPlane";}
        inline GeometryType geometryType() const {return PLANAR_GEOMETRY;}
        inline FaceType faceType() const {return TRI;}
        inline bool faceCrdsIncluded() const {return true;}
};

class QuadRectSeed : public MeshSeed<2> {
    public:
        QuadRectSeed() : MeshSeed("quadRectSeed.dat", 9, 13, 12, 4, 4) {initFromFile();}
        index_type nFacesAfterUniformRefinement(const index_type maxRecursion) const;
        index_type nEdgesAfterUniformRefinement(const index_type nverts, const index_type nfaces, const index_type recursionLevel) const;
        index_type nVerticesAfterUniformRefinement(const index_type maxRecursion) const;
        inline index_type nRootFaces() const {return 4;}
        inline std::string idString() const {return "quadRectPlane";}
        inline GeometryType geometryType() const {return PLANAR_GEOMETRY;}
        inline FaceType faceType() const {return QUAD;}
        inline bool faceCrdsIncluded() const {return true;}
};

class IcosTriSphereSeed : public MeshSeed<3> {
    public:
        IcosTriSphereSeed() : MeshSeed("icosTriSphereSeed.dat", 12, 32, 30, 20, 3) {initFromFile();}
        index_type nFacesAfterUniformRefinement(const index_type maxRecursion) const;
        index_type nEdgesAfterUniformRefinement(const index_type nverts, const index_type nfaces, const index_type recursionLevel) const;
        index_type nVerticesAfterUniformRefinement(const index_type maxRecursion) const;
        inline index_type nRootFaces() const {return 20;}
        inline std::string idString() const {return "icosTriSphere";}
        inline GeometryType geometryType() const {return SPHERICAL_SURFACE_GEOMETRY;}
        inline FaceType faceType() const {return TRI;}
        inline bool faceCrdsIncluded() const {return true;}
};

class CubedSphereSeed : public MeshSeed<3> {
    public:
        CubedSphereSeed() : MeshSeed("cubedSphereSeed.dat", 8, 14, 12, 6, 4) {initFromFile();}
        index_type nFacesAfterUniformRefinement(const index_type maxRecursion) const;
        index_type nEdgesAfterUniformRefinement(const index_type nverts, const index_type nfaces, const index_type recursionLevel) const;
        index_type nVerticesAfterUniformRefinement(const index_type maxRecursion) const;
        inline index_type nRootFaces() const {return 6;}
        inline std::string idString() const {return "cubedSphere";}
        inline GeometryType geometryType() const {return SPHERICAL_SURFACE_GEOMETRY;}
        inline FaceType faceType() const {return QUAD;}
        inline bool faceCrdsIncluded() const {return true;}
};

class QuadCubicSeed : public MeshSeed<2> {
    public:
    QuadCubicSeed() : MeshSeed("quadCubicSeed.dat", 37, 49, 12, 4, 4) {initFromFile();}
    index_type nFacesAfterUniformRefinement(const index_type maxRecursion) const;
    index_type nEdgesAfterUniformRefinement(const index_type nverts, const index_type nfaces, const index_type recursionLevel) const;
    index_type nVerticesAfterUniformRefinement(const index_type maxRecursion) const;
    inline index_type nRootFaces() const {return 6;}
    inline std::string idString() const {return "cubedSphere";}
    inline GeometryType geometryType() const {return SPHERICAL_SURFACE_GEOMETRY;}
    inline FaceType faceType() const {return QUAD;}
    inline bool faceCrdsIncluded() const {return true;}
};

}
}
#endif
