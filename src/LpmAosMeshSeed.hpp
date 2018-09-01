#ifndef LPM_AOS_MESH_SEED_HPP
#define LPM_AOS_MESH_SEED_HPP

#include <string>
#include <vector>
#include "LpmTypeDefs.h"
#include "LpmConfig.h"
#include "LpmAosTypes.hpp"
#include <iostream>

namespace Lpm {

namespace Aos {

class MeshSeed {
    public:
        typedef std::vector<Vec<2>> crd2_vec_type;
        typedef std::vector<Vec<3>> crd3_vec_type;
        typedef std::vector<index_type> ind_vec_type;
    
        //void initFromFile();
        
        virtual ~MeshSeed() {};
        
        virtual index_type nFacesAfterUniformRefinement(const index_type maxRecursion) const = 0;
        virtual index_type nEdgesAfterUniformRefinement(const index_type nverts, const index_type nfaces) const = 0;
        virtual index_type nVerticesAfterUniformRefinement(const index_type maxRecursion) const = 0;
        virtual index_type nRootFaces() const = 0;
        virtual std::string idString() const = 0;
        virtual index_type ndim() const = 0;
    
        std::string infoString() const;
        
        void initFromFile();
    
    protected:
        MeshSeed(const std::string fname, const int ndim, const int ncoords, const int nedges, const int nfaces, 
            const int nedgesperface) : _fname(fname), _ndim(ndim), _nCrds(ncoords), _nEdges(nedges), 
            _nFaces(nfaces), _nEdgesPerFace(nedgesperface) {
            //_nCrds+=_nFaces;
        }
        
        int _nCrds;
        int _nEdges;
        int _nFaces;
        std::string _fname;
        short _ndim;
        short _nEdgesPerFace;
        
        crd2_vec_type r2vertexCrds;
        crd3_vec_type r3vertexCrds;
        ind_vec_type edgeOrigs;
        ind_vec_type edgeDests;
        ind_vec_type edgeLefts;
        ind_vec_type edgeRights;
        std::vector<ind_vec_type> faceVerts;
        std::vector<ind_vec_type> faceEdges;
        
};

class TriHexSeed : public MeshSeed {
    public:
        TriHexSeed() : MeshSeed("triHexSeed.dat", 2, 7, 12, 6, 3) {}
        index_type nFacesAfterUniformRefinement(const index_type maxRecursion) const;
        index_type nEdgesAfterUniformRefinement(const index_type nverts, const index_type nfaces) const;
        index_type nVerticesAfterUniformRefinement(const index_type maxRecursion) const;
        inline index_type nRootFaces() const {return 6;}
        inline std::string idString() const {return "triHexPlane";}
        inline index_type ndim() const {return 2;}
};

class QuadRectSeed : public MeshSeed {
    public:
        QuadRectSeed() : MeshSeed("quadRectSeed.dat", 2, 9, 12, 4, 4) {}
        index_type nFacesAfterUniformRefinement(const index_type maxRecursion) const;
        index_type nEdgesAfterUniformRefinement(const index_type nverts, const index_type nfaces) const;
        index_type nVerticesAfterUniformRefinement(const index_type maxRecursion) const;
        inline index_type nRootFaces() const {return 4;}
        inline std::string idString() const {return "quadRectPlane";}
        inline index_type ndim() const {return 2;}
};

class IcosTriSphereSeed : public MeshSeed {
    public:
        IcosTriSphereSeed() : MeshSeed("icosTriSphereSeed.dat", 3, 12, 30, 20, 3) {}
        index_type nFacesAfterUniformRefinement(const index_type maxRecursion) const;
        index_type nEdgesAfterUniformRefinement(const index_type nverts, const index_type nfaces) const;
        index_type nVerticesAfterUniformRefinement(const index_type maxRecursion) const;
        inline index_type nRootFaces() const {return 20;}
        inline std::string idString() const {return "icosTriSphere";}
        inline index_type ndim() const {return 3;}
};

class CubedSphereSeed : public MeshSeed {
    public:
        CubedSphereSeed() : MeshSeed("cubedSphereSeed.dat", 3, 8, 12, 6, 4) {}
        index_type nFacesAfterUniformRefinement(const index_type maxRecursion) const;
        index_type nEdgesAfterUniformRefinement(const index_type nverts, const index_type nfaces) const;
        index_type nVerticesAfterUniformRefinement(const index_type maxRecursion) const;
        inline index_type nRootFaces() const {return 6;}
        inline std::string idString() const {return "cubedSphere";}
        inline index_type ndim() const {return 3;}
};

}
}
#endif
