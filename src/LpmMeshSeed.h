#ifndef _LPM_MESH_SEED_H_
#define _LPM_MESH_SEED_H_

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmFaces.h"
#include "LpmEdges.h"
#include "LpmCoords.h"
#include "LpmLogger.h"
#include "LpmXyzVector.h"
#include <string>
#include <memory>

namespace Lpm {

class MeshSeed {
    public:
        void initMeshFromSeed(std::shared_ptr<Coords> crds, std::shared_ptr<Edges> edges, std::shared_ptr<Faces> faces, 
            const scalar_type domainRadius = 1.0);
        
        std::string infoString() const;
        
        virtual index_type nFaces(const int recursionLevel) const = 0;
        virtual index_type nVertices(const int recursionLevel) const = 0;
        virtual index_type nEdges(const index_type nVerts, const index_type nFaces) const = 0;
        virtual index_type nRootFaces() const = 0;
        virtual std::string idString() const = 0;

        inline void setLogProc(const int rank) {log->setProcRank(rank);}
    protected:
        MeshSeed(const std::string fname, const int nDim, const index_type nCrds, const index_type nEdges, const index_type nFaces, 
            const int nEdgesPerFace) : 
           _nCoords(nCrds), _nEdges(nEdges), _nFaces(nFaces), _fname(fname), _nDim(nDim), _nEdgesPerFace(nEdgesPerFace) {};
    
        void readSeedFile();
           
        std::vector<XyzVector> vertCrds;
        std::vector<index_type> edgeOrigs;
        std::vector<index_type> edgeDests;
        std::vector<index_type> edgeLefts;
        std::vector<index_type> edgeRights;
        std::vector<std::vector<index_type>> faceEdges;
           
        index_type _nCoords;
        index_type _nEdges;
        index_type _nFaces;
        std::string _fname;
        int _nDim;
        int _nEdgesPerFace;
        
        static std::unique_ptr<Logger> log;
};

class TriHexSeed : public MeshSeed {
    public:
        TriHexSeed() : MeshSeed("triHexSeed.dat", 2, 7, 12, 6, 3) {};
        index_type nFaces(const int recursionLevel) const;
        index_type nVertices(const int recursionLevel) const;
        index_type nEdges(const index_type nverts, const index_type nfaces) const;
        inline index_type nRootFaces() const {return 6;}
        inline std::string idString() const {return "triHexPlane";}
};

class QuadRectSeed : public MeshSeed {
    public:
        QuadRectSeed() : MeshSeed("quadRectSeed.dat", 2, 9, 12, 4, 4) {};
        index_type nFaces(const int recursionLevel) const;
        index_type nVertices(const int recursionLevel) const;
        index_type nEdges(const index_type nverts, const index_type nfaces) const;
        inline index_type nRootFaces() const {return 4;}
        inline std::string idString() const {return "quadRectPlane";}
};

class IcosTriSphereSeed : public MeshSeed {
    public:
        IcosTriSphereSeed() : MeshSeed("icosTriSphereSeed.dat", 3, 12, 30, 20, 3) {};
        index_type nFaces(const int recursionLevel) const;
        index_type nVertices(const int recursionLevel) const;
        index_type nEdges(const index_type nverts, const index_type nfaces) const;
        inline index_type nRootFaces() const {return 20;}
        inline std::string idString() const {return "icosTriSphere";}
};

class CubedSphereSeed : public MeshSeed {
    public:
        CubedSphereSeed() : MeshSeed("cubedSphereSeed.dat", 3, 8, 12, 6, 4) {};
        index_type nFaces(const int recursionLevel) const;
        index_type nVertices(const int recursionLevel) const;
        index_type nEdges(const index_type nverts, const index_type nfaces) const;
        inline index_type nRootFaces() const {return 6;}
        inline std::string idString() const {return "cubedSphere";}
};

}



#endif
