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
        void initMeshFromSeed(std::shared_ptr<Coords> crds, std::shared_ptr<Edges> edges, std::shared_ptr<Faces> faces);
        
        std::string infoString() const;

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
        TriHexSeed(const scalar_type rFac = 1.0) : MeshSeed("triHexSeed.dat", 2, 7, 12, 6, 3), radialFactor(rFac) {};
        
    protected:
        scalar_type radialFactor;
//         void readSeedFile();
};

class QuadRectSeed : public MeshSeed {
    public:
        QuadRectSeed(const scalar_type xFac = 1.0, const scalar_type yFac = 1.0) : 
            MeshSeed("quadRectSeed.dat", 2, 9, 12, 4, 4), 
            xFactor(xFac), yFactor(yFac) {};
    protected:
//         void readSeedFile();
        scalar_type xFactor;
        scalar_type yFactor;
};

class IcosTriSphereSeed : public MeshSeed {
    public:
        IcosTriSphereSeed(const scalar_type sphereRadius = 1.0) : MeshSeed("icosTriSphereSeed.dat", 3, 12, 30, 20, 3), 
            radialFactor(sphereRadius) {};
    
        void initFromSeed(std::shared_ptr<Coords> crds, std::shared_ptr<Edges> edges, std::shared_ptr<Faces> faces);
        
    protected:
//         void readSeedFile();
        scalar_type radialFactor;
};

class CubedSphereSeed : public MeshSeed {
    public:
        CubedSphereSeed(const scalar_type sphereRadius = 1.0) : MeshSeed("cubedSphereSeed.dat", 3, 8, 12, 6, 4), radialFactor(sphereRadius) {};
    
    protected:
//         void readSeedFile();
        scalar_type radialFactor;
};

}



#endif
