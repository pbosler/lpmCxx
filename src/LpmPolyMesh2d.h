#ifndef _LPM_POLYMESH_2D_H_
#define _LPM_POLYMESH_2D_H_

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmCoords.h"
#include "LpmEdges.h"
#include "LpmFaces.h"
#include "LpmMeshSeed.h"
#include "LpmXyzVector.h"
#include "LpmLogger.h"
#include "LpmVtkFileIO.h"
#include <memory>
#include <cmath>
#include <string>


namespace Lpm {

class PolyMesh2d {
    public:
        friend class MeshedParticles;
    
        PolyMesh2d(MeshSeed& seed, const int maxRecursionLevel, const bool isLagrangian = false, const scalar_type domainRadius = 1.0,
            const int prank = 0);
            
        inline GeometryType geometry() const {return _coords->geometry();}
        
        inline std::shared_ptr<Coords> getPhysCoords() const {return _coords;}
        inline std::shared_ptr<Coords> getLagCoords() const {return _lagCoords;}
        
        inline std::shared_ptr<Faces> getFaces() const {return _faces;}
        inline std::shared_ptr<Coords> makeCoordsFromFaces() const {return _faces->makeCoordsFromCentroids();}
        inline std::shared_ptr<Field> areaWeightField() const {return _faces->centroidAreas();}
    
        inline index_type nVertices() const {return _coords->n();}
        inline index_type nEdges() const {return _edges->n();}
        inline index_type nFaces() const {return _faces->n();}
        
        inline bool faceIsDivided(const index_type faceInd) const {return _faces->hasChildren(faceInd);}
        inline XyzVector faceCentroid(const index_type faceInd, const bool lagrangian = false) const 
            {return _faces->centroid(faceInd, lagrangian);}
        inline scalar_type faceArea(const index_type faceInd) const {return _faces->area(faceInd);}
        
        inline bool isLagrangian() const {return (_lagCoords.use_count() > 0);}
        
        inline index_type nLeafFaces() const {return _faces->nLeaves();}
        inline index_type nLeafEdges() const {return _edges->nLeaves();}
        
        inline scalar_type surfaceArea() const {return _faces->surfaceArea();}
        inline scalar_type avgMeshSize() const {return std::sqrt(surfaceArea() / _faces->n());}
        
        inline scalar_type maxMeshSize() const {return _edges->maxLength();}
        inline scalar_type minMeshSize() const {return _edges->minLength();}
        
        inline std::vector<index_type> faceVertices(const index_type faceInd) const 
            {return _faces->vertexIndices(faceInd);}
        
        inline std::vector<XyzVector> getCoordVecs(const std::vector<index_type> crdInds) const 
            {return _coords->getVectors(crdInds);}
    
        index_type locateFaceContainingPoint(const XyzVector& queryPt) const;
        
        inline std::vector<index_type> adjacentFaceIndices(const index_type faceInd) const 
            {return _faces->ccwAdjacentFaces(faceInd);}
        
        void writeToVTKFile(const std::string& fname, const std::string& desc = "") const;
        
        std::shared_ptr<Field> createFaceAreaField() const;
        
        inline void setLogProc(const int rank) {log->setProcRank(rank);}
    protected:
        
        index_type walkSearch(const XyzVector& queryPt, index_type startIndex) const;
        index_type treeSearch(const XyzVector& queryPt, index_type startIndex) const;
        index_type nearestRootFace(const XyzVector& queryPt) const;
    
        std::shared_ptr<Coords> _coords;
        std::shared_ptr<Coords> _lagCoords;
        std::shared_ptr<Edges> _edges;
        std::shared_ptr<Faces> _faces;
        
        int nRootFaces;
        bool lagrangian;
        
        static std::unique_ptr<Logger> log;

};

}

#endif
