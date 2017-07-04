#ifndef _LPM_HEADER_TEMPLATE_H_
#define _LPM_HEADER_TEMPLATE_H_

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
        typedef std::tuple<index_type, index_type, index_type, index_type> quad_index_type;
        friend class MeshedParticles;
    
        PolyMesh2d(MeshSeed& seed, const int maxRecursionLevel, const bool isLagrangian = false, const scalar_type domainRadius = 1.0);
    
        inline GeometryType geometry() const {return coords->geometry();}
        
        inline std::shared_ptr<Coords> getPhysCoords() const {return coords;}
        inline std::shared_ptr<Coords> getLagCoords() const {return lagCoords;}
    
        inline index_type nVertices() const {return coords->n();}
        inline index_type nEdges() const {return edges->n();}
        inline index_type nFaces() const {return faces->n();}
        
        inline bool faceIsDivided(const index_type faceInd) const {return faces->hasChildren(faceInd);}
        inline XyzVector faceCentroid(const index_type faceInd, const bool lagrangian = false) const 
            {return faces->centroid(faceInd, lagrangian);}
        inline scalar_type faceArea(const index_type faceInd) const {return faces->area(faceInd);}
        
        inline bool isLagrangian() const {return (lagCoords.use_count() > 0);}
        
        inline index_type nLeafFaces() const {return faces->nLeaves();}
        inline index_type nLeafEdges() const {return edges->nLeaves();}
        
        inline scalar_type surfaceArea() const {return faces->surfaceArea();}
        inline scalar_type avgMeshSize() const {return std::sqrt(surfaceArea() / faces->n());}
        
        inline scalar_type maxMeshSize() const {return edges->maxLength();}
        inline scalar_type minMeshSize() const {return edges->minLength();}
        
        inline std::vector<index_type> faceVertices(const index_type faceInd) const 
            {return faces->vertexIndices(faceInd);}
        
        inline std::vector<XyzVector> getCoordVecs(const std::vector<index_type> crdInds) const 
            {return coords->getVectors(crdInds);}
    
        index_type locateFaceContainingPoint(const XyzVector& queryPt) const;
        
        inline std::vector<index_type> adjacentFaceIndices(const index_type faceInd) const 
            {return faces->ccwAdjacentFaces(faceInd);}
        
        void writeToVTKFile(const std::string& fname, const std::string& desc = "") const;
        
        std::shared_ptr<Field> createFaceAreaField() const;
        
    protected:
        
        index_type walkSearch(const XyzVector& queryPt, index_type startIndex) const;
        index_type treeSearch(const XyzVector& queryPt, index_type startIndex) const;
        index_type nearestRootFace(const XyzVector& queryPt) const;
    
        std::shared_ptr<Coords> coords;
        std::shared_ptr<Coords> lagCoords;
        std::shared_ptr<Edges> edges;
        std::shared_ptr<Faces> faces;
        
        int nRootFaces;
        bool lagrangian;
        
        static std::unique_ptr<Logger> log;

};

}

#endif
