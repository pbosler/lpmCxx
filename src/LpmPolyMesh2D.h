#ifndef _LPM_HEADER_TEMPLATE_H_
#define _LPM_HEADER_TEMPLATE_H_

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmCoords.h"
#include "LpmEdges.h"
#include "LpmFaces.h"
#include "LpmMeshSeed.h"
#include "LpmXyzVector.h"
#include <memory>
#include <cmath>

namespace Lpm {

class PolyMesh2d {
    public:
        inline index_type nVertices() const {return coords->n();}
        inline index_type nEdges() const {return edges->n();}
        inline index_type nFaces() const {return faces->n();}
        
        inline scalar_type surfaceArea() const {return faces->surfaceArea();}
        inline scalar_type avgMeshSize() const {return std::sqrt(surfaceArea() / faces->n());}
        
        scalar_type maxMeshSize() const;
        scalar_type minMeshSize() const;
    
        index_type locatePointInFace(const XyzVector& queryPt) const;
        
        std::vector<XyzVector> nNearbyCoordinates(const XyzVector& queryPt, const index_type n) const;
        
    protected:
        index_type estNVertices(const int recursionLevel) const;
        index_type estNFaces(const int recursionLevel) const;
        index_type estNnEdges(const index_type nVerts, const index_type nFaces) const;
        
        index_type walkSearch(const XyzVector& queryPt, index_type startIndex) const;
        index_type treeSearch(const XyzVector& queryPt, index_type startIndex) const;
        index_type nearestRootFace(const XyzVector& queryPt) const;
    
        std::shared_ptr<Coords> coords;
        std::shared_ptr<Coords> lagCoords;
        std::shared_ptr<Edges> edges;
        std::shared_ptr<Faces> faces;
        
        std::unique_ptr<MeshSeed> mSeed;
};

}

#endif
