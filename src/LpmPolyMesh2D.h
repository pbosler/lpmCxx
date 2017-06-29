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
        enum GeometryType {PLANAR_GEOMETRY, SPHERICAL_SURFACE_GEOMETRY};
        typedef std::tuple<index_type, index_type, index_type, index_type> quad_index_type;
    
        PolyMesh2d(MeshSeed& seed, const int maxRecursionLevel, const bool isLagrangian = false, const scalar_type domainRadius = 1.0);
    
        inline index_type nVertices() const {return coords->n();}
        inline index_type nEdges() const {return edges->n();}
        inline index_type nFaces() const {return faces->n();}
        
        inline scalar_type surfaceArea() const {return faces->surfaceArea();}
        inline scalar_type avgMeshSize() const {return std::sqrt(surfaceArea() / faces->n());}
        
        inline scalar_type maxMeshSize() const {return edges->maxLength();}
        inline scalar_type minMeshSize() const {return edges->minLength();}
    
        index_type locateFaceContainingPoint(const XyzVector& queryPt) const;
        
        std::vector<index_type> nNearbyCoordinates(const XyzVector& queryPt, const index_type n) const;
        
        void writeToVTKFile(const std::string& fname, const std::string& desc = "") const;
    protected:
        
        index_type walkSearch(const XyzVector& queryPt, index_type startIndex) const;
        index_type treeSearch(const XyzVector& queryPt, index_type startIndex) const;
        index_type nearestRootFace(const XyzVector& queryPt) const;
    
        std::shared_ptr<Coords> coords;
        std::shared_ptr<Coords> lagCoords;
        std::shared_ptr<Edges> edges;
        std::shared_ptr<Faces> faces;
        
        GeometryType geometry;
        int nRootFaces;
        bool lagrangian;
        
        static std::unique_ptr<Logger> log;

};

}

#endif
