#ifndef _LPM_TREECODE_H_
#define _LPM_TREECODE_H_

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmXyzVector.h"
#include "LpmCoords.h"
#include "LpmLogger.h"
#include <memory>
#include <vector>
#include <string>

namespace Lpm {


struct Box3d {
    Box3d() : xmin(0.0), xmax(0.0), ymin(0.0), ymax(0.0), zmin(0.0), zmax(0.0) {};
    Box3d(const scalar_type x0, const scalar_type x1, const scalar_type y0, const scalar_type y1, const scalar_type z0,
        const scalar_type z1) : xmin(x0), xmax(x1), ymin(y0), ymax(y1), zmin(z0), zmax(z1) {};
    
    inline scalar_type volume() const {return (xmax - xmin) * (ymax - ymin) * (zmax - zmin);}
    inline scalar_type area2d() const {return (xmax - xmin) * (ymax - ymin);}
    
    inline XyzVector centroid() const {return XyzVector(0.5 * (xmax - xmin), 0.5 * (ymax - ymin), 0.5 * (zmax - zmin));}
    
    inline bool containsPoint(const XyzVector& vec) const {return (xmin <= vec.x && vec.x < xmax) &&
                                                                  (ymin <= vec.y && vec.y < ymax) &&
                                                                  (zmin <= vec.z && vec.z < zmax);}
    
    scalar_type longestEdge() const;
    scalar_type shortestEdge() const;
    scalar_type aspectRatio() const;
    
    scalar_type edgelength(const int dim) const;
    
    std::vector<Box3d> bisectAll() const;
    std::vector<Box3d> bisectAlongDims(const bool* dims) const;
    
    scalar_type radius() const;
    
    /// minimum x-coordinate of this box
    scalar_type xmin;
    /// maximum x-coordinate of this box
    scalar_type xmax;
    /// minimum y-coordinate of this box
    scalar_type ymin;
    /// maximum y-coordinate of this box
    scalar_type ymax;
    /// minimum z-coordinate of this box
    scalar_type zmin;
    /// maximum z-coordinate of this box
    scalar_type zmax;
};

struct Treenode {
    Treenode(const Box3d& bbox, const index_type nCrds = 0, const scalar_type maxAspectRatio = 1.0);
    Treenode(const Box3d& bbox, const std::shared_ptr<Treenode>& pparent = NULL, 
             const std::vector<index_type>& crdInds, const scalar_type maxAspectRatio = 1.0);

    inline void setLogProc(const int prank) {log->setProcRank(prank);}

    Box3d box;
    scalar_type maxAspectRatio;
    int level;
    
    std::vector<index_type> coordsContained;
    
    inline index_type nCoords() const {return coordsContained.size();}
    
    void shrinkBox(const std::shared_ptr<Coords>& crds);
    
    /// ptr to parent
    std::shared_ptr<Treenode> parent;
    
    /// ptrs to child boxes
    std::vector<std::shared_ptr<Treenode>> children;
    
    bool hasChildren() const {return (!children.empty());}
    
    void writePoints(std::ofstream& os) const;
    
    static std::unique_ptr<Logger> log;
};

index_type nTreenodes(const std::shared_ptr<Treenode>& node);

void generateTree(std::shared_ptr<Treenode>& node, std::shared_ptr<Coords>& crds, const index_type maxCoordsPerNode);

void writeTreeToVtk(const std::string& filename, const std::string& desc = "", const std::shared_ptr<Treenode>& node);
}

#endif
