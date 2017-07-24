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
    
    inline XyzVector centroid() const {return XyzVector(0.5 * (xmax + xmin), 0.5 * (ymax + ymin), 0.5 * (zmax + zmin));}
    
    inline bool containsPoint(const XyzVector& vec) const {return (xmin <= vec.x && vec.x < xmax) &&
                                                                  (ymin <= vec.y && vec.y < ymax) &&
                                                                  (zmin <= vec.z && vec.z < zmax);}
    
    scalar_type longestEdge() const;
    scalar_type shortestEdge() const;
    scalar_type aspectRatio() const;
    
    scalar_type edgeLength(const int dim) const;
    
    std::vector<Box3d> bisectAll() const;
    std::vector<Box3d> bisectAlongDims(const bool* dims) const;
    
    std::string infoString() const;
    
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

class Octree {
    struct Node {
        Node();
        Node(const Coords* crds, const scalar_type maxAspectRatio = 1.0);
        Node(const Box3d& bbox, const Node* pparent = NULL, 
                 const std::vector<index_type>& crdInds = std::vector<index_type>(), const scalar_type maxAspectRatio = 1.0);

        Box3d box;
        scalar_type maxAspectRatio;
        int level;
        /// ptr to parent
        Node* parent;
        /// ptrs to child boxes
        std::vector<std::unique_ptr_ptr<Node>> kids;
        std::vector<index_type> coordsContained;
    
        inline index_type nCoords() const {return coordsContained.size();}
    
        void shrinkBox(const Coords* crds);
    
        bool hasKids() const {return (!kids.empty());}
    
        void writePoints(std::ofstream& os) const;
        inline void setLogProc(const int prank) {log->setProcRank(prank);}
    
        std::string infoString() const;
    };
    std::unique_ptr<Node> root;
    static std::unique_ptr<Logger> log;
    inline void setLogProc(const int prank) {log->setProcRank(prank);}
    std::string infoString() const;
};

int treeDepth(const std::shared_ptr<Treenode> node);

index_type nTreenodes(const std::shared_ptr<Treenode> node);

void generateTree(std::shared_ptr<Treenode> node, std::shared_ptr<Coords> crds, const index_type maxCoordsPerNode);

void writeTreeToVtk(const std::string& filename, const std::string& desc = "", const std::shared_ptr<Treenode> node = NULL);
void writeVTKPoints(std::ofstream& os, const std::shared_ptr<Treenode> node);
void writeVtkCells(std::ofstream& os, const std::shared_ptr<Treenode> node, index_type& vertIndex);
void writeVtkCellType(std::ofstream& os, const std::shared_ptr<Treenode> node);
void writeVtkLevelData(std::ofstream& os, const std::shared_ptr<Treenode> node);
}

#endif
