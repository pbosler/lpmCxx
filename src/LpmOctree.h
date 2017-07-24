#ifndef _LPM_TREECODE_H_
#define _LPM_TREECODE_H_

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmXyzVector.h"
#include "LpmCoords.h"
#include "LpmLogger.h"
#include "LpmBox3d.h"
#include <memory>
#include <vector>
#include <string>

namespace Lpm {

class Tree {
    Tree(const Coords* crds, const scalar_type maxAspectRatio = 1.0);

    struct Node {
        Node();
        
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
