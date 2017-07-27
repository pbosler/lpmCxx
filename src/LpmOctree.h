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
#include <iostream>

namespace Lpm {

struct Node {        
    Node(const Box3d& bbox, Node* pparent = NULL, 
         const std::vector<index_type>& crdInds = std::vector<index_type>());
    
    inline index_type nCoords() const {return coordsContained.size();}
    bool hasKids() const {return (!kids.empty());}

    void writePoints(std::ofstream& os) const;
    
    inline void setLogProc(const int prank) {log->setProcRank(prank);}
    
    std::string infoString() const;
    
    Box3d box;
    int level;
    /// ptr to parent
    Node* parent;
    /// ptrs to child boxes
    std::vector<std::unique_ptr<Node>> kids;
    std::vector<index_type> coordsContained;
    
    void writePoints(std::ostream& os) const;
    
    static std::unique_ptr<Logger> log;
};

class Tree {
    public:
        Tree(const std::shared_ptr<Coords> crds, const scalar_type maxAspectRatio = 1.0, const int prank = 0);
        
        inline void setLogProc(const int prank) {log->setProcRank(prank); _root->log->setProcRank(prank);}
        
        std::string infoString() const;
        
        inline int depth() {return _depth;}

        inline index_type nNodes() {return _nnodes;}

        virtual void buildTree(const index_type maxCoordsPerNode);
    
        void writeToVtk(const std::string& filename, const std::string& desc = "") const;

    protected:

        void shrinkBox(Node* node);
        
        virtual void generateTree(Node* node, const index_type maxCoordsPerNode);
        int computeTreeDepth(Node* node) const;
        
        void writeVtkPoints(std::ostream& os, Node* node) const;
        void writeVtkCells(std::ostream& os, Node* node, index_type& vertIndex) const;
        void writeVtkCellType(std::ostream& os, Node* node) const;
        void writeLevelDataToVtk(std::ostream& os, Node* node) const;

        int _depth;
        index_type _nnodes;
        std::weak_ptr<Coords> _crds;    
        
        scalar_type _maxAspectRatio;
        static std::unique_ptr<Logger> log;

        std::unique_ptr<Node> _root;
};


}

#endif
