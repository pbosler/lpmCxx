#ifndef _LPM_TREECODE_H_
#define _LPM_TREECODE_H_

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmXyzVector.h"
#include "LpmCoords.h"
#include "LpmSphericalCoords.h"
#include "LpmLogger.h"
#include "LpmBox3d.h"
#include "LpmFaces.h"
#include <memory>
#include <vector>
#include <string>
#include <iostream>

namespace Lpm {



struct Node {        
    Node(const Box3d& bbox, Node* pparent = NULL, 
         const std::vector<index_type>& crdInds = std::vector<index_type>());
    virtual ~Node() {};
    
    inline index_type nCoords() const {return coordsContained.size();}
    bool hasKids() const {return !kids.empty();}
    bool isLeaf() const {return kids.empty();}

    void writePoints(std::ofstream& os) const;
    
    //inline void setLogProc(const int prank) {log->setProcRank(prank);}
    
    virtual std::string infoString() const;
    
    Box3d box;
    int level;
    /// ptr to parent
    Node* parent;
    /// ptrs to child boxes
    std::vector<std::unique_ptr<Node>> kids;
    std::vector<index_type> coordsContained;
    
    void writePoints(std::ostream& os) const;
    
    //static std::unique_ptr<Logger> log;
};

class Tree {
    public:
        enum TREE_DEPTH_CONTROL {MAX_COORDS_PER_NODE, MAX_DEPTH};
    
        Tree(const std::shared_ptr<Coords> crds, const scalar_type maxAspectRatio = 1.0, const int prank = 0);
        
        
        virtual ~Tree() {};
        
        //inline void setLogProc(const int prank) {log->setProcRank(prank);}
        
        std::string infoString() const;
        
        inline int depth() {return _depth;}

        inline index_type nNodes() {return _nnodes;}
        
        index_type recursiveNNodes() const {return recursiveNodeCount(_root.get());}

        virtual void buildTree(const TREE_DEPTH_CONTROL depth_type, const index_type intParam, const bool do_shrink=false);
    
        virtual void writeToVtk(const std::string& filename, const std::string& desc = "") const;

    protected:
        Tree() : _depth(0), _nnodes(0), _maxAspectRatio(1.0) {};

        void shrinkBox(Node* node);
        void shrinkBox(Node* node, std::shared_ptr<Faces> faces);
        
        void generateTreeMaxNodes(Node* node, Coords* crds, const index_type maxCoordsPerNode, const bool do_shrink);
        void generateTreeMaxDepth(Node* node, Coords* crds, const index_type maxDepth, const bool do_shrink);
        
        int computeTreeDepth(Node* node) const;
        
        index_type recursiveNodeCount(Node* node) const;
        
        void writeVtkPoints(std::ostream& os, Node* node) const;
        void writeVtkCells(std::ostream& os, Node* node, index_type& vertIndex) const;
        void writeVtkCellType(std::ostream& os, Node* node) const;
        void writeLevelDataToVtk(std::ostream& os, Node* node) const;

        int _depth;
        index_type _nnodes;
        std::weak_ptr<Coords> _crds;    
        
        scalar_type _maxAspectRatio;
//         static std::unique_ptr<Logger> log;

        std::unique_ptr<Node> _root;
};


}

#endif
