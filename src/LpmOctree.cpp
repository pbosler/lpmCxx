#include "LpmOctree.h"
#include "LpmOutputMessage.h"
#include <limits>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <cmath>

namespace Lpm {

std::unique_ptr<Logger> Tree::log(new Logger(OutputMessage::debugPriority, "Tree_log"));
//std::unique_ptr<Logger> Node::log(new Logger(OutputMessage::debugPriority, "Node_log"));

Tree::Tree(const std::shared_ptr<Coords> crds, const scalar_type maxRatio, const int prank) : _crds(crds), _maxAspectRatio(maxRatio),
    _depth(0), _nnodes(1) {
    log->setProcRank(prank);
    
    Box3d rbox(crds->minX() - BOX_PADDING_FACTOR * crds->minX(), 
      crds->maxX() + BOX_PADDING_FACTOR * crds->maxX(),
      crds->minY() - BOX_PADDING_FACTOR * crds->minY(), 
      crds->maxY() + BOX_PADDING_FACTOR * crds->maxY(), 
      crds->minZ() - BOX_PADDING_FACTOR * crds->minZ(), 
      crds->maxZ() + BOX_PADDING_FACTOR * crds->maxZ());
            
    std::vector<index_type> rinds;
    rinds.reserve(crds->n());
    for (index_type i = 0; i < crds->n(); ++i) {
        rinds.push_back(i);
    }
    
    Node* rparent = NULL;
    
    _root = std::unique_ptr<Node>(new Node(rbox, rparent, rinds));
}

std::string Tree::infoString() const {
    std::stringstream ss;
    ss << "Tree info:" << std::endl;
    ss << "\tnNodes = " << _nnodes << std::endl;
    ss << "\tdepth = " << _depth << std::endl;
    ss << "\tmaxAspectRatio = " << _maxAspectRatio << std::endl;
    //ss << _root->infoString();
    return ss.str();
}

std::string Node::infoString() const {
    std::stringstream ss;
    ss << "Treenode info: " << std::endl;
    ss << "\tbox: " << box.infoString();
    ss << "\tlevel = " << level << std::endl;
    ss << "\t&parent = " << (level > 0 ? parent : 0) << std::endl;
    ss << "\thasKids = " << (hasKids() ? "true" : "false") << ", kids.size() = " << kids.size() << std::endl;
    ss << "\tcoordsContained.size() = " << coordsContained.size() << std::endl;
    if (!coordsContained.empty()) {
        ss << "\tmin coord ind = " << *std::min_element(coordsContained.begin(), coordsContained.end())
            << ", max coord in = " << *std::max_element(coordsContained.begin(), coordsContained.end()) << std::endl;
    }
    return ss.str();
}

Node::Node(const Box3d& bbox, Node* pparent, const std::vector<index_type>& crdInds) :
    box(bbox), parent(pparent), coordsContained(crdInds) {
    if (pparent) {
        level = pparent->level + 1;
    }
    else {
        level = 0;
    }
}

void Tree::buildTree(const TREE_DEPTH_CONTROL tree_depth_type, const index_type intParam, const bool do_shrink) {
    std::shared_ptr<Coords> crd_ptr = _crds.lock();
    if (tree_depth_type == MAX_COORDS_PER_NODE) {
        generateTreeMaxNodes(_root.get(), crd_ptr.get(), intParam, do_shrink);
    }
    else if (tree_depth_type == MAX_DEPTH) {
        generateTreeMaxDepth(_root.get(), crd_ptr.get(), intParam, do_shrink);
    }
    _depth = computeTreeDepth(_root.get());
}

void Tree::generateTreeMaxNodes(Node* node, Coords* crd_ptr, const index_type maxCoordsPerNode, const bool do_shrink) {
    if (node->coordsContained.size() <= maxCoordsPerNode) {
        return;
    }
    else {
        //
        //  determine box dimensions to split
        //
//         bool splitDims[3];
//         int splitCount = 0;
//         const scalar_type edgeThreshold = node->box.longestEdge() / _maxAspectRatio;
//         for (int i = 0; i < 3; ++i) {
//             if (node->box.edgeLength(i) >= edgeThreshold) {
//                 splitDims[i] = true;
//                 splitCount += 1;
//             }
//             else {
//                 splitDims[i] = false;
//             }
//         }
        //
        //  make child boxes
        //
//        std::vector<Box3d> kidboxes = node->box.bisectAlongDims(splitDims);


        std::vector<Box3d> kidboxes = node->box.bisectAll();
        for (int i = 0; i < kidboxes.size(); ++i) {
            //
            //  find coordinates contained by child box
            //
            std::vector<index_type> kidcoords;
            kidcoords.reserve(node->coordsContained.size());
            for (index_type j = 0; j < node->coordsContained.size(); ++j) {
                if (kidboxes[i].containsPoint(crd_ptr->getVec(node->coordsContained[j]))) {
                    kidcoords.push_back(node->coordsContained[j]);
                }
            }
            if (!kidcoords.empty()) {
                kidcoords.shrink_to_fit();
                node->kids.push_back(std::unique_ptr<Node>(new Node(kidboxes[i], node, kidcoords)));
            }
        }
        if (node->kids.empty()) {
            OutputMessage errMsg("All kids are empty, this shouldn't happen.", OutputMessage::errorPriority, "Tree::generateTreeMaxNodes");
            log->logMessage(errMsg);
            return;
        }
        else {
            _nnodes += node->kids.size();
             for (int i = 0; i < node->kids.size(); ++i) {
                if (do_shrink) {
                    shrinkBox(node->kids[i].get());  
                }
                generateTreeMaxNodes(node->kids[i].get(), crd_ptr, maxCoordsPerNode, do_shrink);
            }
        }
    }
}

void Tree::generateTreeMaxDepth(Node* node, Coords* crd_ptr, const index_type maxDepth, const bool do_shrink) {
    if (node->level == maxDepth) {
        return;
    }
    else {
        std::vector<Box3d> kidboxes = node->box.bisectAll();
        for (int i = 0; i < kidboxes.size(); ++i) {
            std::vector<index_type> kidcoords;
            kidcoords.reserve(node->coordsContained.size());
            for (index_type j = 0; j < node->coordsContained.size(); ++j) {
                if (kidboxes[i].containsPoint(crd_ptr->getVec(node->coordsContained[j]))) {
                    kidcoords.push_back(node->coordsContained[j]);
                }
            }
            if (!kidcoords.empty()) {
                kidcoords.shrink_to_fit();
                node->kids.push_back(std::unique_ptr<Node>(new Node(kidboxes[i], node, kidcoords)));
            }
        }
        if (node->kids.empty()) {
            OutputMessage errMsg("All kids are empty, this shouldn't happen.", OutputMessage::errorPriority, "Tree::generateTreeMaxDepth");
            log->logMessage(errMsg);
            return;
        }
        else {
            _nnodes += node->kids.size();
            for (int i=0; i<node->kids.size(); ++i) {
                if (do_shrink) {
                    shrinkBox(node->kids[i].get());
                }
                generateTreeMaxDepth(node->kids[i].get(), crd_ptr, maxDepth, do_shrink);
            }
        }
    }
}


void Tree::shrinkBox(Node* node) {
    scalar_type xmin = std::numeric_limits<scalar_type>::max();
    scalar_type xmax = std::numeric_limits<scalar_type>::lowest();
    scalar_type ymin = xmin;
    scalar_type ymax = xmax;
    scalar_type zmin = xmin;
    scalar_type zmax = xmax;
    std::shared_ptr<Coords> crd_ptr = _crds.lock();
    for (index_type i = 0; i < node->coordsContained.size(); ++i) {
        const XyzVector crdVec = crd_ptr->getVec(node->coordsContained[i]);
        if (crdVec.x < xmin)
            xmin = crdVec.x;
        if (crdVec.x > xmax)
            xmax = crdVec.x;
        if (crdVec.y < ymin)
            ymin = crdVec.y;
        if (crdVec.y > ymax)
            ymax = crdVec.y;
        if (crdVec.z < zmin)
            zmin = crdVec.z;
        if (crdVec.z > zmax)
            zmax = crdVec.z;
    }
    node->box.xmin = xmin;
    node->box.xmax = xmax;
    node->box.ymin = ymin;
    node->box.ymax = ymax;
    node->box.zmin = zmin;
    node->box.zmax = zmax;
}

void Tree::shrinkBox(Node* node, std::shared_ptr<Faces> faces) {
    scalar_type xmin = std::numeric_limits<scalar_type>::max();
    scalar_type xmax = std::numeric_limits<scalar_type>::lowest();
    scalar_type ymin = xmin;
    scalar_type ymax = xmax;
    scalar_type zmin = xmin;
    scalar_type zmax = xmax;
    for (index_type i = 0; i < node->coordsContained.size(); ++i) {
        const XyzVector cntd = faces->centroid(node->coordsContained[i]);
        if (cntd.x < xmin)
            xmin = cntd.x;
        if (cntd.x > xmax)
            xmax = cntd.x;
        if (cntd.y < ymin)
            ymin = cntd.y;
        if (cntd.y > ymax)
            ymax = cntd.y;
        if (cntd.z < zmin)
            zmin = cntd.z;
        if (cntd.z > zmax)
            zmax = cntd.z;
    }
    node->box.xmin = xmin;
    node->box.xmax = xmax;
    node->box.ymin = ymin;
    node->box.ymax = ymax;
    node->box.zmin = zmin;
    node->box.zmax = zmax;
}


void Tree::writeToVtk(const std::string& filename, const std::string& desc) const {
    std::ofstream fs(filename);
    if (!fs.is_open()) {
        OutputMessage errMsg("cannot open .vtk file", OutputMessage::errorPriority, "Lpm::Tree::writeToVtk");
        log->logMessage(errMsg);
        throw std::ios_base::failure("file write error");    
    }
    
    fs << "# vtk DataFile Version 2.0" << std::endl;
    fs << desc << std::endl;
    fs << "ASCII" << std::endl;
    fs << "DATASET UNSTRUCTURED_GRID" << std::endl;
    fs << "POINTS " << 8 * _nnodes << " double" << std::endl;
    writeVtkPoints(fs, _root.get());
    fs << "CELLS " << _nnodes << " " << 9 * _nnodes << std::endl;
    index_type vertIndex = 0;
    writeVtkCells(fs, _root.get(), vertIndex);
    fs << "CELL_TYPES " << _nnodes << std::endl;
    writeVtkCellType(fs, _root.get());
    fs << "CELL_DATA " << _nnodes << std::endl;
    fs << "SCALARS tree_level int 1" << std::endl;
    fs << "LOOKUP_TABLE default" << std::endl;
    writeLevelDataToVtk(fs, _root.get());
}

void Tree::writeVtkPoints(std::ostream& os, Node* node) const {
    node->writePoints(os);
    if (node->hasKids()) {
        for (int i = 0; i < node->kids.size(); ++i) {
            writeVtkPoints(os, node->kids[i].get());
        }
    }
}

void Node::writePoints(std::ostream& os) const {
    os << box.xmin << " " << box.ymin << " " << box.zmin << std::endl;
    os << box.xmax << " " << box.ymin << " " << box.zmin << std::endl;
    os << box.xmax << " " << box.ymax << " " << box.zmin << std::endl;
    os << box.xmin << " " << box.ymax << " " << box.zmin << std::endl;
    os << box.xmin << " " << box.ymin << " " << box.zmax << std::endl;
    os << box.xmax << " " << box.ymin << " " << box.zmax << std::endl;
    os << box.xmax << " " << box.ymax << " " << box.zmax << std::endl;
    os << box.xmin << " " << box.ymax << " " << box.zmax << std::endl;
}


void Tree::writeVtkCells(std::ostream& os, Node* node, index_type& vertIndex) const {
    os << 8 << " ";
    for (int i = 0; i < 8; ++i) {
        os << vertIndex + i << " ";
    }
    os << std::endl;
    vertIndex += 8;
    if (node->hasKids()) {
        for (int i = 0; i < node->kids.size(); ++i) {
             writeVtkCells(os, node->kids[i].get(), vertIndex);
        }
    }
}

void Tree::writeVtkCellType(std::ostream& os, Node* node) const {
    os << "12" << std::endl;
    if (node->hasKids()) {
        for (int i = 0; i < node->kids.size(); ++i) {
            writeVtkCellType(os, node->kids[i].get());
        }
    }
}

void Tree::writeLevelDataToVtk(std::ostream& os, Node* node) const {
    os << node->level << std::endl;
    if (node->hasKids()) {
        for (int i = 0; i < node->kids.size(); ++i) {
            writeLevelDataToVtk(os, node->kids[i].get());
        }
    }
}

index_type Tree::recursiveNodeCount(Node* node) const {
    index_type result = 1;
    result += node->kids.size();
    if (node->hasKids()) {
        for (int i = 0; i < node->kids.size(); ++i) {
            result += recursiveNodeCount(node->kids[i].get());
        }
    }
    return result;
}

int Tree::computeTreeDepth(Node* node) const {
    int result = 0;
    if (node->level > result)
        result = node->level;
    if (node->hasKids()) {
        for (int i = 0; i < node->kids.size(); ++i) {
            result = computeTreeDepth(node->kids[i].get());
        }
    }
    return result;
}

}