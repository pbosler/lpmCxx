#include "LpmOctree.h"
#include "LpmOutputMessage.h"
#include <limits>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <cmath>

#define BOX_PADDING_FACTOR 0.00001

namespace Lpm {

std::unique_ptr<Logger> Tree::log(new Logger(OutputMessage::debugPriority, "Tree_log"));

Tree::Tree(const Coords* crds, const scalar_type maxRatio, const int prank) : _crds(crds), _maxAspectRatio(maxRatio),
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
    
    root = std::unique_ptr<Node>(new Node(rbox, rparent, rinds));
}

std::string Tree::infoString() const {
    std::stringstream ss;
    ss << "Tree info:" << std::endl;
    ss << "\tnNodes = " << _nnodes << std::endl;
    ss << "\tdepth = " << _depth << std::endl;
    ss << "\tmaxAspectRatio = " << _maxAspectRatio << std::endl;
    return ss.str();
}

std::string Tree::Node::infoString() const {
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

Tree::Node::Node(const Box3d& bbox, const Node* pparent, const std::vector<index_type>& crdInds) :
    box(bbox), parent(pparent), coordsContained(crdInds) {
    if (pparent) {
        level = pparent->level + 1;
    }
    else {
        level = 0;
    }
}

void generateTree(std::shared_ptr<Treenode> node, std::shared_ptr<Coords> crds, const index_type maxCoordsPerNode) {
    //
    //  if node has less than maximum allowed points, we're done
    //
    if (node->coordsContained.size() <= maxCoordsPerNode) {
        return;
    }
    else { 
        // 
        //   divide node
        // 
        bool splitDims[3];
        int dimcount = 0;
        const scalar_type edgeThreshold = node->box.longestEdge() / node->maxAspectRatio;
        for (int i = 0; i < 3; ++i) {
            if ( node->box.edgeLength(i) >= edgeThreshold ) {
                splitDims[i] = true;
                dimcount += 1;
            }   
            else {
                splitDims[i] = false;
            }
        }
        //
        //  define child boxes
        //
        std::vector<Box3d> kids = node->box.bisectAlongDims(splitDims);
        //std::vector<Box3d> kids = node->box.bisectAll();
        for (int i = 0; i < kids.size(); ++i) {
            std::vector<index_type> kidcoords;
            //
            //  find coordinates contained by each child
            //
            for (index_type j = 0; j < node->coordsContained.size(); ++j) {
                if (kids[i].containsPoint(crds->getVec(node->coordsContained[j]))) {
                    kidcoords.push_back(node->coordsContained[j]);
                }
            }
            //
            //  only add non-empty children
            //
            if (!kidcoords.empty()) {
                node->kids.push_back(std::shared_ptr<Treenode>(new Treenode(kids[i], node, kidcoords)));
            }
        }
        
        if (node->kids.empty()) {
            OutputMessage errMsg("All kids are empty, this shouldn't happen.", OutputMessage::errorPriority, "Treenode::generateTree");
            node->log->logMessage(errMsg);
            return;
        }
        else {
            for (int i = 0; i < node->kids.size(); ++i) {
                node->kids[i]->shrinkBox(crds);
                generateTree(node->kids[i], crds, maxCoordsPerNode);
            }
        }
    }
}

void Tree::buildTree(std::unique_ptr<Tree::Node> node, const index_type maxCoordsPerNode) {
    if (node->coordsContained.size() <= maxCoordsPerNode) {
        return;
    }
    else {
        //
        //  determine box dimensions to split
        //
        bool splitDims[3];
        int splitCount = 0;
        const scalar_type edgeThreshold = node->box.longestEdge() / _maxAspectRatio;
        for (int i = 0; i < 3; ++i) {
            if (node->box.edgeLength(i) >= edgeThreshold) {
                splitDims[i] = true;
                splitCount += 1;
            }
            else {
                splitDims[i] = false;
            }
        }
        //
        //  make child boxes
        //
        std::vector<Box3d> kidboxes = node->box.bisectAlongDims(splitDims);
        for (int i = 0; i < kidboxes.size(); ++i) {
            //
            //  find coordinates contained by child box
            //
            std::vector<index_type> kidcoords;
            kidcoords.reserve(node->coordsContained.size());
            for (index_type j = 0; j < node->coordsContained.size(); ++j) {
                if (kidboxes[i].containsPoint(_crds->getVec(node->coordsContained[j]))) {
                    kidcoords.push_back(node->coordsContained[j]);
                }
            }
            kidcoords.shrink_to_fit();
            if (!kidcoords.empty()) {
                node->kids.push_back(std::unique_ptr<Node>(new Node(kidboxes[i], node.get(), kidcoords)));
            }
        }
        if (node->kids.empty()) {
            OutputMessage errMsg("All kids are empty, this shouldn't happen.", OutputMessage::errorPriority, "Treenode::generateTree");
            node->log->logMessage(errMsg);
            return;
        }
        else {
            _depth += 1;
            _nnodes += node->kids.size();
        }
        for (int i = 0; i < node->kids.size(); ++i) {
             shrinkBox(node->kids[i]);  
             buildTree(node->kids[i], maxCoordsPerNode);
        }
    }
}


void writeTreeToVtk(const std::string& filename, const std::string& desc, const std::shared_ptr<Treenode> root) {
    std::ofstream fs(filename);
    if (!fs.is_open()) {
        OutputMessage errMsg("cannot open .vtk file", OutputMessage::errorPriority, "Lpm::writeTreeToVtk");
        root->log->logMessage(errMsg);
        throw std::ios_base::failure("file write error");
    }
    const index_type nNodes = nTreenodes(root); 
    fs << "# vtk DataFile Version 2.0" << std::endl;
    fs << desc << std::endl;
    fs << "ASCII" << std::endl;
    fs << "DATASET UNSTRUCTURED_GRID" << std::endl;
    fs << "POINTS " << 8 * nNodes << " double" << std::endl;
    writeVTKPoints(fs, root);
    fs << "CELLS " << nNodes << " " << 9 * nNodes << std::endl;
    index_type vertIndex = 0;
    writeVtkCells(fs, root, vertIndex);
    fs << "CELL_TYPES " << nNodes << std::endl;
    Lpm::writeVtkCellType(fs, root);
    fs << "CELL_DATA " << nNodes << std::endl;
    fs << "SCALARS tree_level int 1" << std::endl;
    fs << "LOOKUP_TABLE default" << std::endl;
    writeVtkLevelData(fs, root);
}

void writeVTKPoints(std::ofstream& os, const std::shared_ptr<Treenode> node) {
    node->writePoints(os);
    if (node->hasKids()) {
        for (int i = 0; i < node->kids.size(); ++i) {
            writeVTKPoints(os, node->kids[i]);
        }
    }
}

void writeVtkCells(std::ofstream& os, const std::shared_ptr<Treenode> node, index_type& vertIndex) {
    os << 8 << " ";
    for (int i = 0; i < 8; ++i) {
        os << vertIndex + i << " ";
    }
    os << std::endl;
    vertIndex += 8;
    if (node->hasKids()) {
        for (int i = 0; i < node->kids.size(); ++i) {
             writeVtkCells(os, node->kids[i], vertIndex);
        }
    }
}

void writeVtkCellType(std::ofstream& os, const std::shared_ptr<Treenode> node) {
    os << "12" << std::endl;
    if (node->hasKids()) {
        for (int i = 0; i < node->kids.size(); ++i) {
            writeVtkCellType(os, node->kids[i]);
        }
    }
}

void writeVtkLevelData(std::ofstream& os, const std::shared_ptr<Treenode> node) {
    os << node->level << std::endl;
    if (node->hasKids()) {
        for (int i = 0; i < node->kids.size(); ++i) {
            writeVtkLevelData(os, node->kids[i]);
        }
    }
}

int treeDepth(const std::shared_ptr<Treenode> node) {
    int result = 0; 
    if (node->level > result)
        result = node->level;
    if (node->hasKids()) {
        for (int i = 0; i < node->kids.size(); ++i) {
            result = treeDepth(node->kids[i]);
        }
    }
    return result;
}

void Treenode::writePoints(std::ofstream& os) const {
    os << box.xmin << " " << box.ymin << " " << box.zmin << std::endl;
    os << box.xmax << " " << box.ymin << " " << box.zmin << std::endl;
    os << box.xmax << " " << box.ymax << " " << box.zmin << std::endl;
    os << box.xmin << " " << box.ymax << " " << box.zmin << std::endl;
    os << box.xmin << " " << box.ymin << " " << box.zmax << std::endl;
    os << box.xmax << " " << box.ymin << " " << box.zmax << std::endl;
    os << box.xmax << " " << box.ymax << " " << box.zmax << std::endl;
    os << box.xmin << " " << box.ymax << " " << box.zmax << std::endl;
}

index_type nTreenodes(const std::shared_ptr<Treenode> node) {
    index_type result = 1;
    if (node->hasKids()) {
        for (int i = 0; i < node->kids.size(); ++i) {
            result += nTreenodes(node->kids[i]);
        }
    }
    return result;
}


void Tree::shrinkBox(std::unique_ptr<Node> node) {
    scalar_type xmin = std::numeric_limits<scalar_type>::max();
    scalar_type xmax = std::numeric_limits<scalar_type>::lowest();
    scalar_type ymin = xmin;
    scalar_type ymax = xmax;
    scalar_type zmin = xmin;
    scalar_type zmax = xmax;
    for (index_type i = 0; i < node->coordsContained.size(); ++i) {
        const XyzVector crdVec = _crds->getVec(node->coordsContained[i]);
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

}