#include "LpmOctree.h"
#include "LpmOutputMessage.h"
#include <limits>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <numeric>

namespace Lpm {

std::unique_ptr<Logger> Treenode::log(new Logger(OutputMessage::debugPriority, "Treenode_log"));

scalar_type Box3d::longestEdge() const {
    scalar_type result = xmax - xmin;
    if (ymax-ymin > result)
        result = ymax - ymin;
    if (zmax-zmin > result)
        result = zmax - zmin;
    return result;
}

scalar_type Box3d::shortestEdge() const {
    const scalar_type dx = xmax - xmin;
    const scalar_type dy = ymax - ymin;
    const scalar_type dz = zmax - zmin;
    scalar_type result = dx;
    if (dy < dx)
        result = dy;
    if (dz > 0.0 && dz < result)
        result = dz;
    return result;
}

scalar_type Box3d::aspectRatio() const {
    return longestEdge() / shortestEdge();
}

scalar_type Box3d::radius() const {
    const XyzVector cntd = centroid();
    std::vector<XyzVector> corners;
    corners.push_back(XyzVector(xmin, ymin, zmin));
    corners.push_back(XyzVector(xmin, ymax, zmin));
    corners.push_back(XyzVector(xmin, ymin, zmax));
    corners.push_back(XyzVector(xmin, ymax, zmax));
    corners.push_back(XyzVector(xmax, ymin, zmin));
    corners.push_back(XyzVector(xmax, ymax, zmin));
    corners.push_back(XyzVector(xmax, ymin, zmax));
    corners.push_back(XyzVector(xmax, ymax, zmax));
    scalar_type result = 0.0;
    for (int i = 0; i < 8; ++i) {
        const scalar_type testDist = distance(cntd, corners[i]);
        if (testDist > result)
            result = testDist;
    }
    return result;
}

std::vector<Box3d> Box3d::bisectAll() const {
    std::vector<Box3d> result;
    result.push_back(Box3d(xmin, 0.5 * (xmin + xmax), ymin, 0.5 * (ymin + ymax), zmin, 0.5 * (zmin + zmax)));
    result.push_back(Box3d(0.5 * (xmin + xmax), xmax, ymin, 0.5 * (ymin + ymax), zmin, 0.5 * (zmin + zmax)));
    result.push_back(Box3d(xmin, 0.5 * (xmin + xmax), 0.5 * (ymin + ymax), ymax, zmin, 0.5 * (zmin + zmax)));
    result.push_back(Box3d(0.5 * (xmin + xmax), xmax, 0.5 * (ymin + ymax), ymax, zmin, 0.5 * (zmin + zmax)));
    result.push_back(Box3d(xmin, 0.5 * (xmin + xmax), ymin, 0.5 * (ymin + ymax), 0.5 * (zmin + zmax), zmax));
    result.push_back(Box3d(0.5 * (xmin + xmax), xmax, ymin, 0.5 * (ymin + ymax), 0.5 * (zmin + zmax), zmax));
    result.push_back(Box3d(xmin, 0.5 * (xmin + xmax), 0.5 * (ymin + ymax), ymax, 0.5 * (zmin + zmax), zmax));
    result.push_back(Box3d(0.5 * (xmin + xmax), xmax, 0.5 * (ymin + ymax), ymax, 0.5 * (zmin + zmax), zmax));
    return result;
}

std::vector<Box3d> Box3d::bisectAlongDims(const bool* dims) const {
    std::vector<Box3d> result;
    int dimCount = 0;
    for (int i = 0; i < 3; ++i) {
        if (dims[i])
            dimCount += 1;
    }
    if (dimCount == 1) {
        if (dims[0]) {
            result.push_back(Box3d(xmin, 0.5 * (xmin + xmax), ymin, ymax, zmin, zmax));
            result.push_back(Box3d(0.5 * (xmin + xmax), xmax, ymin, ymax, zmin, zmax));
        }
        else if (dims[1]) {
            result.push_back(Box3d(xmin, xmax, ymin, 0.5 * (ymin + ymax), zmin, zmax));
            result.push_back(Box3d(xmin, xmax, 0.5 * (ymin + ymax), ymax, zmin, zmax));
        }
        else if (dims[2]) {
            result.push_back(Box3d(xmin, xmax, ymin, ymax, zmin, 0.5 * (zmin + zmax)));
            result.push_back(Box3d(xmin, xmax, ymin, ymax, 0.5 * (zmin + zmax), zmax));
        }
    }
    else if (dimCount == 2) {
        if (dims[0] && dims[1]) {
            result.push_back(Box3d(xmin, 0.5 * (xmin + xmax), ymin, 0.5 * (ymin + ymax), zmin, zmax));
            result.push_back(Box3d(0.5 * (xmin + xmax), xmax, ymin, 0.5 * (ymin + ymax), zmin, zmax));
            result.push_back(Box3d(xmin, 0.5 * (xmin + xmax), 0.5 * (ymin + ymax), ymax, zmin, zmax));
            result.push_back(Box3d(0.5 * (xmin + xmax), xmax, 0.5 * (ymin + ymax), ymax, zmin, zmax));
        }
        else if (dims[0] && dims[2]) {
            result.push_back(Box3d(xmin, 0.5 * (xmin + xmax), ymin, ymax, zmin, 0.5 * (zmin + zmax)));
            result.push_back(Box3d(0.5 * (xmin + xmax), xmax, ymin, ymax, zmin, 0.5 * (zmin + zmax)));
            result.push_back(Box3d(xmin, 0.5 * (xmin + xmax), ymin, ymax, 0.5 * (zmin + zmax), zmax));
            result.push_back(Box3d(0.5 * (xmin + xmax), xmax, ymin, ymax, 0.5 * (zmin + zmax), zmax));
        }
        else if (dims[1] && dims[2]) {
            result.push_back(Box3d(xmin, xmax, ymin, 0.5 * (ymin + ymax), zmin, 0.5 * (zmin + zmax)));
            result.push_back(Box3d(xmin, xmax, 0.5 * (ymin + ymax), ymax, zmin, 0.5 * (zmin + zmax)));
            result.push_back(Box3d(xmin, xmax, ymin, 0.5 * (ymin + ymax), 0.5 * (zmin + zmax), zmax));
            result.push_back(Box3d(xmin, xmax, 0.5 * (ymin + ymax), ymax, 0.5 * (zmin + zmax), zmax));
        }
    }
    else {
        result = bisectAll();
    }
    return result;
}

scalar_type Box3d::edgeLength(const int dim) const {
    scalar_type result = 0.0;
    switch (dim) {
        case 0 : {
            result = xmax - xmin;
            break;
        }
        case 1 : {
            result = ymax - ymin;
            break;
        }
        case 2 : {
            result = zmax - zmin;
            break;
        }
    }
    return result;
}

std::string Box3d::infoString() const {
    std::stringstream ss;
    ss << "x in [" << xmin << ", " << xmax << "], y in [" << ymin << ", " << ymax << "], z in [" << zmin << ", " << zmax << "]" << std::endl;
    return ss.str();
}

Treenode::Treenode() : box(Box3d()), maxAspectRatio(1.0), level(-1), coordsContained(std::vector<index_type>()), 
    parent(NULL), children(std::vector<std::shared_ptr<Treenode>>()) {}

Treenode::Treenode(const std::shared_ptr<Coords> crds, const scalar_type maxRatio) : 
    box(crds->minX(), crds->maxX(), crds->minY(), crds->maxY(), crds->minZ(), crds->maxZ()),
    coordsContained(crds->n(), -1), maxAspectRatio(maxRatio), level(0), 
    parent(NULL), children(std::vector<std::shared_ptr<Treenode>>()) {
    for (index_type i = 0; i < crds->n(); ++i) {
        coordsContained[i] = i;
    }
}

std::string Treenode::infoString() const {
    std::stringstream ss;
    ss << "Treenode info: " << std::endl;
    ss << "\tbox: " << box.infoString();
    ss << "\tmaxAspectRatio = " << maxAspectRatio << std::endl;
    ss << "\tlevel = " << level << std::endl;
    ss << "\tparent = " << (parent.use_count() > 0 ? parent.get() : 0) << std::endl;
    ss << "\thasChildren = " << (hasChildren() ? "true" : "false") << ", children.size() = " << children.size() << std::endl;
    ss << "\tcoordsContained.size() = " << coordsContained.size() << std::endl;
    if (!coordsContained.empty()) {
        ss << "\tmin coord ind = " << *std::min_element(coordsContained.begin(), coordsContained.end())
            << ", max coord in = " << *std::max_element(coordsContained.begin(), coordsContained.end()) << std::endl;
    }
    return ss.str();
}

Treenode::Treenode(const Box3d& bbox, const std::shared_ptr<Treenode> pparent, const std::vector<index_type>& crdInds, const scalar_type maxRatio) :
    box(bbox), parent(pparent), coordsContained(crdInds), maxAspectRatio(maxRatio), level(pparent->level +1) {};

void generateTree(std::shared_ptr<Treenode> node, std::shared_ptr<Coords> crds, const index_type maxCoordsPerNode) {
    //
    //  if node has less than maximum allowed points, we're done
    //
//     OutputMessage trcMsg("Entering generateTree", OutputMessage::tracePriority, "LpmOctree::generateTree");
//     node->log->logMessage(trcMsg);
    if (node->coordsContained.size() <= maxCoordsPerNode) {
//         std::cout << "box contains " << node->coordsContained.size() << ", at level " << node->level << "; returning." << std::endl;
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
//         std::cout << "splitting box along " << dimcount << " dimensions" << std::endl;
        //
        //  define child boxes
        //
        std::vector<Box3d> kids = node->box.bisectAlongDims(splitDims);
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
//                 std::cout << "adding child " << i << " : " << kids[i].infoString() << std::endl;
//                 std::cout << "\tcoordsContained: ";
//                 for (index_type j = 0; j < kidcoords.size(); ++j) 
//                     std::cout << kidcoords[j] << " ";
//                 std::cout << std::endl;
                node->children.push_back(std::shared_ptr<Treenode>(new Treenode(kids[i], node, kidcoords, node->maxAspectRatio)));
            }
        }
        
        if (node->children.empty()) {
            OutputMessage errMsg("All children are empty, this shouldn't happen.", OutputMessage::errorPriority, "Treenode::generateTree");
            node->log->logMessage(errMsg);
            return;
        }
        else {
            for (int i = 0; i < node->children.size(); ++i) {
                node->children[i]->shrinkBox(crds);
                generateTree(node->children[i], crds, maxCoordsPerNode);
            }
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
    if (node->hasChildren()) {
        for (int i = 0; i < node->children.size(); ++i) {
            writeVTKPoints(os, node->children[i]);
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
    if (node->hasChildren()) {
        for (int i = 0; i < node->children.size(); ++i) {
             writeVtkCells(os, node->children[i], vertIndex);
        }
    }
}

void writeVtkCellType(std::ofstream& os, const std::shared_ptr<Treenode> node) {
    os << "12" << std::endl;
    if (node->hasChildren()) {
        for (int i = 0; i < node->children.size(); ++i) {
            writeVtkCellType(os, node->children[i]);
        }
    }
}

void writeVtkLevelData(std::ofstream& os, const std::shared_ptr<Treenode> node) {
    os << node->level << std::endl;
    if (node->hasChildren()) {
        for (int i = 0; i < node->children.size(); ++i) {
            writeVtkLevelData(os, node->children[i]);
        }
    }
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
    if (node->hasChildren()) {
        for (int i = 0; i < node->children.size(); ++i) {
            result += nTreenodes(node->children[i]);
        }
    }
    return result;
}


void Treenode::shrinkBox(const std::shared_ptr<Coords>& crds) {
    scalar_type xmin = std::numeric_limits<scalar_type>::max();
    scalar_type xmax = std::numeric_limits<scalar_type>::lowest();
    scalar_type ymin = xmin;
    scalar_type ymax = xmax;
    scalar_type zmin = xmin;
    scalar_type zmax = xmax;
    for (index_type i = 0; i < coordsContained.size(); ++i) {
        const XyzVector crdVec = crds->getVec(coordsContained[i]);
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
    box.xmin = xmin;
    box.xmax = xmax;
    box.ymin = ymin;
    box.ymax = ymax;
    box.zmin = zmin;
    box.zmax = zmax;
}

}