#include "LpmTreecode.h"
#include "LpmOutputMessage.h"
#include <limits>
#include <iostream>
#include <fstream>
#include <sstream>

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

Treenode::Treenode(const Box3d& bbox, const index_type nCrds, const scalar_type maxRatio) : box(bbox), 
    coordsContained(nCrds, -1), maxAspectRatio(maxRatio), level(0) {
    for (index_type i = 0; i < nCrds; ++i)
        coordsContained.push_back(i);
}

Treenode::Treenode(const Box3d& bbox, const std::shared_ptr<Treenode>& pparent, const std::vector<index_type>& crdInds, const scalar_type maxRatio) :
    box(bbox), parent(pparent), coordsContained(crdInds), maxAspectRatio(maxRatio), level(pparent->level +1) {};

void generateTree(std::shared_ptr<Treenode>& node, std::shared_ptr<Coords>& crds, const index_type maxCoordsPerNode) {
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
        const scalar_type edgeThreshold = node->box.longestEdge() / node->maxAspectRatio;
        for (int i = 0; i < 3; ++i) {
            if ( node->box.edgeLength(i) >= edgeThreshold ) {
                splitDims[i] = true;
            }   
            else {
                splitDims[i] = false;
            }
        }
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
                node->children.push_back(std::shared_ptr<Treenode>(new Treenode(kids[i], node, kidcoords, node->maxAspectRatio)));
            }
        }
        
        if (node->children.empty()) {
            OutputMessage errMsg("All children are empty.", OutputMessage::errorPriority, "Treenode::generateTree");
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

void writeTreeToVtk(const std::string& filename, const std::string& desc, const std::shared_ptr<Treenode>& root) {
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
    Lpm::writeVTKPoints(fs, root);
    fs << "CELLS " << nNodes << " " << 9 * nNodes << std::endl;
    index_type vertIndex = 0;
    Lpm::writeVtkCells(fs, root, vertIndex);
    fs << "CELL_TYPES " << nNodes << std::endl;
    Lpm::writeVtkCellType(fs, root);
    fs << "CELL_DATA " << nNodes << std::endl;
    fs << "SCALARS tree_level int 1" << std::endl;
    writeVtkLevelData(fs, root);
}

void writeVTKPoints(std::ofstream& os, const std::shared_ptr<Treenode>& node) {
    node->writePoints(os);
    if (node->hasChildren()) {
        for (int i = 0; i < node->children.size(); ++i) {
            writeVTKPoints(os, node->children[i]);
        }
    }
}

void writeVtkCells(std::ofstream& os, const std::shared_ptr<Treenode>& node, index_type& vertIndex) {
    os << 8 << " ";
    for (int i = 0; i < 8; ++i) {
        os << vertIndex + i << " ";
    }
    os << std::endl;
    vertIndex += 8;
    if (node->hasChildren()) {
        for (int i = 0; i < node->children.size(); ++i) {
            Lpm::writeVtkCells(os, node->children[i], vertIndex);
        }
    }
}

void writeVtkCellType(std::ofstream& os, const std::shared_ptr<Treenode>& node) {
    os << "VTK_HEXAHEDRON" << std::endl;
    if (node->hasChildren()) {
        for (int i = 0; i < node->children.size(); ++i) {
            writeVtkCellType(os, node->children[i]);
        }
    }
}

void writeVtkLevelData(std::ofstream& os, const std::shared_ptr<Treenode>& node) {
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
    os << box.xmin << " " << box.ymax << " " << box.zmin << std::endl;
    os << box.xmax << " " << box.ymax << " " << box.zmin << std::endl;
    os << box.xmin << " " << box.ymin << " " << box.zmax << std::endl;
    os << box.xmax << " " << box.ymin << " " << box.zmax << std::endl;
    os << box.xmin << " " << box.ymax << " " << box.zmax << std::endl;
    os << box.xmax << " " << box.ymax << " " << box.zmax << std::endl;
}

index_type nTreenodes(const std::shared_ptr<Treenode>& node) {
    index_type result = 0;
    if (node->hasChildren()) {
        for (int i = 0; i < node->children.size(); ++i) {
            result += nTreenodes(node->children[i]);
        }
    }
    else {
        result +=1;
    }
    return result;
}


void Treenode::shrinkBox(const std::shared_ptr<Coords>& crds) {
    scalar_type xmin = std::numeric_limits<scalar_type>::lowest();
    scalar_type xmax = std::numeric_limits<scalar_type>::max();
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