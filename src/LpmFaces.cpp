#include "LpmFaces.h"
#include "LpmOutputMessage.h"
#include <iostream>
#include <sstream>
#include <numeric>
#include <algorithm>

namespace Lpm {

std::unique_ptr<Logger> Faces::log(new Logger(OutputMessage::debugPriority, "Faces_base_log"));

Faces::Faces(const index_type nMax, const index_type nMaxEdgesPerFace, const std::shared_ptr<Edges> edge_ptr, 
    const std::shared_ptr<Coords> crd_ptr,  const std::shared_ptr<Coords> lag_crd_ptr, const bool sim3d) : 
    _nMax(nMax), _nMaxEdges(nMaxEdgesPerFace), edges(edge_ptr), crds(crd_ptr), lagCrds(lag_crd_ptr), _is3d(sim3d) {
    _edgeInds.reserve(nMax);     
    _area.reserve(nMax);
    _hasChildren.reserve(nMax);
    _parent.reserve(nMax);
    _children.reserve(nMax);
    if (sim3d) {
        _positiveCell.reserve(nMax);
        _negativeCell.reserve(nMax);
    }
}

index_type Faces::nDivided() const {
    return std::count(_hasChildren.begin(), _hasChildren.end(), true);
}

void Faces::insert(const std::vector<index_type>& eInds) {
    _edgeInds.push_back(eInds);
    _area.push_back(0.0);
    _hasChildren.push_back(false);
    _parent.push_back(-1);
    std::vector<index_type> kids = {-1, -1, -1, -1};
    _children.push_back(kids);
    if (_is3d) {
        _positiveCell.push_back(-1);
        _negativeCell.push_back(-1);
    }
}

std::vector<index_type> Faces::ccwAdjacentFaces(const index_type ind) const {
    std::vector<index_type> result;
    std::vector<index_type> edgeList = _edgeInds[ind];
    for (index_type i = 0; i < edgeList.size(); ++i) {
        if (edgeIsPositive(ind, edgeList[i])) {
            result.push_back(edges->rightFace(edgeList[i]));
        }
        else {
            result.push_back(edges->leftFace(edgeList[i]));
        }
    }
    return result;
}

scalar_type Faces::surfaceArea() const {
    return std::accumulate(_area.begin(), _area.end(), scalar_type(0.0));
}

XyzVector Faces::centroid(const index_type i, const bool lagrangian) const {
    const std::vector<index_type> verts = vertexIndices(i);
    if (lagrangian)
        return lagCrds->centroid(verts);
    else
        return crds->centroid(verts);    
}

scalar_type Faces::computeArea(const index_type i) {
    const XyzVector cntd = centroid(i);
    const std::vector<index_type> verts = vertexIndices(i);
    scalar_type area = 0.0;
    for (index_type j = 0; j < verts.size(); ++j) {
        area += crds->triArea(cntd, verts[j], verts[(j+1)%verts.size()]);
    }
    _area[i] = area;
    return area;
}

void Faces::resetAreas() {
    scalar_type surfArea = 0.0;
    for (index_type i = 0; i < n(); ++i) {
        if (_hasChildren[i]) {
            _area[i] = 0.0;
        }
        else {
            surfArea += computeArea(i);
        }
    }
    std::stringstream ss;
    ss << "surface area reset, surfArea = " << surfArea;
    OutputMessage statusMsg(ss.str(), OutputMessage::remarkPriority, "Faces::resetAreas");
    log->logMessage(statusMsg);
}

bool Faces::verifyConnectivity(const index_type i) const {
    std::vector<std::string> msgStrings;
//     std::stringstream ss;
    bool edgesConnect = true;
//     ss << "Face " << i << ":\n";
    const std::vector<index_type> edgeList = _edgeInds[i];
    for (index_type j = 0; j < edgeList.size(); ++j) {
        if (!(edges->rightFace(edgeList[j]) == i || edges->leftFace(edgeList[j]) == i) ) {
//             ss << "\tdoes not connect to edge " << edgeList[j] << ", face-relative index " << j << std::endl;
            edgesConnect = false;
        }
    }
    std::vector<index_type> verts;
    index_type vert2;
    for (index_type j = 0; j < _edgeInds[i].size(); ++j) {
        if (edgeIsPositive(i, edgeList[j]) ) {
            verts.push_back(edges->orig(edgeList[j]));
            vert2 = edges->dest(edgeList[j]);
        }
        else {
            verts.push_back(edges->dest(edgeList[j]));
            vert2 = edges->orig(edgeList[j]);
        }
    }
    bool isClosed = (verts[0] == vert2);
//     if (!isClosed) {
//         ss << "\tvertices are not closed\n";
//     }
    bool result = (edgesConnect && isClosed);
//     if (!result) {
//         OutputMessage errMsg(ss.str(), OutputMessage::errorPriority, "Faces::verifyConnectivity");
//         log->logMessage(errMsg);
//         throw std::runtime_error("Faces::verifyConnectivity");
//     }
    return result;
}

std::string Faces::faceRecord(const index_type ind) const {
    std::stringstream ss;
    const std::vector<index_type> edgeList(_edgeInds[ind]);
    ss << "face record " << ind << ":" << std::endl;
    log->startSection();
    for (index_type j = 0; j < edgeList.size(); ++j)
        ss << edges->edgeRecord(_edgeInds[ind][j]) << std::endl;
    ss << std::endl;
    log->endSection();
    return ss.str();
}

std::vector<index_type> Faces::vertexIndices(const index_type i) const {
    const bool correctConnectivity = verifyConnectivity(i);
    if (!correctConnectivity) {
        std::stringstream ss;
        ss << faceRecord(i);
        OutputMessage errMsg(ss.str(), OutputMessage::errorPriority, "Faces::vertexIndices");
        log->logMessage(errMsg);
        throw std::runtime_error("Faces::verifyConnectivity failed");
    }
    std::vector<index_type> verts;
    const std::vector<index_type> edgeList(_edgeInds[i]);
    for (index_type j = 0; j < edgeList.size(); ++j) {
        if (edgeIsPositive(i, edgeList[j]) ) {
            verts.push_back(edges->orig(edgeList[j]));
        }
        else {
            verts.push_back(edges->dest(edgeList[j]));
        }
    }
    return verts;
}


}