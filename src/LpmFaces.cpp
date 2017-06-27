#include "LpmFaces.h"
#include "LpmOutputMessage.h"
#include <iostream>
#include <sstream>

namespace Lpm {

std::unique_ptr<Logger> Faces::log(new Logger(OutputMessage::debugPriority));

Faces::Faces(const index_type nMax, const index_type nMaxEdgesPerFace, const std::shared_ptr<Edges> edge_ptr, 
    const std::shared_ptr<Coords> crd_ptr,  const bool sim3d) : 
    _nMax(nMax), _nLeaves(0), _nMaxEdges(nMaxEdgesPerFace), edges(edge_ptr), crds(crd_ptr) {
    _edgeInds.reserve(nMax);
    for (index_type i = 0; i < nMax; ++i) {
        _edgeInds[i].reserve(nMaxEdgesPerFace);
    }        
    _area.reserve(nMax);
    _hasChildren.reserve(nMax);
    _parent.reserve(nMax);
    _children.reserve(nMax);
    if (sim3d) {
        _positiveCell.reserve(nMax);
        _negativeCell.reserve(nMax);
    }
}

XyzVector Faces::centroid(const index_type i) const {
    const std::vector<index_type> verts = vertexIndices(i);
    return crds->centroid(verts);    
}

bool Faces::verifyConnectivity(const index_type i) const {
    std::vector<std::string> msgStrings;
    std::stringstream ss;
    bool edgesConnect = true;
    ss << "Face " << i << ":\n";
    const std::vector<index_type> edgeList = _edgeInds[i];
    for (index_type j = 0; j < edgeList.size(); ++j) {
        if (!(edges->rightFace(edgeList[j]) == i || edges->leftFace(edgeList[j]) == i) ) {
            ss << "\tdoes not connect to edge " << edgeList[j] << ", face-relative index " << j << std::endl;
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
    if (!isClosed) {
        ss << "\tvertices are not closed\n";
    }
    bool result = (edgesConnect && isClosed);
    if (!result) {
        OutputMessage errMsg(ss.str(), OutputMessage::errorPriority, "Faces::verifyConnectivity");
        log->logMessage(errMsg);
    }
    return result;
}

std::vector<index_type> Faces::vertexIndices(const index_type i) const {
#ifdef DEBUG_ALL
    const bool correctConnectivity = verifyConnectivity(i);
    if (!correctConnectivity) {
        OutputMessage errMsg("connectivity error", OutputMessage::errorPriority, "Faces::vertices");
        log->logMessage(errMsg);
    }
#endif
    std::vector<index_type> verts;
    for (index_type j = 0; j < _edgeInds[i].size(); ++j) {
        if (edgeIsPositive(i, j) ) {
            verts.push_back(edges->orig(_edgeInds[i][j]));
        }
        else {
            verts.push_back(edges->dest(_edgeInds[i][j]));
        }
    }
    return verts;
}


}
