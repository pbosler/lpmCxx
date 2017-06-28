#include "LpmTriFaces.h"

namespace Lpm {

// std::unique_ptr<Logger> Faces::log(new Logger(OutputMessage::debugPriority));

TriFaces::TriFaces(const index_type nMax, const std::shared_ptr<Edges> edge_ptr, 
    const std::shared_ptr<Coords> crd_ptr, const std::shared_ptr<Coords> lag_crd_ptr,
    const bool sim3d) : Faces(nMax, 3, edge_ptr, crd_ptr, lag_crd_ptr, sim3d) {};

// void TriFaces::insert(const index_type indA, const index_type indB, const index_type indC) {
//     if ( n() + 1 >= _nMax) {
//         OutputMessage errMsg("not enough memory", OutputMessage::errorPriority, "TriFaces::insert");
//         log->logMessage(errMsg);
//         return;
//     }
//     std::vector<index_type> edgeList = {indA, indB, indC};
//     _edgeInds.push_back(edgeList);
// }

void TriFaces::divide(const index_type i) {
    bool errStat = false;
    if ( _nMax < n() + 4 ) {
        OutputMessage errMsg("not enough memory", OutputMessage::errorPriority, "TriFaces::divide");
        log->logMessage(errMsg);
        errStat = true;
    }
    if (_hasChildren[i]) {
        OutputMessage errMsg("Faces::divide() called on face that already has children", OutputMessage::errorPriority,
            "TriFaces::divide");
        log->logMessage(errMsg);
        errStat = true;
    }
    if (errStat) {
        return;
    }
    std::vector<std::vector<index_type>> newFaceEdgeInds(4, std::vector<index_type>(3, -1));
    std::vector<std::vector<index_type>> newFaceVertInds(4, std::vector<index_type>(3, -1));
    std::vector<index_type> newFaceInds = {n(), n() + 1, n() + 2, n() + 3};
    
    const std::vector<index_type> parentVerts = vertexIndices(i);
    for (int j = 0; j < 3; ++j)
        newFaceVertInds[j][j] = parentVerts[j];
    
    //
    //  loop over parent edges
    //
    for (int j = 0; j < 3; ++j) {
        const index_type parentEdge = _edgeInds[i][j];
        std::pair<index_type, index_type> edgeChildren;
        if (edges->isDivided(parentEdge)) {
            edgeChildren = edges->children(parentEdge);
        }
        else {
            edgeChildren.first = edges->n();
            edgeChildren.second = edges->n() + 1;
            
            edges->divide(parentEdge);
        }
        
        if (edgeIsPositive(i, parentEdge) ) {
            newFaceEdgeInds[j][j] = edgeChildren.first;
            edges->setLeftFace(edgeChildren.first, newFaceInds[j]);
            
            newFaceEdgeInds[(j+1)%3][j] = edgeChildren.second;
            edges->setLeftFace(edgeChildren.second, newFaceInds[(j+1)%3]);
        }
        else {
            newFaceEdgeInds[j][j] = edgeChildren.second;
            edges->setRightFace(edgeChildren.second, newFaceInds[j]);
            
            newFaceEdgeInds[(j+1)%3][j] = edgeChildren.first;
            edges->setRightFace(edgeChildren.first, newFaceInds[(j+1)%3]);
        }
        
        const index_type vertInd = edges->dest(edgeChildren.first);
        if (j == 0) {
            newFaceVertInds[0][1] = vertInd;
            newFaceVertInds[1][0] = vertInd;
            newFaceVertInds[3][2] = vertInd;
        }
        else if (j == 1) {
            newFaceVertInds[1][2] = vertInd;
            newFaceVertInds[2][1] = vertInd;
            newFaceVertInds[3][0] = vertInd;
        }
        else if ( j == 2) {
            newFaceVertInds[2][0] = vertInd;
            newFaceVertInds[0][2] = vertInd;
            newFaceVertInds[3][1] = vertInd;
        }
    }
 
    //
    //  new interior edges
    //
    const index_type edgeInsertPoint = edges->n();
    newFaceEdgeInds[3] = {edgeInsertPoint, edgeInsertPoint + 1, edgeInsertPoint + 2};
    newFaceEdgeInds[0][1] = edgeInsertPoint + 1;
    newFaceEdgeInds[1][2] = edgeInsertPoint + 2;
    newFaceEdgeInds[2][0] = edgeInsertPoint;
    edges->insert(newFaceVertInds[2][1], newFaceVertInds[2][0], newFaceInds[3], newFaceInds[2]);
    edges->insert(newFaceVertInds[0][2], newFaceVertInds[0][1], newFaceInds[3], newFaceInds[0]);
    edges->insert(newFaceVertInds[1][0], newFaceVertInds[1][2], newFaceInds[3], newFaceInds[1]);
    //
    //  new child faces
    //
    for (int j = 0; j < 4; ++j) {
        insert(newFaceEdgeInds[j]);
    }
    _hasChildren[i] = true;
    _children[i] = quad_index_type(newFaceInds[0], newFaceInds[1], newFaceInds[2], newFaceInds[3]);
    _area[i] = 0.0;
    for (int j = 0; j < 4; ++j)
        computeArea(newFaceInds[j]);
    
    
#ifdef DEBUG_ALL
    std::cout << "dividing face " << i << ":\n";
    std::cout << "\tnewFaceEdges : \n";
    for (int j = 0; j < 4; ++j) {
        std::cout << "\t" << j << ": ";
        std::cout << newFaceEdgeInds[j][0] << ", " << newFaceEdgeInds[j][1] << ", " <<
            newFaceEdgeInds[j][2] << std::endl;
    }
    std::cout << "\tnewFaceVerts: \n";
    for (int j = 0; j < 4; ++j) {
        std::cout << "\t" << j << ": ";
        std::cout << newFaceVertInds[j][0] << ", " << newFaceVertInds[j][1] << ", " <<
            newFaceVertInds[j][2] << std::endl;
    }
    std::cout << "\tverifyingConnectivity: ";
    bool parentConn = verifyConnectivity(i);
    std::cout << (parentConn ? " parent ok " : " connectivity error in parent ");
    std::vector<bool> childConn(4, false);
    for (int j = 0; j < 4; ++j) {
        childConn[j] = verifyConnectivity(newFaceInds[j]);
        std::cout << (childConn[j] ? " child ok " : "connectivity error in child");
    }    
    std::cout << std::endl;
#endif
}



}
