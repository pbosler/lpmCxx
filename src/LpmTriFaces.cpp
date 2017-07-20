#include "LpmTriFaces.h"
#include <exception>
#include <sstream>

namespace Lpm {

// std::unique_ptr<Logger> Faces::log(new Logger(OutputMessage::debugPriority));

TriFaces::TriFaces(const index_type nMax, const std::shared_ptr<Edges> edge_ptr, 
    const std::shared_ptr<Coords> crd_ptr, const std::shared_ptr<Coords> lag_crd_ptr,
    const bool sim3d) : Faces(nMax, 3, edge_ptr, crd_ptr, lag_crd_ptr, sim3d) {};

void TriFaces::divide(const index_type i) {
    if ( _nMax < n() + 4 ) {
        std::stringstream ss;
        ss << "not enough memory to divide face " << i;
        OutputMessage errMsg(ss.str(), OutputMessage::errorPriority, "TriFaces::divide");
        log->logMessage(errMsg);
        throw std::bad_alloc();
    }
    if (_hasChildren[i]) {
        std::stringstream ss;
        ss << "Faces::divide() called on face " << i << ", which already has children";
        OutputMessage errMsg(ss.str(), OutputMessage::errorPriority,
            "TriFaces::divide");
        log->logMessage(errMsg);
        throw std::runtime_error("mesh logic error");
    }
    
    std::shared_ptr<Edges> edge_ptr = edges.lock();

    std::vector<std::vector<index_type>> newFaceEdgeInds(4, std::vector<index_type>(3, -1));
    std::vector<std::vector<index_type>> newFaceVertInds(4, std::vector<index_type>(3, -1));
    const std::vector<index_type> newFaceInds = {n(), n() + 1, n() + 2, n() + 3};
    
    const std::vector<index_type> parentVerts = vertexIndices(i);
    for (int j = 0; j < 3; ++j)
        newFaceVertInds[j][j] = parentVerts[j];
    
    //
    //  loop over parent edges
    //
    for (int j = 0; j < 3; ++j) {
        const index_type parentEdge = _edgeInds[i][j];
        std::pair<index_type, index_type> edgeChildren;
        if (edge_ptr->isDivided(parentEdge)) {
            edgeChildren = edge_ptr->children(parentEdge);
        }
        else {
            edgeChildren.first = edge_ptr->n();
            edgeChildren.second = edge_ptr->n() + 1;
            
            edge_ptr->divide(parentEdge);
        }
        
        if (edgeIsPositive(i, parentEdge) ) {
            newFaceEdgeInds[j][j] = edgeChildren.first;
            edge_ptr->setLeftFace(edgeChildren.first, newFaceInds[j]);
            
            newFaceEdgeInds[(j+1)%3][j] = edgeChildren.second;
            edge_ptr->setLeftFace(edgeChildren.second, newFaceInds[(j+1)%3]);
        }
        else {
            newFaceEdgeInds[j][j] = edgeChildren.second;
            edge_ptr->setRightFace(edgeChildren.second, newFaceInds[j]);
            
            newFaceEdgeInds[(j+1)%3][j] = edgeChildren.first;
            edge_ptr->setRightFace(edgeChildren.first, newFaceInds[(j+1)%3]);
        }
        
        const index_type vertInd = edge_ptr->dest(edgeChildren.first);
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
    const index_type edgeInsertPoint = edge_ptr->n();
    newFaceEdgeInds[3] = {edgeInsertPoint, edgeInsertPoint + 1, edgeInsertPoint + 2};
    newFaceEdgeInds[0][1] = edgeInsertPoint + 1;
    newFaceEdgeInds[1][2] = edgeInsertPoint + 2;
    newFaceEdgeInds[2][0] = edgeInsertPoint;
    edge_ptr->insert(newFaceVertInds[2][1], newFaceVertInds[2][0], newFaceInds[3], newFaceInds[2]);
    edge_ptr->insert(newFaceVertInds[0][2], newFaceVertInds[0][1], newFaceInds[3], newFaceInds[0]);
    edge_ptr->insert(newFaceVertInds[1][0], newFaceVertInds[1][2], newFaceInds[3], newFaceInds[1]);
    //
    //  new child faces
    //
    for (int j = 0; j < 4; ++j) {
        insert(newFaceEdgeInds[j]);
    }
    _hasChildren[i] = true;
    _children[i] = newFaceInds;
    _area[i] = 0.0;
    for (int j = 0; j < 4; ++j)
        computeArea(newFaceInds[j]);
    if (_is3d) {
        const index_type posCellInd = _positiveCell[i];
        const index_type negCellInd = _negativeCell[i];
        for (int j = 0; j < 4; ++j ) {
            setPositiveCell(newFaceInds[j], posCellInd);
            setNegativeCell(newFaceInds[j], negCellInd);
        }
    }
    {
        std::stringstream ss;
        const bool parentConn = verifyConnectivity(i);
        bool errState = false;
        if (!parentConn) {
            ss << "Connectivity error in parent, faceindex " << i << std::endl;
            errState = true;
        }
        std::vector<bool> childConn(4, false);
        for (int j = 0; j < 4; ++j) {
            childConn[j] = verifyConnectivity(newFaceInds[j]);
            if (!childConn[j]) {
                ss << "connectivity error in child " << j << " of faceindex " << i << std::endl;
                errState = true;
            }
        }  
        if (errState){
            OutputMessage errMsg(ss.str(), OutputMessage::errorPriority, "TriFaces::divide");
            log->logMessage(errMsg);
            throw std::runtime_error("mesh connectivity error");
        }
    }
    edge_ptr.reset();
}



}
