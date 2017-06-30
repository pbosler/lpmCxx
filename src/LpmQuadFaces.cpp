#include "LpmQuadFaces.h"
#include "LpmOutputMessage.h"
#include "LpmXyzVector.h"
#include <vector>
#include <sstream>
#include <exception>

namespace Lpm {

QuadFaces::QuadFaces(const index_type nMax, const std::shared_ptr<Edges> edge_ptr, 
    const std::shared_ptr<Coords> crd_ptr, const std::shared_ptr<Coords> lag_crd_ptr, 
        const bool sim3d) : Faces(nMax, 4, edge_ptr, crd_ptr, lag_crd_ptr, sim3d) {};

void QuadFaces::divide(const index_type i) {
    if (_nMax < n() + 4) {
        std::stringstream ss;
        ss << "not enough memory to divide face " << i;
        OutputMessage errMsg(ss.str(), OutputMessage::errorPriority, "QuadFaces::divide");
        log->logMessage(errMsg);
        throw std::bad_alloc();
    }
    if (_hasChildren[i]) {
        std::stringstream ss;
        ss << "Faces::divide() called on face " << i << ", which already has children";
        OutputMessage errMsg(ss.str(), OutputMessage::errorPriority, "QuadFaces::divide");
        log->logMessage(errMsg);
        throw std::runtime_error("mesh logic error");
    }
    std::vector<std::vector<index_type>> newFaceEdgeInds(4, std::vector<index_type>(4,-1));
    std::vector<std::vector<index_type>> newFaceVertInds(4, std::vector<index_type>(4,-1));
    const std::vector<index_type> newFaceInds = {n(), n() + 1, n() + 2, n() + 3};
    
    const std::vector<index_type> parentVerts = vertexIndices(i);
    for (int j = 0; j < 4; ++j) 
        newFaceVertInds[j][j] = parentVerts[j];
    
    //
    //  Loop over parent edges
    //
    for (int j = 0; j < 4; ++j) {
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
        
        if (edgeIsPositive(i, parentEdge)) {
            newFaceEdgeInds[j][j] = edgeChildren.first;
            edges->setLeftFace(edgeChildren.first, newFaceInds[j]);
            
            newFaceEdgeInds[(j+1)%4][j] = edgeChildren.second;
            edges->setLeftFace(edgeChildren.second, newFaceInds[(j+1)%4]);
        }
        else {
            newFaceEdgeInds[j][j] = edgeChildren.second;
            edges->setRightFace(edgeChildren.second, newFaceInds[j]);
            
            newFaceEdgeInds[(j+1)%4][j] = edgeChildren.first;
            edges->setRightFace(edgeChildren.first, newFaceInds[(j+1)%4]);
        }
        
        const index_type vertInd = edges->dest(edgeChildren.first);
        newFaceVertInds[j][(j+1)%4] = vertInd;
        newFaceVertInds[(j+1)%4][j] = vertInd; 
    }
    
    //
    //  new interior edges
    //
    const index_type edgeInsertPoint = edges->n();
    const index_type crdInsertPoint = crds->n();
    const XyzVector cntd = centroid(i);
    crds->insert(cntd);
    if (lagCrds){
        const XyzVector lagCntd = centroid(i, true);
        lagCrds->insert(lagCntd);
    }
    for (int j = 0; j < 4; ++j)
        newFaceVertInds[j][(j+2)%4] = crdInsertPoint;
    
    newFaceEdgeInds[0][1] = edgeInsertPoint;
    newFaceEdgeInds[0][2] = edgeInsertPoint + 3;
    newFaceEdgeInds[1][3] = edgeInsertPoint;
    newFaceEdgeInds[1][2] = edgeInsertPoint + 2;
    newFaceEdgeInds[2][3] = edgeInsertPoint + 1;
    newFaceEdgeInds[2][0] = edgeInsertPoint + 2;
    newFaceEdgeInds[3][1] = edgeInsertPoint + 1;
    newFaceEdgeInds[3][0] = edgeInsertPoint + 3;
    
    edges->insert(newFaceVertInds[0][1], newFaceVertInds[0][2], newFaceInds[0], newFaceInds[1]);
    edges->insert(newFaceVertInds[3][1], newFaceVertInds[3][2], newFaceInds[3], newFaceInds[2]);
    edges->insert(newFaceVertInds[1][2], newFaceVertInds[1][3], newFaceInds[1], newFaceInds[2]);
    edges->insert(newFaceVertInds[0][2], newFaceVertInds[0][3], newFaceInds[0], newFaceInds[3]);
    
    //
    //  new interior faces
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
    
#ifdef DEBUG_ALL
    std::cout << "dividing face " << i << ":\n";
    std::cout << "\tnewFaceEdgeInds = \n";
    for (int j = 0; j < 4; ++j) {
        std::cout << "\t" << j << ": ";
        for (int k = 0; k < 3; ++k)
            std::cout << newFaceEdgeInds[j][k] << ", ";
        std::cout << newFaceEdgeInds[j][3] << std::endl;
    }
    std::cout << "\tnewFaceVertInds = \n";
    for (int j = 0; j < 4; ++j) {
        std::cout << "\t" << j << ": ";
        for (int k = 0; k < 3; ++k)
            std::cout << newFaceVertInds[j][k] << ", ";
        std::cout << newFaceVertInds[j][3] << std::endl;
    }
    std::cout << "verifying connectivity: ";
    std::cout << (verifyConnectivity(i) ? " parent ok" : " connectivity error in parent");
    for (int j = 0; j < 4; ++j) {
        std::cout << (verifyConnectivity(newFaceInds[j]) ? " child ok " : " connectivity error in child");
    }
    std::cout << std::endl;
#endif
}

}
