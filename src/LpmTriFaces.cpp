#include "LpmTriFaces.h"

namespace Lpm {

std::unique_ptr<Logger> Faces::log(new Logger(OutputMessage::debugPriority));

TriFaces::TriFaces(const index_type nMax, const std::shared_ptr<Edges> edge_ptr, const std::shared_ptr<Coords> crd_ptr, 
    const bool sim3d) : Faces(nMax, 3, edge_ptr, crd_ptr, sim3d) {};

void TriFaces::divide(const index_type i) {
    if ( _nMax < n() + 4 ) {
        OutputMessage errMsg("not enough memory", OutputMessage::errorPriority, "TriFaces::divide");
        log->logMessage(errMsg);
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
        
        if (j == 0) {
            newFaceVertInds[0][1] = edges->dest(edgeChildren.first);
            newFaceVertInds[1][0] = edges->dest(edgeChildren.first);
        }
        else if (j == 1) {
            newFaceVertInds[1][2] = edges->dest(edgeChildren.first);
            newFaceVertInds[2][1] = edges->dest(edgeChildren.first);
        }
        else if ( j == 2) {
            newFaceVertInds[2][0] = edges->dest(edgeChildren.first);
            newFaceVertInds[0][2] = edges->dest(edgeChildren.first);
        }
    }
 
    //
    //  new interior edges
    //
    const index_type edgeInsertPoint = edges->n();
    
    
}



}
