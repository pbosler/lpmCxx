#include "LpmEdges.h"

namespace Lpm {

std::unique_ptr<Logger> Edges::log(new Logger(OutputMessage::debugPriority));

Edges::Edges(const index_type nMax) : _nMax(nMax), _nLeaves(0){
    _orig.reserve(nMax);
    _dest.reserve(nMax);
    _rightFace.reserve(nMax);
    _leftFace.reserve(nMax);
    _child0.reserve(nMax);
    _child1.reserve(nMax);
    _hasChildren.reserve(nMax);
    _parent.reserve(nMax);
}

void Edges::insert(const index_type origInd, const index_type destInd, const index_type leftInd, const index_type rightInd) {
    _orig.push_back(origInd);
    _dest.push_back(destInd);
    _leftFace.push_back(leftInd);
    _rightFace.push_back(rightInd);
    _child0.push_back(-1);
    _child1.push_back(-1);
    _hasChildren.push_back(false);
    _parent.push_back(-1);
    _nLeaves += 1;
}

XyzVector Edges::midpoint(const index_type i, const Coords* crds) const {
    return crds->midpoint(_orig[i], _dest[i]);
}

scalar_type Edges::length(const index_type i, const Coords* crds) const {
    return crds->distance(_orig[i], _dest[i]);
}

XyzVector Edges::origCoord(const index_type i, const Coords* crds) const {
    return crds->getVec(_orig[i]);
}

XyzVector Edges::destCoord(const index_type i, const Coords* crds) const {
    return crds->getVec(_dest[i]);
}

XyzVector Edges::edgeVector(const index_type i, const Coords* crds) const {
    const XyzVector origVec = crds->getVec(_orig[i]);
    const XyzVector destVec = crds->getVec(_dest[i]);
    return destVec - origVec;
}

void Edges::divide(const index_type i, Coords* crds, Coords* lagCrds) {
    const XyzVector midpt = this->midpoint(i, crds);
    const index_type crdInsertPoint = crds->n();
    const index_type edgeInsertPoint = _orig.size();
    
    crds->insert(midpt);
    
    _child0[i] = edgeInsertPoint;
    _child1[i] = edgeInsertPoint + 1;
    _hasChildren[i] = true;
    
    insert(_orig[i], crdInsertPoint, _leftFace[i], _rightFace[i]);
    _parent[edgeInsertPoint] = i;
    
    insert(crdInsertPoint, _dest[i], _leftFace[i], _rightFace[i]);
    _parent[edgeInsertPoint+1] = i;
    
    if (lagCrds) {
        const XyzVector lagMidpt = this->midpoint(i, lagCrds);
        lagCrds->insert(lagMidpt);
    }
    
    _nLeaves -= 1;
}



}
