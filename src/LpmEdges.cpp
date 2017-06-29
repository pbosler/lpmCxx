#include "LpmEdges.h"
#include <algorithm>
#include <limits>

namespace Lpm {

std::unique_ptr<Logger> Edges::log(new Logger(OutputMessage::debugPriority));

Edges::Edges(const index_type nMax, const std::shared_ptr<Coords> crd_ptr, const std::shared_ptr<Coords> lag_crd_ptr) : 
    _nMax(nMax), crds(crd_ptr), lagCrds(lag_crd_ptr) {
    _orig.reserve(nMax);
    _dest.reserve(nMax);
    _rightFace.reserve(nMax);
    _leftFace.reserve(nMax);
    _child0.reserve(nMax);
    _child1.reserve(nMax);
    _hasChildren.reserve(nMax);
    _parent.reserve(nMax);
}

index_type Edges::nDivided() const {
    return std::count(_hasChildren.begin(), _hasChildren.end(), true);
} 

index_type Edges::nLeaves() const {
    return n() - nDivided();
}

scalar_type Edges::length(const index_type i, const bool lagrangian) const {
    return (lagrangian ? lagCrds->distance(_orig[i], _dest[i]) : crds->distance(_orig[i], _dest[i]));
}

scalar_type Edges::minLength(const bool lagrangian) const {
    scalar_type result = std::numeric_limits<scalar_type>::max();
    for (index_type i = 0; i < n(); ++i) {
        if (!_hasChildren[i]) {
            const scalar_type ll = length(i, lagrangian);
            if (ll < result) 
                result = ll;
        }
    }
    return result;
}

scalar_type Edges::maxLength(const bool lagrangian) const {
    scalar_type result = 0.0;
    for (index_type i = 0; i < n(); ++i) {
        if (!_hasChildren[i]) {
            const scalar_type ll = length(i, lagrangian);
            if (ll > result)
                result = ll;
        }
    }
    return result;
}

XyzVector Edges::midpoint(const index_type i, const bool lagrangian) const {
    return (lagrangian ? lagCrds->midpoint(_orig[i], _dest[i]) : crds->midpoint(_orig[i], _dest[i]));
}

XyzVector Edges::origCoord(const index_type i, const bool lagrangian) const {
    return (lagrangian ? lagCrds->getVec(_orig[i]) : crds->getVec(_orig[i]));
}

XyzVector Edges::destCoord(const index_type i, const bool lagrangian) const {
    return (lagrangian ? lagCrds->getVec(_dest[i]) : crds->getVec(_dest[i]));
}

XyzVector Edges::edgeVector(const index_type i, const bool lagrangian) const {
    XyzVector origVec;
    XyzVector destVec;
    if (lagrangian) {
        origVec = lagCrds->getVec(_orig[i]);
        destVec = lagCrds->getVec(_dest[i]);
    }
    else {
        origVec = crds->getVec(_orig[i]);
        destVec = crds->getVec(_dest[i]);
    }
    return destVec - origVec;
}

void Edges::divide(const index_type i) {
    const XyzVector midpt = midpoint(i);
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
        const XyzVector lagMidpt = this->midpoint(i, true);
        lagCrds->insert(lagMidpt);
    }
}

void Edges::insert(const index_type origInd, const index_type destInd, const index_type leftInd, const index_type rightInd) {
    if (n() + 1 > _nMax) {
        OutputMessage errMessage("not enough memory", OutputMessage::errorPriority, "Edges::insert");
        log->logMessage(errMessage);
    }
    _orig.push_back(origInd);
    _dest.push_back(destInd);
    _leftFace.push_back(leftInd);
    _rightFace.push_back(rightInd);
    _child0.push_back(-1);
    _child1.push_back(-1);
    _hasChildren.push_back(false);
    _parent.push_back(-1);
}

}
