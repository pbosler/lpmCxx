#include "LpmEdges.h"
#include <algorithm>
#include <limits>
#include <exception>
#include <sstream>

namespace Lpm {

std::unique_ptr<Logger> Edges::log(new Logger(OutputMessage::debugPriority, "Edges_log"));

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

std::string Edges::edgeRecord(const index_type i) const {
    std::stringstream ss;
    ss << "edge record " << i << ": orig = " << _orig[i] << ", dest = " << _dest[i] << ", leftFace = " << _leftFace[i]
        << ", rightFace = " << _rightFace[i] << ", parent = " << _parent[i];
    if (_hasChildren[i]) {
        ss << ", children = (" << _child0[i] << ", " << _child1[i] << ")";
    }
    else {
        ss << ", undivided." << std::endl;
    }
    ss << "\t\torigCoord = " << crds.lock()->getVec(_orig[i]) << ", destCoord = " << crds.lock()->getVec(_dest[i]) << std::endl;
    return ss.str();
}

index_type Edges::nDivided() const {
    return std::count(_hasChildren.begin(), _hasChildren.end(), true);
} 

index_type Edges::nLeaves() const {
    return n() - nDivided();
}

scalar_type Edges::length(const index_type i, const bool lagrangian) const {
    return (lagrangian ? lagCrds.lock()->distance(_orig[i], _dest[i]) : crds.lock()->distance(_orig[i], _dest[i]));
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
    return (lagrangian ? lagCrds.lock()->midpoint(_orig[i], _dest[i]) : crds.lock()->midpoint(_orig[i], _dest[i]));
}

XyzVector Edges::origCoord(const index_type i, const bool lagrangian) const {
    return (lagrangian ? lagCrds.lock()->getVec(_orig[i]) : crds.lock()->getVec(_orig[i]));
}

XyzVector Edges::destCoord(const index_type i, const bool lagrangian) const {
    return (lagrangian ? lagCrds.lock()->getVec(_dest[i]) : crds.lock()->getVec(_dest[i]));
}

XyzVector Edges::edgeVector(const index_type i, const bool lagrangian) const {
    XyzVector origVec;
    XyzVector destVec;
    if (lagrangian) {
        origVec = lagCrds.lock()->getVec(_orig[i]);
        destVec = lagCrds.lock()->getVec(_dest[i]);
    }
    else {
        origVec = crds.lock()->getVec(_orig[i]);
        destVec = crds.lock()->getVec(_dest[i]);
    }
    return destVec - origVec;
}

void Edges::divide(const index_type i) {
    const XyzVector midpt = midpoint(i);
    const index_type crdInsertPoint = crds.lock()->n();
    const index_type edgeInsertPoint = _orig.size();
    
    crds.lock()->insert(midpt);
    
    _child0[i] = edgeInsertPoint;
    _child1[i] = edgeInsertPoint + 1;
    _hasChildren[i] = true;
    
    insert(_orig[i], crdInsertPoint, _leftFace[i], _rightFace[i]);
    _parent[edgeInsertPoint] = i;
    
    insert(crdInsertPoint, _dest[i], _leftFace[i], _rightFace[i]);
    _parent[edgeInsertPoint+1] = i;
    
    if (!lagCrds.expired()) {
        const XyzVector lagMidpt = this->midpoint(i, true);
        lagCrds.lock()->insert(lagMidpt);
    }
}

void Edges::insert(const index_type origInd, const index_type destInd, const index_type leftInd, const index_type rightInd) {
    if (n() + 1 > _nMax) {
        std::stringstream ss;
        ss << "not enough memory to insert edge: nMax = " << _nMax << ", n = " << n();
        OutputMessage errMessage(ss.str(), OutputMessage::errorPriority, "Edges::insert");
        log->logMessage(errMessage);
        throw std::bad_alloc();
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
