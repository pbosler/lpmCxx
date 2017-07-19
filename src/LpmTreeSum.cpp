#include "LpmTreeSum.h"
#include <cmath>

namespace Lpm {

TaylorCoeffs::TaylorCoeffs(const int maxSeriesOrder, const XyzVector& tgtVec, const XyzVector& centroid) {
    for (int k1 = 0; k1 <= maxSeriesOrder; ++k1) {
        for (int k2 = 0; k2 <= maxSeriesOrder; ++k2) {
            for (int k3 = 0; k3 <= maxSeriesOrder; ++k3) {
                const MultiIndex kInd(k1, k2, k3);
                if (kInd.magnitude() <= maxSeriesOrder) {
                    vals.emplace(kInd, 0.0);
                }
            }
        }
    }
};

MultiIndex::vectorPower(const XyzVector& vec) const {
    return XyzVector(std::pow(vec.x, k1), std::pow(vec.y, k2), std::pow(vec.z, k3));
}

TreeSumNode::TreeSumNode(const std::shared_ptr<Coords> crds, const scalar_type maxAspectRatio, const int maxSeriesOrder):
    Treenode(crds, maxAspectRatio), coeffs(maxSeriesOrder) {
    for (auto& elem : coeffs.vals) {
        momentsA.emplace(elem.first, 0.0);
        momentsB.emplace(elem.first, 0.0);
        momentsC.emplace(elem.first, 0.0);
    }
}

TreeSumNode::computeMoments(const std::shared_ptr<Coords> crds, const std::shared_ptr<Field> srcStrength) {
    const XyzVector cntd = centroid();
    for (auto& elem : coeffs.vals) {
        for (index_type j = 0; j < coordsContained.size(); ++j) {
            const XyzVector srcVec = crds->getVec(coordsContained[j]);
            const XyzVector momVec = elem.first.vectorPower(srcVec - cntd);
            
            momentsA.at(elem.first) += srcStrength->getScalar(coordsContained[j]) * srcVec.x * 
        }
    }
}


bool operator < (const MultiIndex& left, const MultiIndex& right) {
    if (left.k1 != right.k1) 
        return left.k1 < right.k1;
    else {
        if (left.k2 != right.k2)
            return left.k2 < right.k2;
        else
            return left.k3 < right.k3;
    }
}

bool operator > (const MultiIndex& left, const MultiIndex& right) {
    return right < left;
}

bool operator >= (const MultiIndex& left, const MultiIndex& right) {
    return !(left < right);
}

bool operator <= (const MultiIndex& left, const MultiIndex& right) {
    return !(right < left);
}

bool operator == (const MultiIndex& left, const MultiIndex& right) {
    return (left.k1 == right.k1 && left.k2 == right.k2) && left.k3 == right.k3;
}

bool operator != (const MultiIndex& left, const MultiIndex& right) {
    return !(left == right);
}

}