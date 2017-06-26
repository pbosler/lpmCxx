#include "LpmEuclideanCoords.h"
#include "LpmXyzVector.h"
#include <vector>

namespace Lpm {

scalar_type EuclideanCoords::distance(const index_type indA, const index_type indB) const {
    const XyzVector vecA(x[indA], y[indA], z[indA]);
    const XyzVector vecB(x[indB], y[indB], z[indB]);
    return Lpm::distance(vecA, vecB);
}

XyzVector EuclideanCoords::midpoint(const index_type indA, const index_type indB) const {
    const XyzVector vecA(x[indA], y[indA], z[indA]);
    const XyzVector vecB(x[indB], y[indB], z[indB]);
    return Lpm::midpoint(vecA, vecB);
}

XyzVector EuclideanCoords::centroid(const std::vector<index_type>& inds) const {
    std::vector<XyzVector> vecs;
    for (index_type i = 0; i < inds.size(); ++i)
        vecs.push_back(XyzVector(x[inds[i]], y[inds[i]], z[inds[i]]));
    return Lpm::centroid(vecs);
}

scalar_type EuclideanCoords::triArea(const index_type indA, const index_type indB, const index_type indC) const {
    const XyzVector vecA(x[indA], y[indA], z[indA]);
    const XyzVector vecB(x[indB], y[indB], z[indB]);
    const XyzVector vecC(x[indC], y[indC], z[indC]);
    return Lpm::triArea(vecA, vecB, vecC);
}

}
