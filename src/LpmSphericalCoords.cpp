#include "LpmSphericalCoords.h"
#include "LpmXyzVector.h"

namespace Lpm {

scalar_type SphericalCoords::distance(const index_type indA, const index_type indB) const {
    const XyzVector vecA(x[indA], y[indA], z[indA]);
    const XyzVector vecB(x[indB], y[indB], z[indB]);
    return Lpm::sphereDistance(vecA, vecB, _radius);
}

XyzVector SphericalCoords::midpoint(const index_type indA, const index_type indB) const {
    const XyzVector vecA(x[indA], y[indA], z[indA]);
    const XyzVector vecB(x[indB], y[indB], z[indB]);
    return Lpm::sphereMidpoint(vecA, vecB, _radius);
}

XyzVector SphericalCoords::centroid(const std::vector<index_type>& inds) const {
    std::vector<XyzVector> vecs;
    for (index_type i = 0; i < inds.size(); ++i)
        vecs.push_back(XyzVector(x[inds[i]], y[inds[i]], z[inds[i]]));
    return Lpm::sphereCentroid(vecs);
}

scalar_type SphericalCoords::triArea(const index_type indA, const index_type indB, const index_type indC) const {
    const XyzVector vecA(x[indA], y[indA], z[indA]);
    const XyzVector vecB(x[indB], y[indB], z[indB]);
    const XyzVector vecC(x[indC], y[indC], z[indC]);
    return Lpm::sphereTriArea(vecA, vecB, vecC, _radius);
}

}
