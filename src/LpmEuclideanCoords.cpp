#include "LpmEuclideanCoords.h"
#include "LpmXyzVector.h"
#include <vector>

namespace Lpm {

scalar_type EuclideanCoords::distance(const index_type indA, const index_type indB) const {
    const XyzVector vecA(x[indA], y[indA], (geometry != PLANAR_GEOMETRY ? z[indA] : 0.0));
    const XyzVector vecB(x[indB], y[indB], (geometry != PLANAR_GEOMETRY ? z[indB] : 0.0));
    return Lpm::distance(vecA, vecB);
}

scalar_type EuclideanCoords::distance(const XyzVector& vec, const index_type ind) const {
    return Lpm::distance(vec, XyzVector(x[ind], y[ind], (geometry != PLANAR_GEOMETRY ? z[ind] : 0.0)));
}

scalar_type EuclideanCoords::distance(const XyzVector& v0, const XyzVector& v1) const {
    return Lpm::distance(v0, v1);
}

XyzVector EuclideanCoords::midpoint(const index_type indA, const index_type indB) const {
    const XyzVector vecA(x[indA], y[indA], (geometry != PLANAR_GEOMETRY ? z[indA] : 0.0));
    const XyzVector vecB(x[indB], y[indB], (geometry != PLANAR_GEOMETRY ? z[indB] : 0.0));
    return Lpm::midpoint(vecA, vecB);
}

XyzVector EuclideanCoords::centroid(const std::vector<index_type>& inds) const {
    std::vector<XyzVector> vecs;
    for (index_type i = 0; i < inds.size(); ++i)
        vecs.push_back(XyzVector(x[inds[i]], y[inds[i]], (geometry != PLANAR_GEOMETRY ? z[inds[i]] : 0.0)));
    return Lpm::centroid(vecs);
}

scalar_type EuclideanCoords::triArea(const index_type indA, const index_type indB, const index_type indC) const {
    const XyzVector vecA(x[indA], y[indA], (geometry != PLANAR_GEOMETRY ? z[indA] : 0.0));
    const XyzVector vecB(x[indB], y[indB], (geometry != PLANAR_GEOMETRY ? z[indB] : 0.0));
    const XyzVector vecC(x[indC], y[indC], (geometry != PLANAR_GEOMETRY ? z[indC] : 0.0));
    return Lpm::triArea(vecA, vecB, vecC);
}

scalar_type EuclideanCoords::triArea(const XyzVector& v0, const index_type indA, const index_type indB) const {
    const XyzVector v1(x[indA], y[indA], (geometry != PLANAR_GEOMETRY ? z[indA] : 0.0));
    const XyzVector v2(x[indB], y[indB], (geometry != PLANAR_GEOMETRY ? z[indB] : 0.0));
    return Lpm::triArea(v0, v1, v2);
}

}
