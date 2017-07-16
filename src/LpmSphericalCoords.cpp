#include "LpmSphericalCoords.h"
#include "LpmXyzVector.h"
#include <cmath>
#include <random>
#include <chrono>

namespace Lpm {

scalar_type SphericalCoords::distance(const index_type indA, const index_type indB) const {
    const XyzVector vecA(x[indA], y[indA], z[indA]);
    const XyzVector vecB(x[indB], y[indB], z[indB]);
    return Lpm::sphereDistance(vecA, vecB, _radius);
}

scalar_type SphericalCoords::distance(const XyzVector& vec, const index_type ind) const {
    return Lpm::sphereDistance(vec, XyzVector(x[ind], y[ind], z[ind]), _radius);
}

scalar_type SphericalCoords::distance(const XyzVector& v0, const XyzVector& v1) const {
    return Lpm::sphereDistance(v0, v1);
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

scalar_type SphericalCoords::triArea(const XyzVector& v0, const index_type indA, const index_type indB) const {
    const XyzVector v1(x[indA], y[indA], z[indA]);
    const XyzVector v2(x[indB], y[indB], z[indB]);
    return Lpm::sphereTriArea(v0, v1, v2, _radius);
}

scalar_type SphericalCoords::latitude(const index_type ind) const {
    return std::atan2(z[ind], std::sqrt(x[ind]*x[ind] + y[ind]*y[ind]));    
}

scalar_type SphericalCoords::longitude(const index_type ind) const {
    return atan4(y[ind], x[ind]);
}

void SphericalCoords::initRandom(const bool useTimeSeed, const scalar_type domainRadius) {
    unsigned seed(1234);
    if (useTimeSeed) {
        seed = std::chrono::system_clock::now().time_since_epoch().count();
    }
    std::default_random_engine generator(seed);
    
    std::uniform_real_distribution<scalar_type> randDist(-1.0, 1.0);
    for (index_type i = 0; i < _nMax; ++i) {
        scalar_type uu = randDist(generator);
        scalar_type vv = randDist(generator);
        while (uu * uu + vv * vv > 1.0) {
            uu = randDist(generator);
            vv = randDist(generator);
        }
        const scalar_type uvrad = std::sqrt(1.0 - uu * uu - vv * vv);
        const scalar_type xx = 2.0 * uu * uvrad * domainRadius;
        const scalar_type yy = 2.0 * vv * uvrad * domainRadius;
        const scalar_type zz = (1.0 - 2.0 * (uu * uu + vv * vv)) * domainRadius;
        insert(xx, yy, zz);
    }
}

}
