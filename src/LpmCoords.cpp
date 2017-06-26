#include "LpmCoords.h"
#include "LpmXyzVector.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <sstream>
#include <string>

namespace Lpm {

Coords::Coords(const index_type nMax) : _nMax(nMax) {
    x.reserve(nMax);
    y.reserve(nMax);
    z.reserve(nMax);
}

XyzVector Coords::crossProduct(const index_type indA, const index_type indB) const {
    const XyzVector vecA(x[indA], y[indA], z[indA]);
    const XyzVector vecB(x[indB], y[indB], z[indB]);
    return vecA.crossProduct(vecB);
}

scalar_type Coords::dotProduct(const index_type indA, const index_type indB) const {
    const XyzVector vecA(x[indA], y[indA], z[indA]);
    const XyzVector vecB(x[indB], y[indB], z[indB]);
    return vecA.dotProduct(vecB);
}

void Coords::scaleAll(const scalar_type multiplier) {
    for (index_type i = 0; i < x.size(); ++i) {
        x[i] *= multiplier;
        y[i] *= multiplier;
        z[i] *= multiplier;
    }
}

void Coords::normalizeAll() {
    for (index_type i = 0; i < x.size(); ++i) {
        const scalar_type norm = std::sqrt(x[i] * x[i] + y[i] * y[i] + z[i] * z[i]);
        x[i] /= norm;
        y[i] /= norm;
        z[i] /= norm;
    }
}

scalar_type Coords::magnitude(const index_type ind) const {
    const XyzVector vec(x[ind], y[ind], z[ind]);
    return vec.magnitude();
}

void Coords::replace(const index_type ind, const scalar_type nx, const scalar_type ny, const scalar_type nz) {
    x[ind] = nx;
    y[ind] = ny;
    z[ind] = nz;
}

void Coords::replace(const index_type ind, const XyzVector& vec) {
    x[ind] = vec.x;
    y[ind] = vec.y;
    z[ind] = vec.z;
}

void Coords::insert(const scalar_type nx, const scalar_type ny, const scalar_type nz) {
    x.push_back(nx);
    y.push_back(ny);
    z.push_back(nz);
}

void Coords::insert(const XyzVector& vec) {
    x.push_back(vec.x);
    y.push_back(vec.y);
    z.push_back(vec.z);
}

std::string Coords::listAllCoords() const {
    std::ostringstream ss;
    for (index_type i = 0; i < x.size(); ++i)
        ss << i << ": " << XyzVector(x[i], y[i], z[i]) << std::endl;
    return ss.str();
}

}
