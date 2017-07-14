#include "LpmCoords.h"
#include "LpmXyzVector.h"
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <exception>

namespace Lpm {

std::unique_ptr<Logger> Coords::log(new Logger(OutputMessage::debugPriority, "Coords_base_log"));

Coords::Coords(const index_type nMax, const bool sim3d) : _nMax(nMax) {
    x.reserve(nMax);
    y.reserve(nMax);
    if (sim3d)
        z.reserve(nMax);
}

XyzVector Coords::crossProduct(const index_type indA, const index_type indB) const {
    XyzVector vecA;
    XyzVector vecB;
    if (_geometry == PLANAR_GEOMETRY) {
        vecA = XyzVector(x[indA], y[indA]);
        vecB = XyzVector(x[indB], y[indB]);
    }
    else {
        vecA = XyzVector(x[indA], y[indA], z[indA]);
        vecB = XyzVector(x[indB], y[indB], z[indB]);
    }
    return vecA.crossProduct(vecB);
}

scalar_type Coords::dotProduct(const index_type indA, const index_type indB) const {
    XyzVector vecA;
    XyzVector vecB;
    if (_geometry == PLANAR_GEOMETRY) {
        vecA = XyzVector(x[indA], y[indA]);
        vecB = XyzVector(x[indB], y[indB]);
    }
    else {
        vecA = XyzVector(x[indA], y[indA], z[indA]);
        vecB = XyzVector(x[indB], y[indB], z[indB]);
    }    return vecA.dotProduct(vecB);
}

void Coords::scaleAll(const scalar_type multiplier) {
    if (_geometry == PLANAR_GEOMETRY) {
        for (index_type i = 0; i < x.size(); ++i) {
            x[i] *= multiplier;
            y[i] *= multiplier;
        }
    }
    else{
        for (index_type i = 0; i < x.size(); ++i) {
            x[i] *= multiplier;
            y[i] *= multiplier;
            z[i] *= multiplier;
        }
    }
}

std::vector<XyzVector> Coords::getVectors(const std::vector<index_type> inds) const {
    std::vector<XyzVector> result;
    for (index_type i = 0; i < inds.size(); ++i) {
        result.push_back(XyzVector(x[inds[i]], y[inds[i]], (_geometry == PLANAR_GEOMETRY ? 0.0 : z[inds[i]])));
    }
    return result;
}


void Coords::normalizeAll() {
    if (_geometry == PLANAR_GEOMETRY) {
            for (index_type i = 0; i < x.size(); ++i) {
                const scalar_type norm = std::sqrt(x[i] * x[i] + y[i] * y[i]);
                x[i] /= norm;
                y[i] /= norm;
            }
    }
    else {
        for (index_type i = 0; i < x.size(); ++i) {
            const scalar_type norm = std::sqrt(x[i] * x[i] + y[i] * y[i] + z[i] * z[i]);
            x[i] /= norm;
            y[i] /= norm;
            z[i] /= norm;
        }
    }
}

scalar_type Coords::magnitude(const index_type ind) const {
    const XyzVector vec(x[ind], y[ind], (_geometry == PLANAR_GEOMETRY ? 0.0 : z[ind]));
    return vec.magnitude();
}

void Coords::replace(const index_type ind, const scalar_type nx, const scalar_type ny, const scalar_type nz) {
    x[ind] = nx;
    y[ind] = ny;
    if (_geometry != PLANAR_GEOMETRY)
        z[ind] = nz;
}

void Coords::replace(const index_type ind, const XyzVector& vec) {
    x[ind] = vec.x;
    y[ind] = vec.y;
    if (_geometry != PLANAR_GEOMETRY)
        z[ind] = vec.z;
}

void Coords::swap(const index_type i, const index_type j) {
    const scalar_type tmpx = x[i];
    const scalar_type tmpy = y[i];
    x[i] = x[j];
    y[i] = y[j];
    x[j] = tmpx;
    y[j] = tmpy;
    if (_geometry != PLANAR_GEOMETRY) {
        const scalar_type tmpz = z[i];
        z[i] = z[j];
        z[j] = tmpz;
    }
}

void Coords::writeCoords(std::ostream& os) const {
    if ( _geometry == PLANAR_GEOMETRY) {
        for (index_type i = 0; i < n(); ++i) {
            os << x[i] << " " << y[i] << std::endl;
        }
    }
    else {
        for (index_type i = 0; i < n(); ++i) {
            os << x[i] << " " << y[i] << " " << z[i] << std::endl;
        }
    }
}

void Coords::writeCoordsCSV(std::ostream& os) const {
    if ( _geometry == PLANAR_GEOMETRY) {
        os << "x,y" << std::endl;
        for (index_type i = 0; i < n(); ++i) {
            os << x[i] << "," << y[i] << std::endl;
        }
    }
    else {
        os << "x,y,z" << std::endl;
        for (index_type i = 0; i < n(); ++i) {
            os << x[i] << "," << y[i] << "," << z[i] << std::endl;
        }
    }
}

void Coords::insert(const scalar_type nx, const scalar_type ny, const scalar_type nz) {
    if (n() + 1 > _nMax) {
        std::stringstream ss;
        ss << "not enough memory to insert coordinate " << XyzVector(nx, ny, nz);
        OutputMessage errMsg(ss.str(), OutputMessage::errorPriority, "Coords::insert");
        log->logMessage(errMsg);
        throw std::bad_alloc();
    }  
    x.push_back(nx);
    y.push_back(ny);
    if (_geometry != PLANAR_GEOMETRY) 
        z.push_back(nz);    
}

void Coords::insert(const XyzVector& vec) {
    if (n() + 1 > _nMax) {
        std::stringstream ss;
        ss << "not enough memory to insert coordinate " << vec;
        OutputMessage errMsg(ss.str(), OutputMessage::errorPriority, "Coords::insert");
        log->logMessage(errMsg);
        throw std::bad_alloc();
    }
    x.push_back(vec.x);
    y.push_back(vec.y);
    if (_geometry != PLANAR_GEOMETRY)
        z.push_back(vec.z);
}

std::string Coords::listAllCoords() const {
    std::ostringstream ss;
    for (index_type i = 0; i < x.size(); ++i)
        ss << i << ": " << XyzVector(x[i], y[i], (_geometry == PLANAR_GEOMETRY ? 0.0 : z[i]) ) << std::endl;
    return ss.str();
}

}
