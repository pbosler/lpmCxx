#ifndef LPM_AOS_TYPES_HPP
#define LPM_AOS_TYPES_HPP

#include <cmath>
#include <iostream>
#include <vector>
#include <exception>
#include <array>
#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmUtilities.h"

namespace Lpm {

// scalar_type atan4(const scalar_type y, const scalar_type x);

typedef std::vector<index_type> ind_vec;

template<int ndim=3> struct Vec {
    scalar_type x[ndim];

    Vec() {
        for (int i=0; i<ndim; ++i)
            this->x[i] = 0.0;
    }
    
    template <int ndim2> Vec(const Vec<ndim2>& other){
        if (ndim < ndim2) {
            throw std::runtime_error("Vec<ndim>(Vec<ndim2>&) : cannot convert to smaller type.");
        }
        for (int i=0; i<std::min(ndim, ndim2); ++i) {
            this->x[i] = other.x[i];
        }
        for (int i=std::min(ndim, ndim2)-1; i<ndim; ++i) {
            this->x[i] = 0.0;
        }
    }

    Vec(const Vec<ndim>& other) {
        for (int i=0; i<ndim; ++i)
            this->x[i] = other.x[i];
    }

    Vec(const scalar_type xx, const scalar_type yy) {
        this->x[0] = xx;
        this->x[1] = yy;
    }

    Vec(const scalar_type xx, const scalar_type yy, const scalar_type zz) {
        if (ndim == 2) {
            this->x[0] = xx;
            this->x[1] = yy;
        }
        else if (ndim == 3) {
            this->x[0] = xx;
            this->x[1] = yy;
            this->x[2] = zz;
        }
    }

    Vec(const std::array<scalar_type, ndim>& arr) {
        for (int i=0; i<ndim; ++i)
            this->x[i] = arr[i];
    }
    
    Vec& operator= (const Vec<ndim>& other) {
        if (this != &other) {
            for (int i=0; i<ndim; ++i)
                this->x[i] = other.x[i];
        }
        return *this;
    }

    inline std::array<scalar_type, ndim> toArray() const {
        std::array<scalar_type, ndim> result;
        for (int i=0; i<ndim; ++i)
            result[i] = this->x[i];
        return result;
    }

    inline std::vector<scalar_type> toStdVec() const {
        std::vector<scalar_type> result(ndim);
        for (int i=0; i<ndim; ++i)
            result[i] = this->x[i];
    }

    inline Vec(const std::vector<scalar_type>& xx) {
        for (int i=0; i<ndim; ++i)
            this->x[i] = xx[i];
    }

    inline Vec(const scalar_type* xx) {
        for (int i=0; i<ndim; ++i)
            this->x[i] = xx[i];
    }

    inline Vec<ndim>& operator += (const Vec<ndim>& other) {
        for (int i=0; i<ndim; ++i)
            this->x[i] += other.x[i];
        return *this;
    }

    inline Vec<ndim>& operator -= (const Vec<ndim>& other) {
        for (int i=0; i<ndim; ++i)
            this->x[i] -= other.x[i];
        return *this;
    }

   inline  Vec<ndim>& operator *= (const Vec<ndim>& other) {
        for (int i=0; i<ndim; ++i)
            this->x[i] *= other.x[i];
        return *this;
    }

    inline Vec<ndim>& operator /= (const Vec<ndim>& other) {
        for (int i=0; i<ndim; ++i)
            this->x[i] /= other.x[i];
        return *this;
    }

    inline const Vec<ndim> operator + (const Vec<ndim>& other) const {
        return Vec<ndim>(*this) += other;
    }

    inline const Vec<ndim> operator - (const Vec<ndim>& other) const {
        return Vec<ndim>(*this) -= other;
    }

    inline const Vec<ndim> operator * (const Vec<ndim>& other) const {
        return Vec<ndim>(*this) *= other;
    }

    inline const Vec<ndim> operator / (const Vec<ndim>& other) const {
        return Vec<ndim>(*this) /= other;
    }

    inline scalar_type magSq() const {
        scalar_type sumsq(0.0);
        for (int i=0; i<ndim; ++i)
            sumsq += this->x[i]*this->x[i];
        return sumsq;
    }

    inline scalar_type mag() const {
        return std::sqrt(this->magSq());
    }

    inline scalar_type dotProd(const Vec<ndim>& other) const {
        scalar_type dp(0.0);
        for (int i=0; i<ndim; ++i)
            dp += this->x[i]*other.x[i];
        return dp;
    }

    inline scalar_type crossProdComp3(const Vec<2>& other) const {
        return this->x[0]*other.x[1] - this->x[1]*other.x[0];
    }

    inline const Vec<ndim> crossProd(const Vec<ndim>& other) const {
        const scalar_type cp[3] = {this->x[1]*other.x[2] - this->x[2]*other.x[1],
                                   this->x[2]*other.x[0] - this->x[0]*other.x[2],
                                   this->x[0]*other.x[1] - this->x[1]*other.x[0]};
        return Vec<ndim>(cp);
    }

    inline void scaleInPlace(const scalar_type mult) {
        for (int i=0; i<ndim; ++i)
            this->x[i] *= mult;
    }

    inline const Vec<ndim> scale(const scalar_type mult) const {
        scalar_type sm[ndim];
        for (int i=0; i<ndim; ++i)
            sm[i] = this->x[i]*mult;
        return Vec<ndim>(sm);
    }

    inline void normalizeInPlace() {
        const scalar_type len = this->mag();
        this->scaleInPlace(1.0/len);
    }

    inline const Vec<ndim> normalize() const {
        const scalar_type len = this->mag();
        return this->scale(1.0/len);
    }

    inline scalar_type longitude() const {return atan4(this->x[1], this->x[0]);}

    inline scalar_type latitude() const {
        const scalar_type xy2 = this->x[0]*this->x[0] + this->x[1]*this->x[1];
        return std::atan2(this->x[2] , std::sqrt(xy2));
    }

    inline const Vec<ndim> midpoint(const Vec<ndim>& other) const {
        Vec<ndim> result = *this + other;
        result.scaleInPlace(0.5);
        return result;
    }

    inline const Vec<ndim> sphereMidpoint(const Vec<ndim>& other, const scalar_type radius = 1.0) const {
        Vec<ndim> result = this->midpoint(other);
        result.normalizeInPlace();
        result.scaleInPlace(radius);
        return result;
    }

    inline scalar_type dist(const Vec<ndim>& other) const {
        return (*this - other).mag();
    }

    inline scalar_type sphereDist(const Vec<ndim>& other, const scalar_type radius=1.0) const {
        const Vec<ndim> cp = this->crossProd(other);
        const scalar_type dp = this->dotProd(other);
        return std::atan2(cp.mag(), dp) * radius;
    }

    inline const bool operator == (const Vec<ndim>& other) const {
        return this->dist(other) < ZERO_TOL;
    }
};

template <int ndim> Vec<ndim> pointAlongChord(Vec<ndim> a, Vec<ndim> b, const scalar_type s) {
    a.scaleInPlace(1.0-s);
    b.scaleInPlace(1.0+s);
    Vec<ndim> result = a + b;
    result.scaleInPlace(0.5);
    return result;
}

template <int ndim> Vec<ndim> pointAlongCircle(const Vec<ndim>& a, const Vec<ndim>& b, const scalar_type s, 
    const scalar_type radius=1.0) {
    Vec<ndim> result = pointAlongChord(a,b,s);
    result.normalizeInPlace();
    result.scaleInPlace(radius);
    return result;
}

template <int ndim> const Vec<ndim> baryCenter(const std::vector<Vec<ndim>>& vecs) {
    Vec<ndim> result;
    for (int i=0; i<vecs.size(); ++i)
        result += vecs[i];
    result.scaleInPlace(1.0/vecs.size());
    return result;
}

inline const Vec<3> sphereBaryCenter(const std::vector<Vec<3>>& vecs, const scalar_type radius = 1.0) {
    Vec<3> result;
    for (int i=0; i<vecs.size(); ++i)
        result += vecs[i];
    result.scaleInPlace(1.0/vecs.size());
    result.normalizeInPlace();
    result.scaleInPlace(radius);
    return result;
}

inline scalar_type triArea(const Vec<3>& vecA, const Vec<3>& vecB, const Vec<3>& vecC) {
    const Vec<3> s1 = vecB - vecA;
    const Vec<3> s2 = vecC - vecA;
    return 0.5*s1.crossProd(s2).mag();
}

inline scalar_type triArea(const Vec<2>& vecA, const Vec<2>& vecB, const Vec<2>& vecC) {
    const Vec<2> s1 = vecB - vecA;
    const Vec<2> s2 = vecC - vecA;
    return 0.5*s1.crossProdComp3(s2);
}

template <int ndim> scalar_type triArea(const std::vector<Vec<ndim>>& vecs) {
    return triArea(vecs[0], vecs[1], vecs[2]);
}

scalar_type sphereTriArea(const Vec<3>& a, const Vec<3>& b, const Vec<3>& c, const scalar_type radius = 1.0);

inline scalar_type sphereTriArea(const std::vector<Vec<3>>& vecs) {
    return sphereTriArea(vecs[0], vecs[1], vecs[2]);
}

scalar_type atan4(const scalar_type y, const scalar_type x);
// template <int ndim> std::ostream& operator << (std::ostream& os, const Vec<ndim>& vec);

std::ostream& operator << (std::ostream& os, const Vec<2>& vec);
std::ostream& operator << (std::ostream& os, const Vec<3>& vec);

}
#endif
