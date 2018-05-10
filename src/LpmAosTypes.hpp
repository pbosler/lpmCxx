#ifndef LPM_AOS_TYPES_HPP
#define LPM_AOS_TYPES_HPP

#include <cmath>
#include <iostream>
//#include <vector>
#include <exception>
#include "LpmConfig.h"
#include "LpmTypeDefs.h"

namespace Lpm {

scalar_type atan4(const scalar_type y, const scalar_type x);

template<int ndim=3> struct Vec {
    scalar_type x[ndim];

    Vec() {
        for (int i=0; i<ndim; ++i)
            this->x[i] = 0.0;
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

    Vec(const std::vector<scalar_type>& xx) {
        for (int i=0; i<ndim; ++i)
            this->x[i] = xx[i];
    }

    Vec(const scalar_type* xx) {
        for (int i=0; i<ndim; ++i)
            this->x[i] = xx[i];
    }

    Vec<ndim>& operator += (const Vec<ndim>& other) {
        for (int i=0; i<ndim; ++i)
            this->x[i] += other.x[i];
        return *this;
    }

    Vec<ndim>& operator -= (const Vec<ndim>& other) {
        for (int i=0; i<ndim; ++i)
            this->x[i] -= other.x[i];
        return *this;
    }

    Vec<ndim>& operator *= (const Vec<ndim>& other) {
        for (int i=0; i<ndim; ++i)
            this->x[i] *= other.x[i];
        return *this;
    }

    Vec<ndim>& operator /= (const Vec<ndim>& other) {
        for (int i=0; i<ndim; ++i)
            this->x[i] /= other.x[i];
        return *this;
    }

    const Vec<ndim> operator + (const Vec<ndim>& other) const {
        return Vec<ndim>(*this) += other;
    }

    const Vec<ndim> operator - (const Vec<ndim>& other) const {
        return Vec<ndim>(*this) -= other;
    }

    const Vec<ndim> operator * (const Vec<ndim>& other) const {
        return Vec<ndim>(*this) *= other;
    }

    const Vec<ndim> operator / (const Vec<ndim>& other) const {
        return Vec<ndim>(*this) /= other;
    }

    scalar_type magSq() const {
        scalar_type sumsq(0.0);
        for (int i=0; i<ndim; ++i)
            sumsq += this->x[i]*this->x[i];
        return sumsq;
    }

    scalar_type mag() const {
        return std::sqrt(this->magSq());
    }

    scalar_type dotProd(const Vec<ndim>& other) const {
        scalar_type dp(0.0);
        for (int i=0; i<ndim; ++i)
            dp += this->x[i]*other.x[i];
        return dp;
    }

    scalar_type crossProd(const Vec<2>& other) const {
        return this->x[0]*other.x[1] - this->x[1]*other.x[0];
    }

    const Vec<3> crossProd(const Vec<3>& other) const {
        const scalar_type cp[3] = {this->x[1]*other.x[2] - this->x[2]*other.x[1],
                                   this->x[2]*other.x[0] - this->x[0]*other.x[2],
                                   this->x[0]*other.x[1] - this->x[1]*other.x[0]};
        return Vec<3>(cp);
    }

    void scaleInPlace(const scalar_type mult) {
        for (int i=0; i<ndim; ++i)
            this->x[i] *= mult;
    }

    const Vec<ndim> scale(const scalar_type mult) const {
        scalar_type sm[ndim];
        for (int i=0; i<ndim; ++i)
            sm[i] = this->x[i]*mult;
        return Vec<ndim>(sm);
    }

    void normalizeInPlace() {
        const scalar_type len = this->mag();
        this->scaleInPlace(1.0/len);
    }

    const Vec<ndim> normalize() const {
        const scalar_type len = this->mag();
        return this->scale(1.0/len);
    }

    scalar_type longitude() const {return atan4(this->x[1], this->x[0]);}

    scalar_type latitude() const {
        const scalar_type xy2 = this->x[0]*this->x[0] + this->x[1]*this->x[1];
        return std::atan2(this->x[2] , std::sqrt(xy2));
    }

    const Vec<ndim> midpoint(const Vec<ndim>& other) const {
        Vec<ndim> result = *this + other;
        result.scaleInPlace(0.5);
        return result;
    }

    const Vec<3> sphereMidpoint(const Vec<3>& other, const scalar_type radius = 1.0) const {
        Vec<3> result = this->midpoint(other);
        result.normalizeInPlace();
        result.scaleInPlace(radius);
        return result;
    }

    scalar_type dist(const Vec<ndim>& other) const {
        return (*this - other).mag();
    }

    scalar_type sphereDist(const Vec<3>& other, const scalar_type radius=1.0) const {
        const Vec<3> cp = this->crossProd(other);
        const scalar_type dp = this->dotProd(other);
        return std::atan2(cp.mag(), dp) * radius;
    }

    const bool operator == (const Vec<ndim>& other) const {
        return (*this - other).mag() < ZERO_TOL;
    }
};

template <int ndim> const Vec<ndim> baryCenter(const std::vector<Vec<ndim>>& vecs) {
    Vec<ndim> result;
    for (int i=0; i<vecs.size(); ++i)
        result += vecs[i];
    result.scaleInPlace(1.0/vecs.size());
    return result;
}

const Vec<3> sphereBaryCenter(const std::vector<Vec<3>>& vecs, const scalar_type radius = 1.0) {
    Vec<3> result;
    for (int i=0; i<vecs.size(); ++i)
        result += vecs[i];
    result.scaleInPlace(1.0/vecs.size());
    result.normalizeInPlace();
    result.scaleInPlace(radius);
    return result;
}

scalar_type triArea(const Vec<3>& vecA, const Vec<3>& vecB, const Vec<3>& vecC) {
    const Vec<3> s1 = vecB - vecA;
    const Vec<3> s2 = vecC - vecA;
    return 0.5*s1.crossProd(s2).mag();
}

scalar_type triArea(const Vec<2>& vecA, const Vec<2>& vecB, const Vec<2>& vecC) {
    const Vec<2> s1 = vecB - vecA;
    const Vec<2> s2 = vecC - vecA;
    return 0.5*s1.crossProd(s2);
}

template <int ndim> scalar_type triArea(const std::vector<Vec<ndim>>& vecs) {
    return triArea(vecs[0], vecs[1], vecs[2]);
}

scalar_type sphereTriArea(const Vec<3>& a, const Vec<3>& b, const Vec<3>& c, const scalar_type radius = 1.0) {
    const scalar_type s1 = a.sphereDist(b, radius);
    const scalar_type s2 = b.sphereDist(c, radius);
    const scalar_type s3 = c.sphereDist(a, radius);
    const scalar_type halfPerim = 0.5*(s1 + s2 + s3);
    const scalar_type zz = std::tan(0.5*halfPerim) * std::tan(0.5*(halfPerim-s1)) * std::tan(0.5*(halfPerim-s2)) *
        std::tan(0.5*(halfPerim-s3));
    return 4.0 * std::atan(std::sqrt(zz)) * radius*radius;
}

scalar_type sphereTriArea(const std::vector<Vec<3>>& vecs) {
    return sphereTriArea(vecs[0], vecs[1], vecs[2]);
}

scalar_type atan4(const scalar_type y, const scalar_type x) {
    scalar_type result = 0.0;
	if ( x == 0.0 )
	{
		if ( y > 0.0 )
			result = 0.5 * PI;
		else if ( y < 0.0 )
			result = 1.5 * PI;
		else if ( y == 0.0 )
			result = 0.0;
	}
	else if ( y == 0 )
	{
		if ( x > 0.0 )
			result = 0.0;
		else if ( x < 0.0 )
			result = PI;
	}
	else
	{
		scalar_type theta = std::atan2( std::abs(y), std::abs(x) );
		if ( x > 0.0 && y > 0.0 )
			result = theta;
		else if ( x < 0.0 && y > 0.0 )
			result = PI - theta;
		else if ( x < 0.0 && y < 0.0 )
			result = PI + theta;
		else if ( x > 0.0 && y < 0.0 )
			result = 2.0 * PI - theta;
	}
	return result;
}

std::ostream& operator << (std::ostream& os, const Vec<2>& vec) {
    os << "(" << vec.x[0] << ", " << vec.x[1] << ")";
    return os;
}

std::ostream& operator << (std::ostream& os, const Vec<3>& vec) {
        os << "(" << vec.x[0] << ", " << vec.x[1] << ", " << vec.x[2] << ")";
}

}
#endif
