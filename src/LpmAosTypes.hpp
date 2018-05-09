#ifndef LPM_AOS_TYPES_HPP
#define LPM_AOS_TYPES_HPP

#include <cmath>
#include <iostream>
#include "LpmConfig.h"
#include "LpmTypeDefs.h"

namespace Lpm {

scalar_type atan4(const scalar_type y, const scalar_type x);

template<int ndim=3> struct Vec {
    scalar_type x[ndim];
    
    Vec(scalar_type* xx) {
        for (int i=0; i<ndim; ++i) 
            x[i] = xx[i];
    } 
    
    Vec<ndim> operator += (const Vec<ndim>& other) {
        for (int i=0; i<ndim; ++i)
            this->x[i] += other.x[i];
    }
    
    Vec<ndim> operator -= (const Vec<ndim>& other) {
        for (int i=0; i<ndim; ++i)
            this->x[i] -= other.x[i];
    }
    
//     inline Vec3 operator *= (const Vec3& other) {
//         x *= other.x;
//         y *= other.y;
//         z *= other.z;
//     }
//     
//     inline Vec3 operator /= (const Vec3& other) {
//         x /= other.x;
//         y /= other.y;
//         z /= other.z;
//     }
//     
//     inline Vec3 plus(const Vec3& other) const {return Vec3(x+other.x, y+other.y, z+other.z);}
//     inline Vec3 minus(const Vec3& other) const {return Vec3(x-other.x, y-other.y, z-other.z);}
//     inline Vec3 times(const Vec3& other) const {return Vec3(x*other.x, y*other.y, z*other.z);}
//     inline Vec3 divide(const Vec3& other) const {return Vec3(x/other.x, y/other.y, z/other.z);}
//     
//     inline scalar_type mag() const {return std::sqrt(x*x + y*y + z*z);}
//     inline scalar_type magSq() const {return x*x + y*y + z*z;}
//     
//     inline scalar_type dotProd(const Vec3& other) const {return x*other.x + y*other.y + z*other.z;}
//     
//     inline Vec3 crossProd(const Vec3& other) const 
//         {return Vec3(y*other.z - z*other.y, z*other.x - x*other.z, x*other.y - y*other.x);}
//     
//     inline void scaleInPlace(const scalar_type mult) {x*=mult; y*=mult; z*=mult;}
//     
//     inline Vec3 scale(const scalar_type mult) const 
//         {return Vec3(x*mult, y*mult, z*mult);}
//     
//     inline void normalizeInPlace() { 
//         const scalar_type len = mag();
//         x /= len;
//         y /= len;
//         z /= len;}
//     
//     inline Vec3 normalize() const {
//         const scalar_type len = mag();
//         return Vec3(x/len, y/len, z/len):}
//         
//     inline scalar_type longitude() const {return atan4(y,x);}
//     inline scalar_type latitude() const {return std::atan2(z / std::sqrt(x*x + y*y));}
//     
//     inline Vec3 midpoint(const Vec3& other) const {return (this->plus(other)).scale(0.5);}
//     inline Vec3 sphereMidpoint(const Vec3& other, const scalar_type radius = 1.0) const {
//         return midpoint(other).scale(1.0/radius);
//     }
//     
//     inline scalar_type dist(const Vec3& other) const {
//         return this->minus(other).mag();
//     }
//     
//     inline scalar_type sphereDist(const Vec3& other, const scalar_type radius = 1.0) const {
//         const Vec3 cp = this->crossProd(other);
//         const scalar_type dp = this->dotProd(other);
//         return std::atan2(cp.mag(), dotProd) * radius;
//     }
};

// inline bool operator == (const Vec3& a, const Vec3& b) {
//     return ((b-a).mag() < ZERO_TOL);
// }
// 
// inline Vec3 operator + (const Vec3& a, const Vec3& b) {
//     return a.plus(b);
// }
// 
// inline Vec3 operator - (const Vec3& a, const Vec3& b) {
//     return a.minus(b);
// }
// 
// inline Vec3 operator * (const Vec3& a, const Vec3& b) {
//     return a.times(b);
// }
// 
// inline Vec3 operator / (const Vec3& a, const Vec3& b) {
//     return a.divide(b);
// }
// 
// inline Vec3 midpoint(const Vec3& a, const Vec3& b) {
//     return a.midpoint(b);
// }
// 
// inline Vec3 sphereMidpoint(const Vec3& a, const Vec3& b, const scalar_type radius = 1.0) {
//     return a.sphereMidpoint(b,radius);
// }
// 
// inline Vec3 baryCen(const std::vector<Vec3>& vecs) {
//     Vec3 result;
//     for (index_type i=0; i<vecs.size(); ++i)
//         result += vecs[i]
//     return result.scale(1.0/vecs.size());
// }
// 
// inline Vec3 sphereBaryCen(const std::vector<Vec3>& vecs, const scalar_type radius = 1.0) {
//     Vec3 result = baryCen(vecs);
//     return result.scale(1.0/radius);
// }
// 
// inline scalar_type dist(const Vec3& a, const Vec3& b) {return a.dist(b);}
// 
// inline scalar_type sphereDist(const Vec3&a, const Vec3& b, const scalar_type radius = 1.0) {
//     return a.sphereDist(b, radius);}
// 
// inline scalar_type triArea(const Vec3& a, const Vec3& b, const Vec3& c) {
//     Vec3 d1 = b-a;
//     Vec3 d2 = c-a;
//     return 0.5*(d1.crossProd(d2).mag());
// }
// 
// inline scalar_type sphereTriArea(const Vec3& a, const Vec3& b, const Vec3& c, const scalar_type radius = 1.0) {
//     const scalar_type s1 = sphereDist(a,b);
//     const scalar_type s2 = sphereDist(b,c);
//     const scalar_type s3 = sphereDist(c,a);
//     const scalar_type halfPerim = 0.5*(s1 + s2 + s3);
//     const scalar_type zz = std::tan(0.5*halfPerim) * std::tan(0.5*(halfPerim-s1)) * std::tan(0.5*(halfPerim-s2)) *
//         std::tan(0.5*(halfPerim-s3));
//     return 4.0 * std::atan(std::sqrt(zz)) * radius*radius;
// }
// 
// scalar_type atan4(const scalar_type y, const scalar_type x) {
//     scalar_type result = 0.0;
// 	if ( x == 0.0 )
// 	{
// 		if ( y > 0.0 ) 
// 			result = 0.5 * PI;
// 		else if ( y < 0.0 )
// 			result = 1.5 * PI;
// 		else if ( y == 0.0 )
// 			result = 0.0;
// 	}
// 	else if ( y == 0 )
// 	{
// 		if ( x > 0.0 )
// 			result = 0.0;
// 		else if ( x < 0.0 )
// 			result = PI;
// 	}
// 	else
// 	{
// 		scalar_type theta = std::atan2( std::abs(y), std::abs(x) );
// 		if ( x > 0.0 && y > 0.0 )
// 			result = theta;
// 		else if ( x < 0.0 && y > 0.0 ) 
// 			result = PI - theta;
// 		else if ( x < 0.0 && y < 0.0 ) 
// 			result = PI + theta;
// 		else if ( x > 0.0 && y < 0.0 )
// 			result = 2.0 * PI - theta;
// 	}
// 	return result;
// }
// 
// }

std::ostream& operator << (std::ostream& os, const Vec<2>& vec) {
    os << "(" << vec.x[0] << ", " << vec.x[1] << ")";
    return os;
}

std::ostream& operator << (std::ostream& os, const Vec<3>& vec) {
        os << "(" << vec.x[0] << ", " << vec.x[1] << ", " << vec.x[2] << ")";
}

#endif