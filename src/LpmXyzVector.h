#ifndef __LPM_XyzVector__
#define __LPM_XyzVector__

#include "LpmTypeDefs.h"
#include "LpmConfig.h"
#include <iostream>
#include <cmath>
#include <vector>

namespace Lpm {

class XyzVector
{
	public :
		scalar_type x;
		scalar_type y;
		scalar_type z;
		
		XyzVector( const scalar_type nx = 0.0, const scalar_type ny = 0.0, const scalar_type nz = 0.0 ) :
			x(nx), y(ny), z(nz) {};
		
		~XyzVector() {};
		
		inline XyzVector operator += (const XyzVector& other) 
			{ x += other.x;
			  y += other.y;
			  z += other.z;
			  return *this;	}
		inline XyzVector operator += (XyzVector& other)
			{ x += other.x;
			  y += other.y;
			  z += other.z;
			  return *this;	}
		inline XyzVector operator -= (const XyzVector& other) 
			{ x -= other.x;
			  y -= other.y;
			  z -= other.z;
			  return *this;	}
		inline XyzVector operator -= (XyzVector& other)
			{ x -= other.x;
			  y -= other.y;
			  z -= other.z;
			  return *this;	}
		
		inline scalar_type magnitude() const { return std::sqrt( x*x + y*y + z*z); }
		
		inline scalar_type dotProduct( const XyzVector& other) const { return x*other.x + y*other.y + z*other.z	; }
		
		inline XyzVector crossProduct( const XyzVector& other) const 
			{ return XyzVector( y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x); }
		
		inline void scale( const scalar_type multiplier ){ x *= multiplier; y *= multiplier; z *= multiplier; }
		
		inline void normalize() { const scalar_type norm = magnitude(); x /= norm; y /= norm; z /= norm; }
			
};


}
#endif
