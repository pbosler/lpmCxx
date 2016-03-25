#ifndef __LPM_XyzVector__
#define __LPM_XyzVector__

#include <iostream>
#include <cmath>
#include <vector>
#include "GlobalConstants.h"

template < typename scalarType > class XyzVector
{
	public :
		scalarType x;
		scalarType y;
		scalarType z;
		
		XyzVector( const scalarType nx = 0.0, const scalarType ny = 0.0, const scalarType nz = 0.0 ) :
			x(nx), y(ny), z(nz) {};
		
		~XyzVector() {};
		
		XyzVector operator += (const XyzVector& other) 
			{ x += other.x;
			  y += other.y;
			  z += other.z;
			  return *this;	}
		XyzVector operator += (XyzVector& other)
			{ x += other.x;
			  y += other.y;
			  z += other.z;
			  return *this;	}
		XyzVector operator -= (const XyzVector& other) 
			{ x -= other.x;
			  y -= other.y;
			  z -= other.z;
			  return *this;	}
		XyzVector operator -= (XyzVector& other)
			{ x -= other.x;
			  y -= other.y;
			  z -= other.z;
			  return *this;	}
		
		scalarType magnitude() const { return std::sqrt( x*x + y*y + z*z); }
		
		scalarType dotProduct( const XyzVector& other) const { return x*other.x + y*other.y + z*other.z	; }
		
		XyzVector crossProduct( const XyzVector& other) const 
			{ return XyzVector( y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x); }
		
		void scale( const scalarType multiplier ){ x *= multiplier; y *= multiplier; z *= multiplier; }
		
		void normalize() { const scalarType norm = magnitude(); x /= norm; y /= norm; z /= norm; }
			
};

template <typename scalarType> XyzVector<scalarType> operator + ( const XyzVector<scalarType>& vecA, const XyzVector<scalarType>& vecB ) {
	return XyzVector<scalarType>( vecA.x + vecB.x, vecA.y + vecB.y, vecA.z + vecB.z); }

template <typename scalarType> XyzVector<scalarType> operator - ( const XyzVector<scalarType>& vecA, const XyzVector<scalarType>& vecB ) {
	return XyzVector<scalarType>( vecA.x - vecB.x, vecA.y - vecB.y, vecA.z - vecB.z); }
	
template <typename scalarType> XyzVector<scalarType> operator * ( const XyzVector<scalarType>& vecA, const XyzVector<scalarType>& vecB ) {
	return XyzVector<scalarType>( vecA.x * vecB.x, vecA.y * vecB.y, vecA.z * vecB.z); }

template <typename scalarType> std::ostream& operator << ( std::ostream& os, const XyzVector<scalarType>& vec ) {
	os << "(" << vec.x << ", " << vec.y << ", " << vec.z << ")\n" ; return os; }
	
template <typename scalarType> scalarType atan4( const scalarType y, const scalarType x ) {
	static const scalarType PI = GlobalConstants::Instance()->Pi(); 
	scalarType result = 0.0;
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
		scalarType theta = std::atan2( std::abs(y), std::abs(x) );
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

template <typename scalarType> bool operator == ( const XyzVector<scalarType>& vecA, const XyzVector<scalarType>& vecB ){
	static const scalarType	ZERO_TOL = GlobalConstants::Instance()->ZeroTol();
	return ( ( std::abs( vecA.x - vecB.x ) < ZERO_TOL ) && ( std::abs(vecA.y - vecB.y) < ZERO_TOL)) && ( std::abs(vecA.z - vecB.z) < ZERO_TOL ) ;
}

template <typename scalarType> XyzVector<scalarType> midpoint( const XyzVector<scalarType>& vecA, const XyzVector<scalarType>& vecB ) {
	return XyzVector<scalarType>( 0.5 * (vecA.x + vecB.x), 0.5 * (vecA.y + vecB.y), 0.5 * (vecA.z + vecB.z) );
}

template <typename scalarType> XyzVector<scalarType> centroid( const std::vector<XyzVector<scalarType> > vecs ) {
	XyzVector<scalarType> cntd(0.0, 0.0, 0.0);
	for ( int i = 0; i < vecs.size(); ++i)
		cntd += vecs[i];
	cntd.scale( 1.0 / vecs.size() );
	return cntd;
}

template <typename scalarType> scalarType distance( const XyzVector<scalarType>& vecA, const XyzVector<scalarType>& vecB){
	return std::sqrt( (vecB.x - vecA.x) * (vecB.x - vecA.x) + (vecB.y - vecA.y) * ( vecB.y - vecA.y) + (vecB.z - vecA.z) * (vecB.z - vecA.z));
}

//template <typename scalarType> scalarType triArea( const std::vector<XyzVector<scalarType> > vecs ) { }

template <typename scalarType> scalarType triArea( const XyzVector<scalarType>& vecA, const XyzVector<scalarType>& vecB, const XyzVector<scalarType>& vecC){
	XyzVector<scalarType> diff1 = vecB - vecA;
	XyzVector<scalarType> diff2 = vecC - vecA;
	return 0.5 * diff1.crossProduct(diff2).magnitude();
}

template <typename scalarType> scalarType sphereDistance( const XyzVector<scalarType>& vecA, const XyzVector<scalarType>& vecB, 
	const scalarType radius = 1.0 ){
	XyzVector<scalarType> cProd = vecA.crossProduct(vecB);
	const scalarType cProdNorm = cProd.magnitude();
	const scalarType dotProd = vecA.dotProduct(vecB);
	return std::atan2( cProdNorm, dotProd ) * radius;
}

template <typename scalarType> scalarType sphereTriArea( const XyzVector<scalarType>& vecA, 
	const XyzVector<scalarType>& vecB, const XyzVector<scalarType>& vecC, const scalarType radius = 1.0){
	const scalarType side1 = sphereDistance(vecA, vecB);
	const scalarType side2 = sphereDistance(vecB, vecC);
	const scalarType side3 = sphereDistance(vecC, vecA);
	const scalarType halfPerim = 0.5 * ( side1 + side2 + side3 );
	const scalarType zz = std::tan( 0.5 * halfPerim ) * std::tan( 0.5 * (halfPerim - side1) ) *
		std::tan( 0.5 * ( halfPerim - side2 ) ) * std::tan( 0.5 * ( halfPerim - side3 ) );
	return 4.0 * std::atan2( zz, 1.0 ) * radius * radius;
}

template <typename scalarType> XyzVector<scalarType> sphereCentroid( const std::vector<XyzVector<scalarType> > vecs, 
	const scalarType radius = 1.0 ) {
	XyzVector<scalarType> cntd(0.0, 0.0, 0.0);
	for ( int i = 0; i < vecs.size(); ++i)
		cntd += vecs[i];
	cntd.scale( 1.0 / vecs.size() );
	cntd.normalize();
	cntd.scale( radius );
	return cntd;
}

template <typename scalarType> XyzVector<scalarType> sphereMidpoint( const XyzVector<scalarType>& vecA, 
	const XyzVector<scalarType>& vecB, const scalarType radius = 1.0 ) {
	XyzVector<scalarType> midpt( 0.5 * (vecA.x + vecB.x), 0.5 * (vecA.y + vecB.y), 0.5 * (vecA.z + vecB.z ) );
	midpt.normalize();
	midpt.scale(radius);
	return midpt;
}

#endif
