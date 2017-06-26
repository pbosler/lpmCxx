#include "LpmXyzVector.h"

namespace Lpm {

XyzVector operator + ( const XyzVector& vecA, const XyzVector& vecB ) {
	return XyzVector( vecA.x + vecB.x, vecA.y + vecB.y, vecA.z + vecB.z); }

XyzVector operator - ( const XyzVector& vecA, const XyzVector& vecB ) {
	return XyzVector( vecA.x - vecB.x, vecA.y - vecB.y, vecA.z - vecB.z); }
	
XyzVector operator * ( const XyzVector& vecA, const XyzVector& vecB ) {
	return XyzVector( vecA.x * vecB.x, vecA.y * vecB.y, vecA.z * vecB.z); }

std::ostream& operator << ( std::ostream& os, const XyzVector& vec ) {
	os << "(" << vec.x << ", " << vec.y << ", " << vec.z << ")" ; return os; }
	
scalar_type atan4( const scalar_type y, const scalar_type x ) {
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

bool operator == ( const XyzVector& vecA, const XyzVector& vecB ){
	return ( ( std::abs( vecA.x - vecB.x ) < ZERO_TOL ) && ( std::abs(vecA.y - vecB.y) < ZERO_TOL)) && ( std::abs(vecA.z - vecB.z) < ZERO_TOL ) ;
}

XyzVector midpoint( const XyzVector& vecA, const XyzVector& vecB ) {
	return XyzVector( 0.5 * (vecA.x + vecB.x), 0.5 * (vecA.y + vecB.y), 0.5 * (vecA.z + vecB.z) );
}

XyzVector centroid( const std::vector<XyzVector>& vecs ) {
	XyzVector cntd(0.0, 0.0, 0.0);
	for ( int i = 0; i < vecs.size(); ++i)
		cntd += vecs[i];
	cntd.scale( 1.0 / vecs.size() );
	return cntd;
}

scalar_type distance( const XyzVector& vecA, const XyzVector& vecB){
	return std::sqrt( (vecB.x - vecA.x) * (vecB.x - vecA.x) + (vecB.y - vecA.y) * ( vecB.y - vecA.y) + (vecB.z - vecA.z) * (vecB.z - vecA.z));
}

scalar_type triArea( const XyzVector& vecA, const XyzVector& vecB, const XyzVector& vecC){
	XyzVector diff1 = vecB - vecA;
	XyzVector diff2 = vecC - vecA;
	return 0.5 * diff1.crossProduct(diff2).magnitude();
}

scalar_type sphereDistance( const XyzVector& vecA, const XyzVector& vecB, const scalar_type radius ){
	const XyzVector cProd = vecA.crossProduct(vecB);
	const scalar_type dotProd = vecA.dotProduct(vecB);
	return std::atan2( cProd.magnitude(), dotProd ) * radius;
}

scalar_type sphereTriArea( const XyzVector& vecA, const XyzVector& vecB, const XyzVector& vecC, const scalar_type radius){
	const scalar_type side1 = sphereDistance(vecA, vecB);
	const scalar_type side2 = sphereDistance(vecB, vecC);
	const scalar_type side3 = sphereDistance(vecC, vecA);
	const scalar_type halfPerim = 0.5 * ( side1 + side2 + side3 );
	const scalar_type zz = std::tan( 0.5 * halfPerim ) * std::tan( 0.5 * (halfPerim - side1) ) *
		std::tan( 0.5 * ( halfPerim - side2 ) ) * std::tan( 0.5 * ( halfPerim - side3 ) );
	return 4.0 * std::atan2( std::sqrt(zz), 1.0 ) * radius * radius;
}

XyzVector sphereCentroid( const std::vector<XyzVector>& vecs, const scalar_type radius) {
	XyzVector cntd(0.0, 0.0, 0.0);
	for ( int i = 0; i < vecs.size(); ++i)
		cntd += vecs[i];
	cntd.scale( 1.0 / vecs.size() );
	cntd.normalize();
	cntd.scale( radius );
	return cntd;
}

XyzVector sphereMidpoint( const XyzVector& vecA, const XyzVector& vecB, const scalar_type radius) {
	XyzVector midpt( 0.5 * (vecA.x + vecB.x), 0.5 * (vecA.y + vecB.y), 0.5 * (vecA.z + vecB.z ) );
	midpt.normalize();
	midpt.scale(radius);
	return midpt;
}

}