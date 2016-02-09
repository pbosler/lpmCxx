#include "xyzVector.h"
#include <cmath>

/** @file xyzVector.cpp
	@author Peter Bosler, Sandia National Laboratories
	@brief xyzVector class implementation
*/

/** @brief overloaded addition operator. Component-by-component addition.
 */
xyzVector xyzVector::operator + ( xyzVector& other)
{
	return xyzVector( x + other.x, y + other.y, z + other.z);
}

/** @brief overloaded subtraction operator. Component-by-component subtraction.
 */
xyzVector xyzVector::operator - ( xyzVector& other)
{
	return xyzVector( x - other.x, y - other.y, z - other.z);
}

/** @brief overloaded addition-assignment operator. Component-by-component addition.
 */
xyzVector xyzVector::operator += ( const xyzVector& other)
{
	x += other.x;
	y += other.y;
	z += other.z;
	return *this;
}

/** @brief overloaded addition-assignment operator. Component-by-component addition.
 */
xyzVector xyzVector::operator += ( xyzVector& other)
{
	x += other.x;
	y += other.y;
	z += other.z;
	return *this;
}
    
/** @brief overloaded subtraction-assignment operator. Component-by-component subtraction.
 */
xyzVector xyzVector::operator -= (const xyzVector& other )
{
	x -= other.x;
	y -= other.y;
	z -= other.z;
	return *this;
}

/** @brief overloaded subtraction-assignment operator. Component-by-component subtraction.
 */
xyzVector xyzVector::operator -= (xyzVector& other )
{
	x -= other.x;
	y -= other.y;
	z -= other.z;
	return *this;
}
    
/** @brief scalar multiplication of a position vector
	@param multiplier
 */
void xyzVector::scale( const double multiplier)
{
	x *= multiplier;
	y *= multiplier;
	z *= multiplier;
};
    
/** @brief returns the magnitude of a position vector
 */
double xyzVector::magnitude() const
{
	return std::sqrt( x*x + y*y + z*z);
}
    
/** @brief vector inner product
	@param other
 */
double xyzVector::dotProduct( const xyzVector& other) const
{
	return x*other.x + y*other.y + z*other.z;
};
    
/** @brief vector cross product
 @param other
 */
xyzVector xyzVector::crossProduct( const xyzVector& other) const
{
	return xyzVector( y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x );
}
    
/** @brief normalizes a vector so that its length is 1.
 */
void xyzVector::normalize()
{
	const double mag = magnitude();
	x /= mag;
	y /= mag;
	z /= mag;
};

/** @brief overloaded addition operator
	@param vecA
	@param vecB
 */
xyzVector operator + ( const xyzVector& vecA, const xyzVector& vecB)
{
    return xyzVector( vecA.x + vecB.x, vecA.y + vecB.y, vecA.z + vecB.z);
};

/** @brief overloaded subtraction operator
	@param vecA
	@param vecB
 */
xyzVector operator - ( const xyzVector& vecA, const xyzVector& vecB)
{
    return xyzVector( vecA.x - vecB.x, vecA.y - vecB.y, vecA.z - vecB.z);
}

/** @brief basic console output
@param os
@param vec
 */
std::ostream& operator<<( std::ostream& os, const xyzVector& vec )
{
    os << "( " << vec.x << ", " << vec.y << ", " << vec.z << " )" ;
    return os;
}

/** @brief the inverse tangent function, with output in range [0,2*pi) instead of (-pi, pi) as in std::atan2.
 @param y
 @param x
 */
 double atan4( const double y, const double x)
{
    GlobalConstants* physConsts = GlobalConstants::Instance();
    double result = 0.0;
    if ( x == 0.0)
    {
        if ( y > 0.0)
            result = 0.5 * physConsts->Pi();
        else if ( y < 0.0 )
            result = 1.5 * physConsts->Pi();
        else if ( y == 0.0 )
            result = 0.0;
    }
    else if ( y == 0 )
    {
        if ( x > 0.0 )
            result = 0.0;
        else if ( x < 0.0 )
            result = physConsts->Pi();
    }
    else
    {
        double theta = std::atan2( std::abs(y), std::abs(x));
        if ( x > 0.0 && y > 0.0 )
            result = theta;
        else if ( x < 0.0 && y > 0.0 )
            result = physConsts->Pi() - theta;
        else if ( x < 0.0 && y < 0.0 )
            result = physConsts->Pi() + theta;
        else if ( x > 0.0 && y < 0.0 )
            result = 2.0 * physConsts->Pi() - theta;
    }
    return result;
}

/** equality comparison operator for position vectors
 @param vecA
 @param vecB
 */
 bool operator == ( const xyzVector& vecA, const xyzVector& vecB)
{
    return ( (vecA.x == vecB.x && vecA.y == vecB.y) && ( vecA.z == vecB.z));
}

/** @brief returns the Euclidean distance between two vectors
    @param vecA
    @param vecB
 */
double cartDistance( const xyzVector& vecA, const xyzVector& vecB)
{
    return std::sqrt( (vecA.x - vecB.x)*(vecA.x - vecB.x) +
                      (vecA.y - vecB.y)*(vecA.y - vecB.y) +
                      (vecA.z - vecB.z)*(vecA.z - vecB.z));
}

/** @brief returns the midpoint of a chord connecting two vectors
 @param vecA
 @param vecB
 */
xyzVector cartMidpoint( const xyzVector& vecA, const xyzVector& vecB)
{
    return xyzVector( 0.5 * ( vecA.x + vecB.x), 0.5 * (vecA.y + vecB.y), 0.5 * (vecA.z + vecB.z));
}

/** returns the vector centroid of a polygon whose vertices are given as a std::vector of xyzVectors.
 @param vecs
 */
xyzVector cartCentroid( const std::vector<xyzVector> vecs)
{
    xyzVector cntd(0.0, 0.0);
    for ( int i = 0; i < vecs.size(); ++i)
    {
        cntd += vecs[i];
    }
    cntd.scale(1.0/vecs.size());
    return cntd;
}

/** @brief returns the area of a the triangle whose vertices are given by the vectors vecA, vecB, and vecC.
 @param vecA
 @param vecB
 @param vecC
*/
double cartTriArea( const xyzVector& vecA, const xyzVector& vecB, const xyzVector& vecC)
{
    xyzVector diff1 = vecB - vecA;
    xyzVector diff2 = vecC - vecA;
    return  0.5 * diff1.crossProduct(diff2).magnitude();
}



/** @brief returns the great-circle distance between two vectors on a sphere with radius 
 defined by @ref GlobalConstants::_earthRadius.
 @param vecA
 @param vecB
 */
double sphereDistance( const xyzVector& vecA, const xyzVector& vecB, const double radius)
{
    return atan2( (vecA.crossProduct(vecB)).magnitude(), vecA.dotProduct(vecB)) * radius;
}

/** @brief returns the central angle between two vectors on a sphere.
 @param vecA
 @param vecB
 */
double sphereAngle( const xyzVector& vecA, const xyzVector& vecB)
{
    return atan2( vecA.crossProduct(vecB).magnitude(), vecA.dotProduct(vecB));
}

/** @brief returns the midpoint of a great-circle arc connecting two vectors on the sphere.
 @warning does not check input to make sure they have norm = sphere radius.
 @param vecA
 @param vecB
 */
xyzVector sphereMidpoint( const xyzVector& vecA, const xyzVector& vecB, const double radius )
{
    xyzVector result( 0.5 * ( vecA.x + vecB.x), 0.5 * (vecA.y + vecB.y), 0.5 * (vecA.z + vecB.z));
    result.normalize();
    result.scale( radius );
    return result;
}

/** returns the centroid of a spherical polygon whose vertices are contained in a std::vector of xyzVectors.
 @param vecs
 */
 xyzVector sphereCentroid( const std::vector<xyzVector> vecs, const double radius )
{
    xyzVector result = cartCentroid(vecs);
    result.normalize();
    result.scale( radius );
    return result;
}

/** returns the area of a spherical triangle whose vertices are given by vecA, vecB, and vecC
 @warning does not check input to make sure they have norm = sphere radius.
 @param vecA
 @param vecB
 @param vecC
 */
double sphereTriArea( const xyzVector& vecA, const xyzVector& vecB, const xyzVector& vecC, const double radius)
{
    double side1 = sphereAngle(vecA, vecB);
    double side2 = sphereAngle(vecB, vecC);
    double side3 = sphereAngle(vecC, vecA);
    
    double halfPerim = 0.5 * ( side1 + side2 + side3);
    
    double zz = std::tan( 0.5 * halfPerim ) * std::tan( 0.5 * (halfPerim - side1) ) *
    std::tan( 0.5 * (halfPerim - side2) ) * std::tan( 0.5 * (halfPerim - side3));
    
    return 4.0 * std::atan2( std::sqrt(zz), 1.0 ) * radius * radius;
}

/** @brief returns the longitude of an xyzVector
 */
double longitude( const xyzVector vec)
{
    return atan4( vec.y, vec.x);
}

/** @brief returns the latitude of an xyzVector
 */
double latitude( const xyzVector& vec)
{
    return std::atan2( vec.z, std::sqrt( vec.x * vec.x + vec.y * vec.y));
}

double distance( const xyzVector& vecA, const xyzVector& vecB, const xyzVector::geometryKind gk)
{ 
	if ( gk == xyzVector::SphericalGeometry )
		return sphereDistance( vecA, vecB ); 
	else
		return cartDistance(vecA, vecB);
}

xyzVector midpoint( const xyzVector& vecA, const xyzVector& vecB, const xyzVector::geometryKind gk)
{ 
	if ( gk == xyzVector::SphericalGeometry )
		return sphereMidpoint( vecA, vecB); 
	else
		return cartMidpoint(vecA, vecB);
}

xyzVector centroid( const std::vector<xyzVector> vecs, const xyzVector::geometryKind gk)
{ 
	if ( gk == xyzVector::SphericalGeometry )
		return sphereCentroid( vecs ); 
	else
		return cartCentroid( vecs );
}
double triArea( const xyzVector& vecA, const xyzVector& vecB, const xyzVector& vecC, const xyzVector::geometryKind gk)
{ 
	if ( gk == xyzVector::SphericalGeometry )
		return sphereTriArea( vecA, vecB, vecC ); 
	else
		return cartTriArea( vecA, vecB, vecC );
}
