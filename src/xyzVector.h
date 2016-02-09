//
//  xyzVector.h
//  LPM
//
//  Created by Peter Bosler on 11/2/14.
//  Copyright (c) 2014 Peter Bosler. All rights reserved.
//
/**
    @file xyzVector.h
    @author Peter Bosler, Sandia National Laboratories
    @brief header and implementation for xyzVector class.
*/

/**    
    @class xyzVector
    @brief A worker class for vector arithmetic in @f$ \mathbb{R}^3 @f$ ( and @f$\mathbb{R}^2  @f$ with z component = 0 ).
 
    Contains methods for Euclidean geometry as well as some spherical functions for a sphere with radius given by
    @ref GlobalConstants::EarthRadius().
    To use spherical geometry, initialize the instance with geomKind=geometryKind::SphericalGeometry. 
    This will set the functions `distance`, `centroid`, `midpoint`, and `triArea` to return the spherical functions
    instead of the default Euclidean functions.
*/
#ifndef __LPM__xyzVector__
#define __LPM__xyzVector__

#include <iostream>
#include <cmath>
#include <vector>
#include "GlobalConstants.h"

double atan4( const double y, const double x);

class xyzVector
{
public:
	enum geometryKind { EuclideanGeometry, SphericalGeometry};
	
	double x; ///< x coordinate
    double y; ///< y coordinate
    double z; ///< z coordinate
	
    /** @brief Default constructor
        @param nx
        @param ny
        @param nz
        @param gt set to SphericalGeometry to compute geodetic distance instead of Euclidean distance
     */
    xyzVector(  const double nx = 0.0, const double ny = 0.0, const double nz = 0.0) :
    			x(nx), y(ny), z(nz){};

    /// returns the longitude of an xzyVector
    inline double longitude() const { return atan4( y, x); }
    
    /// returns the latitude of an xyzVector
    inline double latitude() const { return std::atan2( z, std::sqrt( x*x + y*y)); }
    
    /** @brief overloaded addition operator. Component-by-component addition.
     */
    xyzVector operator + ( xyzVector& other);
    
    /** @brief overloaded subtraction operator. Component-by-component subtraction.
     */
    xyzVector operator - ( xyzVector& other);

    /** @brief overloaded addition-assignment operator. Component-by-component addition.
     */
    xyzVector operator += ( const xyzVector& other);
    
    /** @brief overloaded addition-assignment operator. Component-by-component addition.
     */
    xyzVector operator += ( xyzVector& other);
    
    /** @brief overloaded subtraction-assignment operator. Component-by-component subtraction.
     */
    xyzVector operator -= (const xyzVector& other );
    
    /** @brief overloaded subtraction-assignment operator. Component-by-component subtraction.
     */
    xyzVector operator -= ( xyzVector& other );
    
    /** @brief scalar multiplication of a position vector
        @param multiplier
     */
    void scale( const double multiplier);
    
    /** @brief returns the magnitude of a position vector
     */
    double magnitude() const;
    
    /** @brief vector inner product
        @param other
     */
    double dotProduct( const xyzVector& other) const;
    
    /** @brief vector cross product
     @param other
     */
	xyzVector crossProduct( const xyzVector& other) const;
    
    /** @brief normalizes a vector so that its length is 1.
     */
    void normalize();
};

/// basic console output
std::ostream& operator<<( std::ostream& os, const xyzVector& vec );
/** equality operator 
	@todo implement with GlobalConstants::ZeroTol() to compare real numbers
*/
bool operator == ( const xyzVector& vecA, const xyzVector& vecB);

/// Euclidean distance between two xyzVector objects given in Cartesian coordinates
double cartDistance( const xyzVector& vecA, const xyzVector& vecB);
/// Midpoint of two xyzVector objects given in Cartesian coordinates
xyzVector cartMidpoint( const xyzVector& vecA, const xyzVector& vecB);
/// Centroid of several xyzVector objects given in Cartesian coordinates
xyzVector cartCentroid( const std::vector<xyzVector> vecs);
/// Area of a planar triangle whose 3 vertices are xyzVector objects given in Cartesian coordinates
double cartTriArea( const xyzVector& vecA, const xyzVector& vecB, const xyzVector& vecC);

xyzVector crossProduct( const xyzVector& vecA, const xyzVector& vecB);

/** @brief Great-circle distance between two xyzVector objects given in Cartesian coordinates, 
	each with magnitude = GlobalConstants::EarthRadius()
*/
double sphereDistance( const xyzVector& vecA, const xyzVector& vecB, const double radius = 1.0 );
/** @brief Angular separation between two xyzVector objects given in Cartesian coordinates, 
	each with magnitude = GlobalConstants::EarthRadius()
*/
double sphereAngle( const xyzVector& vecA, const xyzVector& vecB);
/** @brief Great-circle distance between two xyzVector objects given in Cartesian coordinates
*/
xyzVector sphereMidpoint( const xyzVector& vecA, const xyzVector& vecB, const double radius = 1.0 );
/** @brief Centroid of several xyzVector objects given in Cartesian coordinates
*/
xyzVector sphereCentroid( const std::vector<xyzVector> vecs, const double radius = 1.0 );
/** @brief Area of a spherical triangle with 3 vertices defined by xyzVector objects given in Cartesian coordinates 
*/
double sphereTriArea( const xyzVector& vecA, const xyzVector& vecB, const xyzVector& vecC, const double radius = 1.0);
/// longitude of an xyzVector with magnitude = GlobalConstants::EarthRadius()
double longitude( const xyzVector vec);
/// latitude of an xyzVector with magnitude = GlobalConstants::EarthRadius()
double latitude( const xyzVector& vec);

double distance( const xyzVector& vecA, const xyzVector& vecB, const xyzVector::geometryKind gk);
xyzVector midpoint( const xyzVector& vecA, const xyzVector& vecB, const xyzVector::geometryKind gk);
xyzVector centroid( const std::vector<xyzVector> vecs, const xyzVector::geometryKind gk);
double triArea( const xyzVector& vecA, const xyzVector& vecB, const xyzVector& vecC, const xyzVector::geometryKind gk);
  
#endif /* defined(__LPM__xyzVector__) */
