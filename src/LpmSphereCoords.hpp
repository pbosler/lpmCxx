#ifndef __LPM_SphereCoords__
#define __LPM_SphereCoords__

#include <iostream>
#include <cmath>
#include <vector>
#include "GlobalConstants.h"
#include "LpmXyzVector.hpp"
#include "LpmCoords.hpp"

template <typename scalarType> class LpmSphereCoords : public LpmCoords<scalarType> {
	public :
		LpmSphereCoords( const int nDim, const int nMax, const scalarType radius = 1.0) : 
			LpmCoords<scalarType>(nDim, nMax), _radius(radius) {};
		
		scalarType distance( const size_t indexA, const size_t indexB ) const {
			XyzVector<scalarType> vecA = *this[indexA];
			XyzVector<scalarType> vecB = *this[indexB];
			return sphereDistance(vecA, vecB, _radius);
		};
		
		XyzVector<scalarType> midpoint( const size_t indexA, const size_t indexB ) const {
			XyzVector<scalarType> vecA = *this[indexA];
			XyzVector<scalarType> vecB = *this[indexB];
			return sphereMidpoint(vecA, vecB, _radius);
		};
		
		XyzVector<scalarType> centroid( const std::vector<size_t> indices ) const {
			std::vector< XyzVector<scalarType> > vecs;
			for ( size_t i = 0; i < indices.size(); ++i )
				vecs.push_back( XyzVector<scalarType>( *this[indices[i]] ) );
			return sphereCentroid(vecs, _radius);
		};
		
		scalarType triArea( const size_t indexA, const size_t indexB, const size_t indexC ) const {
			XyzVector<scalarType> vecA = *this[indexA];
			XyzVector<scalarType> vecB = *this[indexB];
			XyzVector<scalarType> vecC = *this[indexC];
			return sphereTriArea( vecA, vecB, vecC, _radius);
		};
		
		scalarType Latitude( const size_t index ) const {
			XyzVector<scalarType> vec = *this[index];
			return std::atan2( vec.z, std::sqrt( vec.x * vec.x + vec.y * vec.y ) );
		};
		
		scalarType Longitude( const size_t index ) const {
			XyzVector<scalarType> vec = *this[index];
			return std::atan2( vec.y, vec.x );
		};
		
	protected : 
		scalarType _radius;
};

#endif

