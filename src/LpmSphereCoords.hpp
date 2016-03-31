#ifndef __LPM_SphereCoords__
#define __LPM_SphereCoords__

#include <iostream>
#include <cmath>
#include <vector>
#include "GlobalConstants.h"
#include "LpmXyzVector.hpp"
#include "LpmCoords.hpp"

using LpmXyzVector::XyzVector;

template <typename scalarType> class LpmSphereCoords : public LpmCoords<scalarType> {
	public :
		LpmSphereCoords( const int nMax, const scalarType radius = 1.0) : 
			LpmCoords<scalarType>( 3, nMax), _radius(radius) {};
		
		scalarType distance( const size_t indexA, const size_t indexB ) const {
			XyzVector<scalarType> vecA = LpmCoords<scalarType>::coordVector(indexA);
			XyzVector<scalarType> vecB = LpmCoords<scalarType>::coordVector(indexB);
			return LpmXyzVector::sphereDistance(vecA, vecB, _radius);
		};
		
		XyzVector<scalarType> midpoint( const size_t indexA, const size_t indexB ) const {
			XyzVector<scalarType> vecA = LpmCoords<scalarType>::coordVector(indexA);
			XyzVector<scalarType> vecB = LpmCoords<scalarType>::coordVector(indexB);
			return LpmXyzVector::sphereMidpoint(vecA, vecB, _radius);
		};
		
		XyzVector<scalarType> centroid( const std::vector<size_t> indices ) const {
			std::vector< XyzVector<scalarType> > vecs;
			for ( size_t i = 0; i < indices.size(); ++i )
				vecs.push_back( LpmCoords<scalarType>::coordVector(indices[i]) );
			return LpmXyzVector::sphereCentroid(vecs, _radius);
		};
		
		scalarType triArea( const size_t indexA, const size_t indexB, const size_t indexC ) const {
			XyzVector<scalarType> vecA = LpmCoords<scalarType>::coordVector(indexA);
			XyzVector<scalarType> vecB = LpmCoords<scalarType>::coordVector(indexB);
			XyzVector<scalarType> vecC = LpmCoords<scalarType>::coordVector(indexC);
			return LpmXyzVector::sphereTriArea( vecA, vecB, vecC, _radius);
		};
		
		scalarType Latitude( const size_t index ) const {
			XyzVector<scalarType> vec = LpmCoords<scalarType>::coordVector(index);
			return std::atan2( vec.z, std::sqrt( vec.x * vec.x + vec.y * vec.y ) );
		};
		
		scalarType Longitude( const size_t index ) const {
			XyzVector<scalarType> vec = LpmCoords<scalarType>::coordVector(index);
			return std::atan2( vec.y, vec.x );
		};
		
	protected : 
		scalarType _radius;
};

#endif

