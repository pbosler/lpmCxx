#ifndef __LPM_EUCLIDEAN_COORDS__
#define __LPM_EUCLIDEAN_COORDS__

#include <iostream>
#include <cmath>
#include <vector>
#include "LpmXyzVector.hpp"
#include "LpmCoords.hpp"

using LpmXyzVector::XyzVector;

template <typename scalarType> class LpmEuclideanCoords : public LpmCoords<scalarType>
{
	public :
		LpmEuclideanCoords( const int nDim = 3, const size_t nMax = 0, const int pRank = 0, const int nProcs = 1) : 
			LpmCoords<scalarType>( nDim, nMax, pRank, nProcs ) {};			
		
		scalarType distance( const size_t indexA, const size_t indexB ) const {
			XyzVector<scalarType> vecA = LpmCoords<scalarType>::coordVector(indexA);
			XyzVector<scalarType> vecB = LpmCoords<scalarType>::coordVector(indexB);
			return LpmXyzVector::distance( vecA, vecB );
		};
		
		XyzVector<scalarType> midpoint( const size_t indexA, const size_t indexB ) const {
			XyzVector<scalarType> vecA = LpmCoords<scalarType>::coordVector(indexA);
			XyzVector<scalarType> vecB = LpmCoords<scalarType>::coordVector(indexB);
			return LpmXyzVector::midpoint( vecA, vecB );
		};
		
		XyzVector<scalarType> centroid( const std::vector<size_t> indices ) const {
			std::vector< XyzVector<scalarType> > vecs;
			for ( size_t i = 0; i < indices.size(); ++i ) {
				vecs.push_back( XyzVector<scalarType>( LpmCoords<scalarType>::coordVector(indices[i] ) ) );
			}
			return LpmXyzVector::centroid( vecs );
		};
		
		scalarType triArea( const size_t indexA, const size_t indexB, const size_t indexC ) const {
			XyzVector<scalarType> vecA = LpmCoords<scalarType>::coordVector(indexA);
			XyzVector<scalarType> vecB = LpmCoords<scalarType>::coordVector(indexB);
			XyzVector<scalarType> vecC = LpmCoords<scalarType>::coordVector(indexC);
			return LpmXyzVector::triArea( vecA, vecB, vecC);
		};
};



#endif 
