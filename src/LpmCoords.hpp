#ifndef __LPM_Coords__
#define __LPM_Coords__

#include <iostream>
#include <cmath>
#include <vector>
#include "GlobalConstants.h"
#include "LpmXyzVector.hpp"

template <typename scalarType> class LpmCoords
{
	public :
		LpmCoords( const int nDim, const int nMax ) : _nDim(nDim), _nMax(nMax) { 
			switch (nDim) {
				case (2) : {
					x.reserve(nMax);
					y.reserve(nMax);
					}
					break;
				case (3) : {
					x.reserve(nMax);
					y.reserve(nMax);
					z.reserve(nMax);
					}
					break;
				default : {
					std::cerr << "LpmCoords::LpmCoords ERROR: nDim must be 2 or 3; nDim = " << nDim << " on input." << std::endl;
					break;
				}					
			}
		}			
		
		size_t size() const { return x.size(); }
		
		XyzVector<scalarType>& operator[](std::size_t i) { 
			XyzVector<scalarType> result;
			switch (_nDim ) {
				case (2) : 
					result = XyzVector<scalarType>(x[i], y[i]); 
					break;
				case (3) : 
					result = XyzVector<scalarType>(x[i], y[i], z[i]); 
					break;
			}
			return result;		
		}
		
		const XyzVector<scalarType>& operator[](std::size_t i) const { 
			XyzVector<scalarType> result;
			switch (_nDim ) {
				case (2) : 
					result = XyzVector<scalarType>(x[i], y[i]); 
					break;
				case (3) : 
					result = XyzVector<scalarType>(x[i], y[i], z[i]); 
					break;
			}
			return result;
		}
		
		void insert( const scalarType nx = 0.0, const scalarType ny = 0.0, const scalarType nz = 0.0 ) {
			x.push_back(nx);
			y.push_back(ny);
			if ( _nDim != 2 )
				z.push_back(nz);
		};
		void insert( const XyzVector<scalarType>& vec ) {
			x.push_back( vec.x );
			y.push_back( vec.y );
			if ( _nDim != 2 )
				z.push_back( vec.z );
		};
		
		void replace( const scalarType nx, const scalarType ny, const scalarType nz, const size_t index ) {
			x[index] = nx;
			y[index] = ny;
			if ( _nDim != 2 )
				z[index] = nz;
		};
		void replace( const XyzVector<scalarType>& vec, const size_t index ) {
			x[index] = vec.x;
			y[index] = vec.y;
			if ( _nDim != 2 )
				z[index] = vec.z;
		};
		
		scalarType dotProduct( const size_t indexA, const size_t indexB ) const {
			XyzVector<scalarType> vecA = *this[indexA];
			XyzVector<scalarType> vecB = *this[indexB];
			return vecA.dotProduct(vecB);
		};
		
		XyzVector<scalarType> crossProduct( const size_t indexA, const size_t indexB ) const {
			XyzVector<scalarType> vecA = *this[indexA];
			XyzVector<scalarType> vecB = *this[indexB];
			return vecA.crossProduct(vecB);
		};
		
		void normalizeAll() {
			for ( size_t i = 0; i < x.size(); ++i ) {
				XyzVector<scalarType> vec = *this[i];
				vec.normalize();
				x[i] = vec.x;
				y[i] = vec.y;
				if ( _nDim != 2 )
					z[i] = vec.z;
			}
		};
		
		void scalarMultiplyAll( const scalarType multiplier ) {
			switch (_nDim ) {
				case (2) : {
					for ( size_t i = 0; i < x.size(); ++i ) {
						x[i] *= multiplier;
						y[i] *= multiplier;
						z[i] *= multiplier;
					}
				}
				break;
				case (3) : {
					for ( size_t i = 0; i < x.size(); ++i ) {
						x[i] *= multiplier;
						y[i] *= multiplier;
						z[i] *= multiplier;
					}
				break;
				}
			}
		};
		
		virtual scalarType distance( const size_t indexA, const size_t indexB ) const {
			XyzVector<scalarType> vecA = *this[indexA];
			XyzVector<scalarType> vecB = *this[indexB];
			return distance( vecA, vecB );
		};
		
		virtual XyzVector<scalarType> midpoint( const size_t indexA, const size_t indexB ) const {
			XyzVector<scalarType> vecA = *this[indexA];
			XyzVector<scalarType> vecB = *this[indexB];
			return midpoint( vecA, vecB );
		};
		
		virtual XyzVector<scalarType> centroid( const std::vector<size_t> indices ) const {
			std::vector< XyzVector<scalarType> > vecs;
			for ( size_t i = 0; i < indices.size(); ++i ) {
				vecs.push_back( XyzVector<scalarType>( *this[indices[i]] ) );
			}
			return centroid( vecs );
		};
		
		virtual scalarType triArea( const size_t indexA, const size_t indexB, const size_t indexC ) const {
			XyzVector<scalarType> vecA = *this[indexA];
			XyzVector<scalarType> vecB = *this[indexB];
			XyzVector<scalarType> vecC = *this[indexC];
			return triArea( vecA, vecB, vecC);
		};
		
	protected :
		int _nMax;
		int _nDim;
		
		std::vector<scalarType> x;
		std::vector<scalarType> y;
		std::vector<scalarType> z;
};

#endif 
