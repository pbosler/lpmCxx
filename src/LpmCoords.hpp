#ifndef __LPM_COORDS_HPP__
#define __LPM_COORDS_HPP__

#include <iostream>
#include <cmath>
#include <vector>
#include <sstream>
#include "OutputMessage.h"
#include "Logger.h"
#include "GlobalConstants.h"
#include "LpmXyzVector.hpp"
/*
Base class for LPM Coordinates; will be subclassed for Euclidean or spherical coords.
*/
using LpmXyzVector::XyzVector;

template <typename scalarType> class LpmCoords {
	public :
		XyzVector<scalarType> coordVector( const size_t i ) const {
			const scalarType xx = x[i];
			const scalarType yy = y[i];
			const scalarType zz = (_nDim == 3) ? z[i] : 0.0;
			return XyzVector<scalarType>(xx,yy,zz);
		};
		
		size_t size() const { return x.size(); }
		size_t nMax() const { return _nMax; }
		int nDim() const { return _nDim; }
		
		void insert( const scalarType newx = 0.0, const scalarType newy = 0.0, const scalarType newz = 0.0) {
			x.push_back(newx);
			y.push_back(newy);
			if ( _nDim == 3 )
				z.push_back(newz);
		};
		
		void insert( const XyzVector<scalarType>& vec ) {
			x.push_back( vec.x );
			y.push_back( vec.y );
			if ( _nDim == 3 )
				z.push_back( vec.z );
		};
		
		void replace( const size_t index, const scalarType newx, const scalarType newy, const scalarType newz = 0.0) {
			x[index] = newx;
			y[index] = newy;
			if ( _nDim == 3 )
				z[index] = newz;
		};
		
		void replace( const size_t index, const XyzVector<scalarType>& vec ) {
			x[index] = vec.x;
			y[index] = vec.y;
			if ( _nDim == 3 )
				z[index] == vec.z;
		};
		
		scalarType magnitude( const size_t i ) const { 
			XyzVector<scalarType> vec = coordVector(i);
			return vec.magnitude();
		};
		
		void normalizeAll() {
			for ( size_t i = 0; i < x.size(); ++i ) {
				XyzVector<scalarType> vec = coordVector(i);
				vec.normalize();
				x[i] = vec.x;
				y[i] = vec.y;
				if ( _nDim == 3 )
					z[i] = vec.z;
			}
		};
		
		void scalarMultiplyAll( const scalarType multiplier ) {
			for ( size_t i = 0; i < x.size(); ++i ) {
				x[i] *= multiplier;
				y[i] *= multiplier;
				if (_nDim == 3 )
					z[i] *= multiplier;
			}
		};
		
		scalarType dotProduct( const size_t indexA, const size_t indexB ) {
			XyzVector<scalarType> vecA = coordVector(indexA);
			XyzVector<scalarType> vecB = coordVector(indexB);
			return vecA.dotProduct(vecB);
		};
		
		XyzVector<scalarType> crossProduct( const size_t indexA, const size_t indexB ) const {
			XyzVector<scalarType> vecA = coordVector(indexA);
			XyzVector<scalarType> vecB = coordVector(indexB);
			return vecA.crossProduct(vecB);
		};
		
		virtual scalarType distance( const size_t indexA, const size_t indexB ) const = 0;
		virtual XyzVector<scalarType> midpoint( const size_t indexA, const size_t indexB ) const = 0;
		virtual XyzVector<scalarType> centroid( std::vector<size_t> indices ) const = 0;
		virtual scalarType triArea( const size_t indexA, const size_t indexB, const size_t indexC ) const = 0;
		virtual ~LpmCoords() {};
		
	protected :
		LpmCoords( const int nDim, const size_t nMax, const int pRank, const int nProcs ) : _nMax(nMax), _nDim(nDim) {
			log = Logger::Instance( OutputMessage::warningPriority, pRank, nProcs );
			switch (nDim) {
				case (2) : {
					x.reserve(nMax);
					y.reserve(nMax);
					break;
					}
				case (3) : {
					x.reserve(nMax);
					y.reserve(nMax);
					z.reserve(nMax);
					break;
					}
				default : {
					std::stringstream ss;
					ss << "ERROR: nDim must be 2 or 3; input nDim = " << nDim ;
					OutputMessage errMsg( ss.str(), OutputMessage::errorPriority, "LpmCoords::LpmCoords" );
					log->logMessage(errMsg);
					break;
					}
			}
		};
	
		size_t _nMax;
		int _nDim;
		
		std::vector<scalarType> x;
		std::vector<scalarType> y;
		std::vector<scalarType> z;
		
		Logger* log;
};

template <typename scalarType> std::ostream& operator << ( std::ostream& os, const LpmCoords<scalarType>& coords ) {
	 os << "coords._nMax = " << coords.nMax() << "; current size = " << coords.size() << std::endl; 
	 os << "coords._nDim = " << coords.nDim() << std::endl;
	 for ( size_t i = 0; i < coords.size(); ++i ) {
	 	XyzVector<scalarType> vec = coords.coordVector(i);
	 	os << i << " : " << vec << std::endl;
	 }
	 return os; 
}

#endif 
