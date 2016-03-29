#ifndef __LPM_PARTICLES_HPP__
#define __LPM_PARTICLES_HPP__

#include <iostream>
#include <vector>
#include <string>
#include "OutputMessage.h"
#include "Logger.h"
#include "LpmXyzVector.hpp"
#include "LpmCoords.hpp"

template <typename scalarType> class LpmParticles {
	public :
		LpmParticles( const int nDim = 3, const int nMax = 0, const int nProcs, const int pRank, 
					  const OutputMessage::priority logLevel = OutputMessage::debugPriority )
			 : _nDim(nDim), _nMax(nMax), physCoords( nDim, nMax ), lagCoords( nDim, nMax ), _numProcs(nProcs), _procRank(pRank)
			{ 
				log = Logger::Instance( logLevel, pRank, nProcs );
				switch ( nDim ) {
					case (2) : {
						area.reserve(nMax);
					}
					break;
					case (3) : {
						volume.reserve(nMax);
					}
					break;
					default : {
						OutputMessage errMsg("Constructor ERROR: invalid nDim", OutputMessage::errorPriority, "LpmParticles::LpmParticles");
						log->logMessage(errMsg);
						break;
					}
				}
			};
			
		size_t size() const { return physCoords.size(); }	
		
		scalarType totalVolume() const {
			scalarType result = 0.0;
			for (size_t i = 0; i < volume.size(); ++i )
				result += volume[i];
			return result;
		}
		
		scalarType totalArea() const {
			scalarType result = 0.0;
			for ( size_t i = 0; i < area.size(); ++i )
				result += area[i];
			return result;
		}
	
		void insert( const XyzVector<scalarType> physX, const XyzVector<scalarType> lagX, const scalarType areaOrVolume = 0.0 ) {
			physCoords.insert(physX);
			lagCoords.insert(physX);
			switch ( _nDim ) {
				case (2) : {
					area.push_back( areaOrVolume );
				}
				break;
				case (3) : {
					volume.push_back(areaOrVolume);
				}
				break;
			}
		};
		
		void insert( const scalarType nx = 0.0, const scalarType ny = 0.0, const scalarType nz = 0.0, const scalarType areaOrVolume = 0.0 ) {
			physCoords.insert( nx, ny, nz );
			lagCoords.insert(nx, ny, nz );
			switch ( _nDim ) {
				case (2) : {
					area.push_back(areaOrVolume);
				}
				break;
				case (3) : {
					volume.push_back(areaOrVolume);
				}
				break;
			}
		};
		
		void replace( const XyzVector<scalarType> physX, const XyzVector<scalarType> lagX, const size_t index, const scalarType areaOrVolume = 0.0 ) {
			physCoords.replace( physX, index );
			lagCoords.replace( lagX, index );
			switch ( _nDim ) {
				case (2) : {
					area[index] = areaOrVolume;
				}
				break;
				case (3) : {
					volume[index] = areaOrVolume;
				}
				break;
			}
		};
		
		scalarType dotProduct( const size_t indexA, const size_t indexB ) const {
			return physCoords.dotProduct( indexA, index B);
		};
		
		XyzVector<scalarType> crossProduct( const size_t indexA, const size_t indexB ) const {
			return physCoords.crossProduct( indexA, index B );
		};
		
		scalarType distance( const size_t indexA, const size_t indexB ) const {
			return physCoords.distance( indexA, indexB );
		};
		
		XyzVector<scalarType> midpoint( const size_t indexA, const size_t indexB ) const {
			return physCoords.midpoint(indexA, indexB);
		};
		
		XyzVector<scalarType> centroid( std::vector<size_t> indices ) const {
			return physCoords.centroid( indices
		};
		
	protected :
		LpmCoords<scalarType> physCoords;
		LpmCoords<scalarType> lagCoords;
		std::vector<scalarType> area;
		std::vector<scalarType> volume;
		
		int _nMax;
		int _nDim;
		int _numProcs;
		int _procRank;
		
		Logger* log;
};

#endif
