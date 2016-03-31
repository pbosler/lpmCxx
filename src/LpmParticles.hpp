#ifndef __LPM_PARTICLES_HPP__
#define __LPM_PARTICLES_HPP__

#include <iostream>
#include <vector>
#include <string>
#include "OutputMessage.h"
#include "Logger.h"
#include "LpmXyzVector.hpp"
#include "LpmCoords.hpp"

using LpmXyzVector::XyzVector;

template <typename scalarType> class LpmParticles {
	public :
		enum geomType { PlanarCartesian, Euclidean3d, SphericalSurface };
		
		LpmParticles( const geomType, const size_t nMax, const int procRank, const int nProcs ) {
			log = Logger::Instance( OutputMessage::debugPriority, procRank, nProcs );
			geometry = geomType;
			switch ( geomType ) {
				case ( PlanarCartesian ) : {
					physCoords = new LpmEuclideanCoords<scalarType>( 2, nMax, procRank, nProcs );
					lagCoords = new LpmEuclideanCoords<scalarType>(2, nMax, procRank, nProcs );
					area.reserve(nMax);
					break;
				}
				case ( Euclidean3d ) : {
					physCoords = new LpmEuclideanCoords<scalarType>( 3, nMax, procRank, nProcs );
					lagCoords = new LpmEuclideanCoords<scalarType>(3, nMax, procRank, nProcs );
					volume.reserve(nMax);
					break;
				}
				case ( SphericalSurface ) : {
					physCoords = new LpmSphereCoords<scalarType>( 3, nMax, procRank, nProcs );
					lagCoords = new LpmSphereCoords<scalarType>( 3, nMax, procRank, nProcs );
					area.reserve(nMax);
					break;
				}
			} 
		};
			
		size_t size() const { return physCoords.size(); }	
		int nDim() const { return physCoords.nDim(); }
		size_t nMax() const { return physCoords.nMax(); }
		
		scalarType totalVolume() const {
			scalarType result = 0.0;
			for (size_t i = 0; i < volume.size(); ++i )
				result += volume[i];
			return result;
		}

		scalarType totalVolume( std::vector<bool> isActive ) const {
			scalarType result = 0.0;
			for (size_t i = 0; i < volume.size(); ++i ) {
				if ( isActive[i] )
					result += volume[i];
			}
			return result;
		}
		
		scalarType totalArea() const {
			scalarType result = 0.0;
			for ( size_t i = 0; i < area.size(); ++i )
				result += area[i];
			return result;
		}
		
		scalarType totalArea( std::vector<bool> isActive ) const {
			scalarType result = 0.0;
			for ( size_t i = 0; i < area.size(); ++i ) {
				if ( isActive[i] ) 
					result += area[i];
			}
			return result;
		}
	
	protected :	
		geomType geometry;
		int _numProcs;
		int _procRank;
		
		Logger* log;
	
		LpmCoords<scalarType>* physCoords;
		LpmCoords<scalarType>* lagCoords;
		std::vector<scalarType> area;
		std::vector<scalarType> volume;
	
};

#endif
