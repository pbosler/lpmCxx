#ifndef __LPM_PARTICLES_HPP__
#define __LPM_PARTICLES_HPP__

#include <iostream>
#include <vector>
#include <string>
#include "OutputMessage.h"
#include "Logger.h"
#include "LpmXyzVector.hpp"
#include "LpmCoords.hpp"
#include "LpmEuclideanCoords.hpp"
#include "LpmSphereCoords.hpp"

using LpmXyzVector::XyzVector;

template <typename scalarType> class LpmParticles {
	public :
		enum geomType { PlanarCartesian, Euclidean3d, SphericalSurface };
		
		LpmParticles( const geomType geomKind, const size_t nMax, const int procRank = 0, const int nProcs = 1) {
			log = Logger::Instance( OutputMessage::debugPriority, procRank, nProcs );
			geometry = geomKind;
			switch ( geomKind ) {
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
				default : {
					OutputMessage errMsg( "ERROR: invalid geometry kind", OutputMessage::errorPriority,
						 "LpmParticles::LpmParticles" );
					log->logMessage(errMsg);
					break;
				}
			} 
		};
		
		~LpmParticles() { 
			delete physCoords;
			delete lagCoords;
		};
		
		geomType geometryKind() const { return geometry; }
		
		void insert( const XyzVector<scalarType>& physX, const XyzVector<scalarType>& lagX, 
					 const scalarType& areaOrVolume = 0.0 ) {
			physCoords->insert(physX);
			lagCoords->insert(lagX);
			switch ( geometry ) {
				case ( PlanarCartesian ) : {
					area.push_back(areaOrVolume);
					break;
				}
				case ( Euclidean3d ) : {
					volume.push_back(areaOrVolume);
					break;
				}
				case ( SphericalSurface ) : {
					area.push_back(areaOrVolume);
					break;
				};
			}
		};
		
		void insert( const scalarType xx, const scalarType yy, const scalarType zz = 0.0, 
					 const scalarType areaOrVolume = 0.0 ) {
			physCoords->insert( xx, yy, zz );
			lagCoords->insert( xx, yy, zz );
			switch ( geometry ) {
				case ( PlanarCartesian ) : {
					area.push_back(areaOrVolume);
					break;
				}
				case ( Euclidean3d ) : {
					volume.push_back(areaOrVolume);
					break;
				}
				case ( SphericalSurface ) : {
					area.push_back(areaOrVolume);
					break;
				};
			}
		};
		
		void replace( const size_t i, const XyzVector<scalarType>& physX, const XyzVector<scalarType> lagX,
				      const scalarType areaOrVolume = 0.0 ) {
			physCoords->replace(i, physX);
			lagCoords->repalce(i, lagX);
			switch ( geometry ) {
				case ( PlanarCartesian ) : {
					area[i] = areaOrVolume;
					break;
				}
				case ( Euclidean3d ) : {
					volume[i] = areaOrVolume;
					break;
				}
				case ( SphericalSurface ) : {
					area[i] = areaOrVolume;
					break;
				};
			}	
		};
			
		XyzVector<scalarType> physCoordVec( const size_t i ) const {
			return physCoords->coordVector( i );
		};
		
		XyzVector<scalarType> lagCoordVec( const size_t i ) const {
			return lagCoords->coordVector( i );
		};
			
		void setArea( const size_t i, const scalarType newArea ) { area[i] = newArea; }
		void setVolume( const size_t i, const scalarType newVol ) { volume[i] = newVol; }
		
		scalarType area( const size_t i ) const { return area[i]; }
		scalarType volume( const size_t i ) const { return volume[i]; }
		
		scalarType Latitude( const size_t i ) const { 
			return dynamic_cast<LpmSphereCoords<scalarType>*>(physCoords)->Latitude(i);
		}
		
		scalarType Longitude( const size_t i ) const {
			return dynamic_cast<LpmSphereCoords<scalarType>*>(physCoords)->Longitude(i);
		}
			
		size_t size() const { return physCoords->size(); }	
		int nDim() const { return physCoords->nDim(); }
		size_t nMax() const { return physCoords->nMax(); }
		
		scalarType distance( const size_t indexA, const size_t indexB ) const {	
			return physCoords->distance(indexA, indexB); 
		};
		
		XyzVector<scalarType> midpoint( const size_t indexA, const size_t indexB ) const {
			return physCoords->midpoint(indexA, indexB);}
		};
		
		scalarType triArea( const size_t indexA, const size_t indexB, const size_t indexC ) const {
			return physCoords->triArea(indexA, indexB, indexC);
		};
		
		XyzVector<scalarType> lagMidpoint( const size_t indexA, const size_t indexB ) const {
			return lagCoords->midpoint(indexA, indexB); }
		
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
		
		void PrintStats( const std::string& codeLoc ) const { 
			std::stringstream ss;
			ss << "Particles Stats : " << std::endl 
			   << "\t nMax = " << nMax() << std::endl
			   << "\t nDim = " << nDim() << std::endl
			   << "\t size = " << size() << std::endl
			   << "\t geomKind = " << geomString();
			OutputMessage statusMsg( ss.str(), OutputMessage::remarkPriority, codeLoc );
			log->logMessage(statusMsg);
		};
		
		void writePhysCoordsToMatlab( std::ostream& fs ) const {
			if ( geometry != SphericalSurface ) {
				XyzVector<scalarType> physX;
				fs << "x = [";
				for ( size_t i = 0; i < size() - 1; ++i ) {
					physX = physCoordVec(i);
					fs << physX.x << ", ";
				}
				physX = physCoordVec( size() - 1);
				fs << physX.x << "];" << std::endl << std::endl;
			
				fs << "y = [";
				for ( size_t i = 0; i < size() - 1; ++i ) {
					physX = physCoordVec(i);
					fs << physX.y << ", ";
				}
				physX = physCoordVec(size() - 1);
				fs << physX.y << "];" << std::endl << std::endl;
				if ( geometry != PlanarCartesian ) {
					fs << "z = ["; 
					for ( size_t i = 0; i < size() - 1; ++i ) {
						physX = physCoordVec(i);
						fs << physX.z << ", ";
					}
					physX = physCoordVec(size() - 1);
					fs << physX.z << "];" << std::endl << std::endl;	
				}
			}
			else {
				fs << "lat = [";
				for ( size_t i = 0; i < size() - 1; ++i ) {
					fs << dynamic_cast<LpmSphereCoords<scalarType>*>(physCoords)->Latitude(i) << ", ";
				}
				fs << dynamic_cast<LpmSphereCoords<scalarType>*>(physCoords)->Latitude(size() - 1) << "];" << std::endl << std::endl;
				
				fs << "lon = [";
				for ( size_t i = 0; i < size() - 1; ++i ) {
					fs << dynamic_cast<LpmSphereCoords<scalarType>*>(physCoords)->Longitude(i) << ", ";
				}
				fs << dynamic_cast<LpmSphereCoords<scalarType>*>(physCoords)->Longitude(size()-1) << "];" << std::endl << std::endl;
			}
		};
		
		void writeLagCoordsToMatlab( std::ostream& fs ) const {
			XyzVector<scalarType> lagX;
			fs << "x0 = [";
			for ( size_t i = 0; i < size() - 1; ++i ) {
				lagX = lagCoordVec(i);
				fs << lagX.x << ", ";
			}
			lagX = lagCoordVec( size() - 1);
			fs << lagX.x << "];" << std::endl << std::endl;
			
			fs << "y0 = [";
			for ( size_t i = 0; i < size() - 1; ++i ) {
				lagX = lagCoordVec(i);
				fs << lagX.y << ", ";
			}
			lagX = lagCoordVec(size() - 1);
			fs << lagX.y << "];" << std::endl << std::endl;
			if ( geometry != PlanarCartesian ) {
				fs << "z0 = ["; 
				for ( size_t i = 0; i < size() - 1; ++i ) {
					lagX = lagCoordVec(i);
					fs << lagX.z << ", ";
				}
				lagX = lagCoordVec(size() - 1);
				fs << lagX.z << "];" << std::endl << std::endl;	
			}
		};
		
		void writeAreaToMatlab( std::ostream& fs ) const {
			fs << "area = [";
			for ( size_t i = 0; i < size() - 1; ++i )
				fs << area[i] << ", ";
			fs << area[size()-1] << "];" << std::endl << std::endl;
		};

		void writeVolumeToMatlab( std::ostream& fs ) const {
			fs << "volume = [";
			for ( size_t i = 0; i < size() - 1; ++i )
				fs << volume[i] << ", ";
			fs << volume[size()-1] << "];" << std::endl << std::endl;
		};	
		
	protected :	
		geomType geometry;
		int _numProcs;
		int _procRank;
		
		std::string geomString() const {
			std::string result; 
			switch ( geometry ) {
				case ( PlanarCartesian ) : {
					result = "PlanarCartesian";
					break;
				}
				case ( Euclidean3d ) : {
					result = "Euclidean3d";
					break;
				}
				case ( SphericalSurface ) : {
					result = "SphericalSurface";
					break;
				};
			}
			return result;
		};
		
		Logger* log;
	
		LpmCoords<scalarType>* physCoords;
		LpmCoords<scalarType>* lagCoords;
		std::vector<scalarType> area;
		std::vector<scalarType> volume;
	
};

#endif
