#ifndef __LPM_SCALARFIELD_HPP__
#define __LPM_SCALARFIELD_HPP__

#include <iostream>
#include <vector>
#include <string>
#include "OutputMessage.h"
#include "Logger.h"
#include "LpmXyzVector.hpp"
#include "LpmParticles.hpp"

using LpmXyzVector::XyzVector;

template <typename scalarType> class LpmScalarField {
	public :
		LpmScalarField( const size_t nMax, const std::string name = " ", const std::string units = " ", 
						const int pRank = 0, const int nProcs = 1) : _nMax(nMax), _name(name), _units(units) {
			log = Logger::Instance( OutputMessage::debugPriority, pRank, nProcs );
			vals.reserve(nMax);				
		};
		
		LpmScalarField( const LpmParticles<scalarType>& particles, const std::string name = " ", 
					    const std::string units = " ", const int pRank = 0, const int nProcs = 1 ) :
					    _nMax(particles.size()), _name(name), _units(units) {
			log = Logger::Instance( OutputMessage::debugPriority, pRank, nProcs );
			vals.reserve(particles.size());
		};
		
		virtual ~LpmScalarField() {};
		
		std::string name() const { return _name; }
		std::string units() const { return _units; }
		
		virtual void insert( const scalarType newVal ) { vals.push_back(newVal); }
		
		virtual void replace( const size_t i, const scalarType newVal ) { vals[i] = newVal; }
		
		virtual scalarType getVal( const size_t i ) const { return vals[i]; }
	
	protected :
		size_t _nMax;
		
		std::string _name;
		std::string _units;
		
		Logger* log;
		
		std::vector<scalarType> vals;
	
};

#endif
