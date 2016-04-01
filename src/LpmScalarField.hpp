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
		LpmScalarField( const size_t nMax, const std::string name = "noName", const std::string units = "null", 
						const int pRank = 0, const int nProcs = 1) : _nMax(nMax), _name(name), _units(units) {
			log = Logger::Instance( OutputMessage::debugPriority, pRank, nProcs );
			vals.reserve(nMax);				
		};
		
		LpmScalarField( const LpmParticles<scalarType>& particles, const std::string name = "noName", 
					    const std::string units = "null", const int pRank = 0, const int nProcs = 1 ) :
					    _nMax(particles.size()), _name(name), _units(units) {
			log = Logger::Instance( OutputMessage::debugPriority, pRank, nProcs );
			vals.reserve(particles.size());
		};
		
		virtual ~LpmScalarField() {};
		
		std::string name() const { return _name; }
		std::string units() const { return _units; }
		size_t nMax() const { return _nMax; }
		size_t size() const { return vals.size(); }
		
		virtual void insert( const scalarType newVal ) { vals.push_back(newVal); }
		
		virtual void replace( const size_t i, const scalarType newVal ) { vals[i] = newVal; }
		
		scalarType getScalarVal( const size_t i ) const { return vals[i]; }
		
		virtual void writeFieldToMatlab( std::ostream& fs ) const {
			fs << _name << "_" << _units << " = [";
			for ( size_t i = 0; i < vals.size() - 1; ++i ) {
				fs << vals[i] << ", ";
			}
			fs << vals[ vals.size() - 1] << "];" << std::endl << std::endl;
		};
	
		virtual void PrintStats(const std::string& codeLoc ) const { 
			std::stringstream ss;
			ss << "ScalarField " << _name << " Stats : " << std::endl
			   << "\t units = " << _units << std::endl
			   << "\t nMax = " << _nMax << std::endl
			   << "\t size = " << vals.size();
			OutputMessage statusMsg( ss.str(), OutputMessage::remarkPriority, codeLoc);
			log->logMessage(statusMsg);
			
		};
	protected :
		size_t _nMax;
		
		std::string _name;
		std::string _units;
		
		Logger* log;
		
		std::vector<scalarType> vals;
	
};

#endif
