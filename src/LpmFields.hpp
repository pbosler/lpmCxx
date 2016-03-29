#ifndef __LPM_FIELDS_HPP__
#define __LPM_FIELDS_HPP__

#include <iostream>
#include <vector>
#include <string>
#include "OutputMessage.h"
#include "Logger.h"
#include "LpmXyzVector.hpp"
#include "LpmCoords.hpp"
#include "LpmParticles.hpp"

template <typename scalarType> class LpmField {
	public :
		LpmField( const int nDim, const int nMax, const std::string name = "", const std::string units = "", 
			const OutputMessage::priority logLevel = OutputMessage::debugPriority, const int nProcs = 1; const int pRank = 0) :	
			_nDim(nDim), _nMax(nMax), _name(name), _units(units) {
				log = Logger::Instance( logLevel, pRank, nProcs );
				switch (nDim) {
					case (1) : {
						scalar.reserve(nMax);
					}
					break;
					case (2) : {
						xComp.reserve(nMax);
						yComp.reserve(nMax);
					}
					break;
					case (3) : {
						xComp.reserve(nMax);
						yComp.reserve(nMax);
						zComp.reserve(nMax);
					}
					break;
					default : {
						OutputMessage errMsg("Constructor ERROR: invalid nDim", OutputMessage::errorPriority, "LpmParticles::LpmParticles");
						log->logMessage(errMsg);
					}
					break;
				}
			}
	
		size_t size() const { _nDim > 1 ? return xComp.size() : return scalar.size(); }
		int nDim() const { return _nDim; }
		int nMax() const { return _nMax; }
		std::string name() const { return _name; }
		std::string units() const {return _units; }
		
		void setName( const std::string newName ) { _name = newName; }
		void setUnits( const std::string newUnits ) { _units = newUnits; }
		
		
	protected :
		std::vector<scalarType> scalar;
		std::vector<scalarType> xComp;
		std::vector<scalarType> yComp;
		std::vector<scalarType> zComp;
		int _nDim;
		int _nMax;		
		std::string _name;
		std::string _units;
		Logger* log;
		OutputMessage::priority logLevel;
};


#endif
