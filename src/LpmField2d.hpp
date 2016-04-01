#ifndef __LPM_FIELD2D_HPP__
#define __LPM_FIELD2D_HPP__

#include <iostream>
#include <vector>
#include <string>
#include "OutputMessage.h"
#include "Logger.h"
#include "LpmXyzVector.hpp"
#include "LpmParticles.hpp"
#include "LpmScalarField.hpp"

using LpmXyzVector::XyzVector;

template <typename scalarType> class LpmField2d : public LpmScalarField<scalarType> {
	public :
		LpmField2d( const size_t nMax, const std::string name = " ", const std::string units = " ", 
			const int pRank = 0, const int nProcs = 1 ) : LpmScalarField<scalarType>(nMax, name, units, pRank, nProcs) {
			vals1.reserve(nMax);
		};
		
		LpmField2d( const LpmParticles<scalarType>& particles, const std::string name = " ", 
			const std::string units = " ", const int pRank = 0, const int nProcs = 1 ) :
			LpmScalarField<scalarType>( particles, name, units, pRank, nProcs ) {
				vals1.reserve(particles.size());
			};
		
		virtual ~LpmField2d() {};
	
		virtual void insert( const scalarType xval, const scalarType yval ) {
			vals.push_back(xval);
			vals1.push_back(yval);	
		};
		
		virtual void insert( const XyzVector<scalarType>& vecVal ) {
			vals.push_back(vecVal.x);
			vals1.push_back(vecVal.y);
		};
		
		virtual void replace( const size_t i, const scalarType xval, const scalarType yval ) {
			vals[i] = xval;
			vals1[i] = yval;
		};
		
		virtual void replace( const size_t i, const XyzVector<scalarType>& vecVal ) {
			vals[i] = vecVal.x;
			vals1[i] = vecVal.y;
		};
		
		virtual scalarType xComp( const size_t i ) const { return vals[i]; }
		virtual scalarType yComp( const size_t i ) const { return vals1[i]; }
		
		virtual XyzVector<scalarType> getVal( const size_t i ) const { 
			return XyzVector<scalarType>(vals[i], vals1[i]); 
		}
		
	protected :
		std::vector<scalarType> vals1;

};


#endif
