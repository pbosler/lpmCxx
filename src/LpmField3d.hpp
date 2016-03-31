#ifndef __LPM_FIELD3D_HPP__
#define __LPM_FIELD3D_HPP__

#include <iostream>
#include <vector>
#include <string>
#include "OutputMessage.h"
#include "Logger.h"
#include "LpmXyzVector.hpp"
#include "LpmParticles.hpp"
#include "LpmField2d.hpp"

using LpmXyzVector::XyzVector;

template <typename scalarType> class LpmField3d : public LpmField2d<scalarType> {
	public :
		LpmField3d( const int nMax, const std::string name = " ", const std::string units = " ", 
			const int pRank = 0, const int nProcs = 1 ) : LpmField2d<scalarType>(nMax, name, units, pRank, nProcs) {
			vals2.reserve(nMax);
		};
		
		LpmField3d( const LpmParticles<scalarType>& particles, const std::string name = " ", 
			const std::string units = " ", const int pRank = 0, const int nProcs = 1 ) :
			LpmField2d<scalarType>( particles, name, units, pRank, nProcs ) {
				vals2.reserve(nMax);
			};
		
		void insert( const scalarType xval, const scalarType yval, const scalarType zval ) {
			vals.push_back(xval);
			vals1.push_back(yval);
			vals2.push_back(zval);
		};
		
		void insert( const XyzVector<scalarType>& vecVal ) {
			vals.push_back(vecVal.x);
			vals1.push_back(vecVal.y);
			vals2.push_back(vecVal.z);
		};
		
		void replace( const size_t i, const scalarType xval, const scalarType yval, const scalarType zval ) {
			vals[i] = xval;
			vals1[i] = yval;
			vals2[i] = zval;
		};
		
		void replace( const size_t i, const XyzVector<scalarType>& vecVal ) {
			vals[i] = vecVal.x;
			vals1[i] = vecVal.y;
			vals2[i] = vecVal.z;
		};
		
		scalarType xComp( const size_t i ) const { return vals[i]; }
		scalarType yComp( const size_t i ) const { return vals1[i]; }
		scalarType zComp( const size_t i ) const { return vals2[i]; }
		
		XyzVector<scalarType> getVal( const size_t i ) const { 
			return XyzVector<scalarType>(vals[i], vals1[i], z[i]); 
		}
		
	protected :
		std::vector<scalarType> vals2;
};

#endif
