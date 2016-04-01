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
		LpmField3d( const size_t nMax, const std::string name = "noName", const std::string units = "null", 
			const int pRank = 0, const int nProcs = 1 ) : LpmField2d<scalarType>(nMax, name, units, pRank, nProcs) {
			vals2.reserve(nMax);
		};
		
		LpmField3d( const LpmParticles<scalarType>& particles, const std::string name = "noName", 
			const std::string units = "null", const int pRank = 0, const int nProcs = 1 ) :
			LpmField2d<scalarType>( particles, name, units, pRank, nProcs ) {
				vals2.reserve(particles.size());
			};
		
		void insert( const scalarType xval, const scalarType yval, const scalarType zval ) {
			LpmScalarField<scalarType>::vals.push_back(xval);
			LpmField2d<scalarType>::vals1.push_back(yval);
			vals2.push_back(zval);
		};
		
		void insert( const XyzVector<scalarType>& vecVal ) {
			LpmScalarField<scalarType>::vals.push_back(vecVal.x);
			LpmField2d<scalarType>::vals1.push_back(vecVal.y);
			vals2.push_back(vecVal.z);
		};
		
		void replace( const size_t i, const scalarType xval, const scalarType yval, const scalarType zval ) {
			LpmScalarField<scalarType>::vals[i] = xval;
			LpmField2d<scalarType>::vals1[i] = yval;
			vals2[i] = zval;
		};
		
		void replace( const size_t i, const XyzVector<scalarType>& vecVal ) {
			LpmScalarField<scalarType>::vals[i] = vecVal.x;
			LpmField2d<scalarType>::vals1[i] = vecVal.y;
			vals2[i] = vecVal.z;
		};
		
		scalarType xComp( const size_t i ) const { return LpmScalarField<scalarType>::vals[i]; }
		scalarType yComp( const size_t i ) const { return LpmField2d<scalarType>::vals1[i]; }
		scalarType zComp( const size_t i ) const { return vals2[i]; }
		
		XyzVector<scalarType> getVal( const size_t i ) const { 
			return XyzVector<scalarType>(LpmScalarField<scalarType>::vals[i], LpmField2d<scalarType>::vals1[i], 
				vals2[i]); 
		}
		
		void PrintStats( const std::string& codeLoc ) const {
			std::stringstream ss;
			ss << "Field3d " << LpmScalarField<scalarType>::_name << " Stats : " << std::endl
			   << "\t units = " << LpmScalarField<scalarType>::_units << std::endl
			   << "\t nMax = " << LpmScalarField<scalarType>::_nMax << std::endl
			   << "\t size = " << vals2.size();
			OutputMessage statusMsg( ss.str(), OutputMessage::remarkPriority, codeLoc );
			LpmScalarField<scalarType>::log->logMessage(statusMsg);
		};
		
		void writeFieldToMatlab( std::ostream& fs ) const {
			fs << LpmScalarField<scalarType>::_name << "_" << LpmScalarField<scalarType>::_units << " = [";
			for ( size_t i = 0; i < vals2.size() - 1; ++i) {
				fs << LpmScalarField<scalarType>::vals[i] << ", " << LpmField2d<scalarType>::vals1[i] << ", "
				   << vals2[i] << "; ";
			}
			fs << LpmScalarField<scalarType>::vals[vals2.size() - 1] << ", " 
			   << LpmField2d<scalarType>::vals1[vals2.size() -1] << ", " 
			   << vals2[vals2.size() - 1] << "];" << std::endl << std::endl;
		};
	protected :
		std::vector<scalarType> vals2;
};

#endif
