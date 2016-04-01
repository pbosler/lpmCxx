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
		LpmField2d( const size_t nMax, const std::string name = "noName", const std::string units = "null", 
			const int pRank = 0, const int nProcs = 1 ) : LpmScalarField<scalarType>(nMax, name, units, pRank, nProcs) {
			vals1.reserve(nMax);
		};
		
		LpmField2d( const LpmParticles<scalarType>& particles, const std::string name = "noName", 
			const std::string units = "null", const int pRank = 0, const int nProcs = 1 ) :
			LpmScalarField<scalarType>( particles, name, units, pRank, nProcs ) {
				vals1.reserve(particles.size());
			};
		
		virtual ~LpmField2d() {};
	
		virtual void insert( const scalarType xval, const scalarType yval ) {
			LpmScalarField<scalarType>::vals.push_back(xval);
			vals1.push_back(yval);	
		};
		
		virtual void insert( const XyzVector<scalarType>& vecVal ) {
			LpmScalarField<scalarType>::vals.push_back(vecVal.x);
			vals1.push_back(vecVal.y);
		};
		
		virtual void replace( const size_t i, const scalarType xval, const scalarType yval ) {
			LpmScalarField<scalarType>::vals[i] = xval;
			vals1[i] = yval;
		};
		
		virtual void replace( const size_t i, const XyzVector<scalarType>& vecVal ) {
			LpmScalarField<scalarType>::vals[i] = vecVal.x;
			vals1[i] = vecVal.y;
		};
		
		virtual scalarType xComp( const size_t i ) const { return LpmScalarField<scalarType>::vals[i]; }
		virtual scalarType yComp( const size_t i ) const { return vals1[i]; }
		
		virtual XyzVector<scalarType> getVectorVal( const size_t i ) const { 
			return XyzVector<scalarType>(LpmScalarField<scalarType>::vals[i], vals1[i]); 
		}
		
		virtual void PrintStats( const std::string& codeLoc ) const {
			std::stringstream ss;
			ss << "Field2d " << LpmScalarField<scalarType>::_name << " Stats : " << std::endl
			   << "\t units = " << LpmScalarField<scalarType>::_units << std::endl
			   << "\t nMax = " << LpmScalarField<scalarType>::_nMax << std::endl
			   << "\t size = " << vals1.size();
			OutputMessage statusMsg( ss.str(), OutputMessage::remarkPriority, codeLoc );
			LpmScalarField<scalarType>::log->logMessage(statusMsg);
		};
		
		virtual void writeFieldToMatlab( std::ostream& fs ) const {
			fs << LpmScalarField<scalarType>::_name << "_" << LpmScalarField<scalarType>::_units << " = [";
			for ( size_t i = 0; i < vals1.size() - 1; ++i) {
				fs << LpmScalarField<scalarType>::vals[i] << ", " << vals1[i] << "; ";
			}
			fs << LpmScalarField<scalarType>::vals[vals1.size() - 1] << ", " 
			   << vals1[vals1.size() -1] << "];" << std::endl << std::endl;
		};
		
	protected :
		std::vector<scalarType> vals1;

};


#endif
