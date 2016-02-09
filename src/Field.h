//
//  Field.h
//  LPM
//
//  Created by Peter Bosler on 11/1/14.
//  Copyright (c) 2014 Peter Bosler. All rights reserved.
//

/** @file Field.h 
	@brief Field class header file
	@author Peter Bosler, Sandia National Laboratories, Multiphysics Applications
*/

#ifndef __LPM__Fields__
#define __LPM__Fields__

#include <iostream>
#include <vector>
#include <string>
#include "OutputMessage.h"
#include "Logger.h"
#include "xyzVector.h"

class Particles;
class Faces;

/** @class Field
	@brief The Field class is a container for holding scalar and vector field data associated with a set of Particles.

	Each Field may be either a scalar field or a vector field.  
	Entries in a Field are stored in a 1-to-1 correspondence with an instance of the Particles class, so that a 
	scalar field entry associated with particle i is given by @ref _scalar[i].  
*/
class Field
{
public:
	std::vector<double> scalar; ///< scalar field data
    std::vector<double> xComp; ///< vector field data
    std::vector<double> yComp; ///< vector field data
    std::vector<double> zComp; ///< vector field data

    Field( const std::string name = " ", const std::string units = " ",
           const int nDim = 1, const int nMax = 20, const int procRank = 0, const int numProcs = 1);

	int nDim() const { return _nDim; }
	int N() const { return _N; }
	int nMax() const { return _nMax; }
	
	std::string name() const { return _name; }
	std::string units() const { return _units; }
	
	void setName( const std::string newName ){ _name = newName; }
	void setUnits( const std::string newUnits ){ _units = newUnits; }
	void setN( const int newN){ _N = newN;}
	
	double scalarFieldValue( const int index ) const { return scalar[index]; }
	xyzVector vectorFieldValue( const int index ) const ;
    
    void insertScalarToField( const int index, const double val);
    void insertVectorToField( const int index, const double vx, const double vy, const double vz = 0.0);
    void insertVectorToField( const int index, const xyzVector& vec);
     
    void initializeToScalarConstant( const double val = 0.0 );
    void initializeToVectorConstant( const xyzVector& vec );
     
    void sumIntoScalarField( const double val);
    void sumIntoVectorField( const xyzVector& vec);
           
    void abs();
    
    void outputForMatlab( const std::string filename, const Particles& particles ) const;
    
    void shrinkMemory();
    
    Field operator+ (Field& other );
    Field operator- (Field& other );
    
    double maximumMagnitude() const;
    
    void writePointDataToVTK( std::ostream& fs ) const;
    void writeCellDataToVTK( std::ostream& fs, const Faces& faces ) const;
    
protected:
    std::string _name;
    std::string _units;
    
    int _N;
    int _nMax;
    int _nDim;
    OutputMessage::priority logLevel;
    Logger* log;
};

//std::ostream& operator<<( std::ostream& os, const Field& afield);


#endif /* defined(__LPM__Fields__) */
