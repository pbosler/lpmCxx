//
//  Fields.cpp
//  LPM
//
//  Created by Peter Bosler on 11/1/14.
//  Copyright (c) 2014 Peter Bosler. All rights reserved.
//

/** @file Field.cpp 
	@brief Field class implementation
	@author Peter Bosler, Sandia National Laboratories, Multiphysics Applications
*/
#include "Field.h"
#include "xyzVector.h"
#include "Particles.h"
#include "Faces.h"
#include <vector>
#include <assert.h>
#include <cmath>
#include <fstream>

Field::Field( const std::string name, const std::string units, const int nDim, const int nMax,
			  const int procRank, const int numProcs)
{
	assert( nDim >=1 && nDim <= 3 );
	assert( procRank >= 0 );
	assert( numProcs >= 1 );
	
	_name = name;
	_units = units;
	_N = 0;
	_nMax = nMax;
	_nDim = nDim;
	
	logLevel = OutputMessage::debugPriority;
	log = Logger::Instance( logLevel, procRank, numProcs);
	
	if ( nDim == 1 )
		scalar = std::vector<double>(nMax, 0.0);
	else 
	{
		xComp = std::vector<double>(nMax, 0.0);
		yComp = std::vector<double>(nMax, 0.0);
		if ( nDim == 3)
			zComp = std::vector<double>(nMax, 0.0);
	}
};

void Field::abs()
{
	if ( _nDim == 1)
	{
		for ( int i = 0; i < _N; ++i)
		{
			scalar[i] = std::abs( scalar[i] );
		}
	}
	else
	{
		for ( int i = 0; i < _N; ++i)
		{
			xComp[i] = std::abs( xComp[i] );
			yComp[i] = std::abs( yComp[i] );
		}
		if ( _nDim == 3 )
		{
			for ( int i = 0; i < _N; ++i)
			{
				zComp[i] = std::abs( zComp[i] );
			}
		}
	}
}

xyzVector Field::vectorFieldValue( const int index ) const 
{
	assert( _nDim > 1);
	assert( index < _nMax );
	switch ( _nDim )
	{
		case 2 :
		{ 
			return xyzVector(xComp[index], yComp[index]);
			break;
		}
		case 3 :
		{
			return xyzVector(xComp[index], yComp[index], zComp[index]);
			break;
		}
		default:
		{
			return xyzVector();
		}
	}
};

void Field::insertScalarToField( const int index, const double val)
{
	assert( index < _nMax);
	scalar[index] = val;
	_N++;
};

void Field::insertVectorToField( const int index, const double vx, const double vy, const double vz)
{
	assert( index < _nMax);
	xComp[index] = vx;
	yComp[index] = vy;
	if ( _nDim == 3) 
		zComp[index] = vz;
	_N++;
};

void Field::insertVectorToField( const int index, const xyzVector& vec)
{
	assert(index < _nMax);
	xComp[index] = vec.x;
	yComp[index] = vec.y;
	if ( _nDim == 3) 
		zComp[index] = vec.z;
	_N++;
};

Field Field::operator+ ( Field& other )
{
	assert( _N == other._N );
	assert( _nDim == other._nDim );
	Field result( "sum", "units", _nDim, _N) ;
	switch ( _nDim )
	{
		case 1:
		{
			for ( int i = 0; i < _N; ++i )
			{
				result.scalar[i] = scalar[i] + other.scalar[i];
			}
		}
		case 2:
		{
			for ( int i = 0; i < _N; ++i )
			{
				result.xComp[i] = xComp[i] + other.xComp[i];
				result.yComp[i] = yComp[i] + other.yComp[i];	
			}
		}
		case 3:
		{
			for ( int i = 0; i < _N; ++i )
			{
				result.xComp[i] = xComp[i] + other.xComp[i];
				result.yComp[i] = yComp[i] + other.yComp[i];	
				result.zComp[i] = zComp[i] + other.zComp[i];
			}
		}
	}
	return result;
};

Field Field::operator- (Field& other )
{
	assert( _N == other._N );
	assert( _nDim == other._nDim );
	Field result( "difference", "units", _nDim, _N) ;
	switch ( _nDim )
	{
		case 1:
		{
			for ( int i = 0; i < _N; ++i )
			{
				result.scalar[i] = scalar[i] - other.scalar[i];
			}
		}
		case 2:
		{
			for ( int i = 0; i < _N; ++i )
			{
				result.xComp[i] = xComp[i] - other.xComp[i];
				result.yComp[i] = yComp[i] - other.yComp[i];	
			}
		}
		case 3:
		{
			for ( int i = 0; i < _N; ++i )
			{
				result.xComp[i] = xComp[i] - other.xComp[i];
				result.yComp[i] = yComp[i] - other.yComp[i];	
				result.zComp[i] = zComp[i] - other.zComp[i];
			}
		}
	}
	return result;
};

double Field::maximumMagnitude() const
{
	double result = 0.0;
	switch ( _nDim )
	{
		case 1:
		{
			for (int i = 0; i < _N; ++i )
			{
				if ( std::abs( scalar[i] ) > result )
					result = std::abs( scalar[i] );
			}
		}
		case 2:
		{
			for ( int i = 0; i < _N; ++i )
			{
				if ( std::sqrt( xComp[i]*xComp[i] + yComp[i]*yComp[i] ) > result )
					result = std::sqrt( xComp[i]*xComp[i] + yComp[i]*yComp[i] );
			}
		}
		case 3:
		{
			for ( int i = 0; i < _N; ++i )
			{
				if ( std::sqrt( xComp[i]*xComp[i] + yComp[i]*yComp[i] + zComp[i]*zComp[i] ) > result )
					result = std::sqrt( xComp[i]*xComp[i] + yComp[i]*yComp[i] + zComp[i]*zComp[i] );
			}
		}
	}
	return result;
};

void Field::initializeToScalarConstant( const double val )
{
	for ( int i = 0; i < _nMax; ++i )
		insertScalarToField(i, val);
};

void Field::initializeToVectorConstant( const xyzVector& vec )
{
	for ( int i = 0; i < _nMax; ++i )
		insertVectorToField(i, vec);
}

void Field::sumIntoScalarField( const double val )
{
	for ( int i = 0; i < _N; ++i )
		scalar[i] += val;
};

void Field::sumIntoVectorField( const xyzVector& vec )
{
	if ( _nDim > 1 )
	{
		for ( int i = 0; i < _N; ++i )
		{
			xComp[i] += vec.x;
			yComp[i] += vec.y;
		}
		if ( _nDim == 3 )
		{
			for ( int i = 0; i < _N; ++i )
				zComp[i] += vec.z;
		}
	}
	else
	{
		OutputMessage statusMsg("cannot add a vector to a scalar.", OutputMessage::errorPriority, "Field::sumIntoVectorField");
		log->logMessage(statusMsg);
		return;
	}
};

void Field::outputForMatlab( const std::string filename, const Particles& particles ) const 
{
	std::ofstream file( filename.c_str() );
	if ( !file )
	{
		OutputMessage statusMsg("cannot open file.", OutputMessage::errorPriority, "Field::outputForMatlab");
		log->logMessage(statusMsg);
		return;
	}
	
	file << "Field_name = '" << _name << "'; " << std::endl;
	file << "units = '" << _units << "';" <<std::endl;
	
	const int n = particles.N();
	
	file << " x = [ ";
	for ( int i = 0; i < n-1; ++i)
	{
		file << particles.x[i] << " , ";
	}
	file << particles.x[n-1] << " ]; " << std::endl << std::endl;
	
	file << " y = [ ";
	for ( int i = 0; i < n-1; ++i)
	{
		file << particles.y[i] << " , ";
	}
	file << particles.y[n-1] << " ]; " << std::endl << std::endl;
	
	file << " x0 = [ ";
	for ( int i = 0; i < n-1; ++i)
	{
		file << particles.x0[i] << " , ";
	}
	file << particles.x0[n-1] << " ]; " << std::endl << std::endl;
	
	file << " y0 = [ ";
	for ( int i = 0; i < n-1; ++i)
	{
		file << particles.y0[i] << " , ";
	}
	file << particles.y0[n-1] << " ]; " << std::endl << std::endl;

	if ( particles.nDim() == 3 )
	{
		file << " z = [ ";
		for ( int i = 0; i < n-1; ++i)
		{
			file << particles.z[i] << " , ";
		}
		file << particles.z[n-1] << " ]; " << std::endl << std::endl;
		
		file << " z0 = [ ";
		for ( int i = 0; i < n-1; ++i)
		{
			file << particles.z0[i] << " , ";
		}
		file << particles.z0[n-1] << " ]; " << std::endl << std::endl;
	}
	
	if ( _nDim == 1 )
	{
		file << "scalarField = [ " ;
		for ( int i = 0; i < n-1; ++i )
		{
			file << scalar[i] << " , ";
		}
		file << scalar[n-1] << " ]; " << std::endl;
	}
	else
	{
		if ( _nDim == 2)
		{
			file << "vectorField = [ ";
			for ( int i = 0; i < n-1; ++i)
			{
				file << xComp[i] << " , " << yComp[i] << " ; " ;
			}
			file << xComp[n-1] << " , " << yComp[n-1] << " ];" << std::endl;
		}
		else
		{
			file << "vectorField = [ ";
			for ( int i = 0; i < n-1; ++i)
			{
				file << xComp[i] << " , " << yComp[i] << " , " << zComp[i] << " ; " ;
			}
			file << xComp[n-1] << " , " << yComp[n-1] << " , " << zComp[n-1] << " ];" << std::endl;
		}
	}
};

void Field::writePointDataToVTK( std::ostream& fs ) const
{
	fs << "POINT_DATA " << _N << std::endl;
	fs << "SCALARS " << _name << " double " << _nDim << std::endl;
	fs << "LOOKUP_TABLE default " << std::endl;
	if ( _nDim == 1 )
	{
		for ( int i = 0; i < _N; ++i )
		{
			fs << scalar[i] << std::endl;
		}
	}
	else if ( _nDim == 2 )
	{
		for ( int i = 0; i < _N; ++i )
		{
			fs << xComp[i] << "     " << yComp[i] << 0.0 << std::endl;
		}
	}
	else
	{
		for ( int i = 0; i < _N; ++i )
		{
			fs << xComp[i] << "     " << yComp[i] << "     " << zComp[i] << std::endl;
		}
	}
};

void Field::writeCellDataToVTK( std::ostream& fs,  const Faces& faces ) const
{
	fs << "CELL_DATA " << faces.nActive() << std::endl;
	fs << "SCALARS " << _name << " double " << _nDim << std::endl;
	if ( _nDim == 1 )
	{
		for (int i = 0; i < faces.N(); ++i )
		{
			if ( faces.isActive(i) )
			{
				std::vector<int> faceVerts = faces.vertices(i);
				for ( int j = 0; j < faces.nVerts(i); ++j)
				{
					double subTriAvg = scalar[ faces.centerParticle(i) ];
					subTriAvg += scalar[ faceVerts[j] ] + scalar[ faceVerts[(j+1)%faces.nVerts(i)] ];
					subTriAvg /= 3.0;
					
					fs << subTriAvg << std::endl;
				}
			}
		}
	}
	else if ( _nDim == 2 )
	{
		for (int i = 0; i < faces.N(); ++i )
		{
			if ( faces.isActive(i) )
			{
				std::vector<int> faceVerts = faces.vertices(i);
				for ( int j = 0; j < faces.nVerts(i); ++j)
				{
					double subTriAvgX = xComp[ faces.centerParticle(i) ];
					double subTriAvgY = yComp[ faces.centerParticle(i) ];
					subTriAvgX += xComp[ faceVerts[j] ] + xComp[ faceVerts[(j+1)%faces.nVerts(i)] ];
					subTriAvgY += yComp[ faceVerts[j] ] + yComp[ faceVerts[(j+1)%faces.nVerts(i)] ];
					subTriAvgX /= 3.0;
					subTriAvgY /= 3.0;
					
					fs << subTriAvgX << "     " << subTriAvgY << std::endl;
				}
			}
		}
	}
	else
	{
		for (int i = 0; i < faces.N(); ++i )
		{
			if ( faces.isActive(i) )
			{
				std::vector<int> faceVerts = faces.vertices(i);
				for ( int j = 0; j < faces.nVerts(i); ++j)
				{
					double subTriAvgX = xComp[ faces.centerParticle(i) ];
					double subTriAvgY = yComp[ faces.centerParticle(i) ];
					double subTriAvgZ = zComp[ faces.centerParticle(i) ];
					subTriAvgX += xComp[ faceVerts[j] ] + xComp[ faceVerts[(j+1)%faces.nVerts(i)] ];
					subTriAvgY += yComp[ faceVerts[j] ] + yComp[ faceVerts[(j+1)%faces.nVerts(i)] ];
					subTriAvgZ += zComp[ faceVerts[j] ] + zComp[ faceVerts[(j+1)%faces.nVerts(i)] ];
					subTriAvgX /= 3.0;
					subTriAvgY /= 3.0;
					subTriAvgZ /= 3.0;
					
					fs << subTriAvgX << "     "<< subTriAvgY << "     "<< subTriAvgZ << std::endl;
				}
			}
		}
	}
};


/** @brief Releases unused memory, if too much was allocated.

	Call this method after a Particles object has been initialized, and only if no additional particles will be added.
*/
void Field::shrinkMemory()
{
	if ( _N < _nMax )
	{
		if ( _nDim == 1) 
		{
			for ( int i = _N+1; i < _nMax; ++i)
			{
				scalar.erase( scalar.cend() );
			}
			
			scalar.shrink_to_fit();
		}
		else
		{
			for ( int i = _N+1; i < _nMax; ++i )
			{
				xComp.erase( xComp.cend() );
				yComp.erase( yComp.cend() );				
			}	
			
			xComp.shrink_to_fit();
			yComp.shrink_to_fit();
			
			if ( _nDim == 3 )
			{
				for ( int i = _N+1; i < _nMax; ++i )
				{
					zComp.erase( zComp.cend() );
				}				
				zComp.shrink_to_fit();
			}
		}
		_nMax = _N;
	}	
};