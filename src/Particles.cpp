//
//  Particles.cpp
//  LPM-Testing
//
//  Created by Bosler, Peter Andrew on 11/3/14.
//  Copyright (c) 2014 Bosler, Peter Andrew. All rights reserved.
//

/** @file Particles.cpp
	@author Peter Bosler, Sandia National Laboratories
	@brief Particles class implementation
*/
#include "Particles.h"
#include <assert.h>
#include "OutputMessage.h"
#include "Logger.h"
#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "xyzVector.h"

/**
	allocates memory and initializes all particle coordinates, areas, and volume to 0.0.
*/
Particles::Particles( const int nDim, const int nMax, 
					  const int procRank, const int numProcs, const xyzVector::geometryKind gkin)
{
    assert( nDim == 2 || nDim == 3 );
    assert( nMax > 0 );
    assert( procRank >= 0);
    assert( numProcs >=1 );
    
    logLevel = OutputMessage::debugPriority;
    log = Logger::Instance(logLevel, procRank, numProcs);
    
//     OutputMessage statusMsg("constructing particles object", OutputMessage::debugPriority, "Particles constructor");
//     log->logMessage(statusMsg);
    
    _N = 0;
    _nDim = nDim;
    _nMax = nMax;
    switch (gkin) 
    {
    	case xyzVector::EuclideanGeometry :
    	{
    		gk = gkin;
    		break;
    	}
    	case xyzVector::SphericalGeometry :
    	{
    		gk = gkin;
    		break;
    	}
    	default :
    	{
    		OutputMessage errorMsg("invalid geometry kind", OutputMessage::errorPriority, "Particles constructor");
    		log->logMessage(errorMsg);
    	}
    }
    
    switch (nDim) {
        case 2:
            x = std::vector<double>(_nMax, 0.0);
            y = std::vector<double>(_nMax, 0.0);
            x0 = std::vector<double>(_nMax, 0.0);
            y0 = std::vector<double>(_nMax, 0.0);
            area = std::vector<double>(_nMax, 0.0);
            incidentEdges = std::vector< std::vector< std::pair<double,int> > > (_nMax);
            break;
        case 3:
            x = std::vector<double>(_nMax, 0.0);
            y = std::vector<double>(_nMax, 0.0);
            z = std::vector<double>(_nMax, 0.0);
            x0 = std::vector<double>(_nMax, 0.0);
            y0 = std::vector<double>(_nMax, 0.0);
            z0 = std::vector<double>(_nMax, 0.0);
            area = std::vector<double>(_nMax, 0.0);
            volume = std::vector<double>(_nMax, 0.0);
            incidentEdges = std::vector< std::vector< std::pair<double,int> > > (_nMax);
            break;
        default:
            OutputMessage errorMsg("invalid number of dimensions", OutputMessage::errorPriority, "Particles constructor");
            log->logMessage(errorMsg);
            break;
    }
    assert( _nMax == x.size());
};

void Particles::setPhysicalCoords( const int index, const xyzVector& physCoords )
{
    x[index] = physCoords.x;
    y[index] = physCoords.y;
    if ( _nDim == 3) 
	    z[index] = physCoords.z;
};

void Particles::setPhysicalCoords(const int index, const double nx, const double ny, const double nz)
{
    x[index] = nx;
    y[index] = ny;
    if ( _nDim == 3)
	    z[index] = nz;
};

xyzVector Particles::physCoord( const int index ) const
{
	if ( _nDim == 3 )
		return xyzVector( x[index], y[index], z[index] );
	else
		return xyzVector( x[index], y[index] );
};

xyzVector Particles::lagCoord( const int index ) const
{
	if ( _nDim == 3 )
		return xyzVector( x0[index], y0[index], z0[index] );
	else
		return xyzVector( x0[index], y0[index] );
}; 

void Particles::setLagrangianCoords( const int index, const xyzVector& lagCoords)
{
    x0[index] = lagCoords.x;
    y0[index] = lagCoords.y;
    if ( _nDim == 3)
	    z0[index] = lagCoords.z;
};

void Particles::setLagrangianCoords( const int index, const double ax, const double ay, const double az)
{
    x0[index] = ax;
    y0[index] = ay;
    if ( _nDim == 3)
    	 z0[index] = az;
};

void Particles::insertParticle( const xyzVector& physCoords, const xyzVector& lagCoords)
{
    assert(_N < _nMax);

	x[_N] = physCoords.x;
	y[_N] = physCoords.y;
	
    x0[_N] = lagCoords.x;
    y0[_N] = lagCoords.y;
    
    if ( _nDim == 3)
    {
    	z[_N] = physCoords.z;
	    z0[_N] = lagCoords.z;
	}
    _N += 1;
};

void Particles::writePointsToVTK( std::ostream& fs, std::string title ) const
{
	fs << "# vtk DataFile Version 2.0 " << std::endl;
	fs << title << std::endl;
	fs << "ASCII" << std::endl;
	fs << "DATASET POLYDATA" << std::endl;
	
	fs << "POINTS " << _N << " double " << std::endl;
	if ( _nDim == 2 )
	{
		for (int i = 0; i < _N; ++i )
		{
			fs << x[i] << "     " << y[i] << "     " << 0.0 << std::endl;
		}
	}
	else
	{
		for (int i = 0; i < _N; ++i )
		{
			fs << x[i] << "     " << y[i] << "     " << z[i] << std::endl;
		}
	}
};

void Particles::writeLagCoordToVTK( std::ostream& fs ) const
{
	fs << "POINT_DATA " << _N << std::endl;
	fs << "SCALARS lagParam double 3 " << std::endl;
	fs << "LOOKUP_TABLE default" << std::endl;
	if ( _nDim == 2 )
	{
		for ( int i = 0; i < _N; ++i)
		{
			fs << x0[i] << "     " << y0[i] << "     " << 0.0 << std::endl;
		}
	
	}
	else
	{
		for ( int i = 0; i < _N; ++i)
		{
			fs << x0[i] << "     " << y0[i] << "     " << z0[i] << std::endl;
		}
	}
};

void Particles::writeVariablesToMatlab( std::ostream& fs ) const 
{
	fs << "x = [";
	for ( int i = 0; i < _N-1; ++i)
	{
		fs << x[i] << " , ";
	}
	fs << x[_N-1] << " ];" << std::endl << std::endl;

	fs << "x0 = [";
	for ( int i = 0; i < _N-1; ++i)
	{
		fs << x0[i] << " , ";
	}
	fs << x0[_N-1] << " ];" << std::endl << std::endl;
	
	fs << "y = [";
	for ( int i = 0; i < _N-1; ++i)
	{
		fs << y[i] << " , ";
	}
	fs << y[_N-1] << " ];" << std::endl << std::endl;
	
	fs << "y0 = [";
	for ( int i = 0; i < _N-1; ++i)
	{
		fs << y0[i] << " , ";
	}
	fs << y0[_N-1] << " ];" << std::endl << std::endl;

	fs << "area = [";
	for ( int i = 0; i < _N-1; ++i )
	{
		fs << area[i] << " , " ;
	}	
	fs << area[_N-1] << " ];" << std::endl << std::endl;	
	
	if ( _nDim == 3)
	{
		fs << "z = [";
		for ( int i = 0; i < _N-1; ++i)
		{
			fs << z[i] << " , ";
		}
		fs << z[_N-1] << " ];" << std::endl << std::endl;

		fs << "z0 = [";
		for ( int i = 0; i < _N-1; ++i)
		{
			fs << z0[i] << " , ";
		}
		fs << z0[_N-1] << " ];" << std::endl << std::endl;
		
		fs << "volume= [";
		for ( int i = 0; i < _N-1; ++i )
		{
			fs << volume[i] << " , " ;
		}	
		fs << volume[_N-1] << " ];" << std::endl << std::endl;
	}
};

std::vector<int> Particles::edgesAtVertex( const int particleIndex ) const
{
	std::vector<int> result;
	for ( int j = 0; j < incidentEdges[particleIndex].size(); ++j )
	{
		result.push_back( incidentEdges[particleIndex][j].second );
	}
	return result;
};

void Particles::printIncidentEdges() const
{
	for (int i = 0; i < _N; ++i)
	{
		std::vector<double> vertAngles;
		std::vector<int> vertEdges;
		for ( int j = 0; j < incidentEdges[i].size(); ++j)
		{
			vertAngles.push_back( incidentEdges[i][j].first );
			vertEdges.push_back( incidentEdges[i][j].second );
		}
		
		if ( vertEdges.empty() )
		{
			std::cout << "i = " << i << ": no incident edges; area = " << area[i];
		}
		else
		{
			std::cout << "i = " << i << ", edges = ";
			for (int j = 0; j < vertEdges.size(); ++j )
			{
				std::cout << "( " << vertEdges[j] << ", " << vertAngles[j] << ")    ";
			}
		}
		std::cout << std::endl;
	}
};

double Particles::minX() const
{
    double result = x[0];
    for ( int i = 1; i < _N; ++i)
    {
        if ( x[i] < result )
            result = x[i];
    }
    return result;
};

double Particles::maxX() const
{
    double result = x[0];
    for ( int i = 1; i < _N; ++i)
    {
        if ( x[i] > result )
            result = x[i];
    }
    return result;
};

double Particles::minY() const
{
    double result = y[0];
    for ( int i = 1; i < _N; ++i)
    {
        if ( y[i] < result )
            result = y[i];
    }
    return result;
};

double Particles::maxY() const
{
    double result = y[0];
    for ( int i = 1; i < _N; ++i)
    {
        if ( y[i] > result )
            result = y[i];
    }
    return result;
};

double Particles::minZ() const
{
    double result = z[0];
    for ( int i = 1; i < _N; ++i)
    {
        if ( z[i] < result )
            result = z[i];
    }
    return result;
};

double Particles::maxZ() const
{
    double result = z[0];
    for ( int i = 1; i < _N; ++i)
    {
        if ( z[i] > result )
            result = z[i];
    }
    return result;
};

double Particles::minX0() const
{
    double result = x0[0];
    for ( int i = 1; i < _N; ++i)
    {
        if ( x0[i] < result )
            result = x0[i];
    }
    return result;
};

double Particles::maxX0() const
{
    double result = x0[0];
    for ( int i = 1; i < _N; ++i)
    {
        if ( x0[i] > result )
            result = x0[i];
    }
    return result;
};

double Particles::minY0() const
{
    double result = y0[0];
    for ( int i = 1; i < _N; ++i)
    {
        if ( y0[i] < result )
            result = y0[i];
    }
    return result;
};

double Particles::maxY0() const
{
    double result = y0[0];
    for ( int i = 1; i < _N; ++i)
    {
        if ( y0[i] > result )
            result = y0[i];
    }
    return result;
};

double Particles::minZ0() const
{
    double result = z0[0];
    for ( int i = 1; i < _N; ++i)
    {
        if ( z0[i] < result )
            result = z0[i];
    }
    return result;
};

double Particles::maxZ0() const
{
    double result = z0[0];
    for ( int i = 1; i < _N; ++i)
    {
        if ( z0[i] > result )
            result = z0[i];
    }
    return result;
};

std::vector< std::string > Particles::getInfo() const
{
	std::vector< std::string > info;
	std::stringstream ss;
	
	std::string s;
	
    ss << "     N = " << _N ;
    info.push_back( ss.str() );

    ss.str(s);
    ss << "  nMax = " << _nMax;  
    info.push_back( ss.str() );
    
    ss.str(s);
    ss << " size of x = " << x.size() ;
    info.push_back( ss.str() );
    
    ss.str(s);
    ss << "      min x = " << minX();
    info.push_back( ss.str() );
    
    ss.str(s);
    ss << "      max x = " << maxX();
    info.push_back( ss.str() );
    
    ss.str(s);
    ss << "   size of y = " << y.size();
    info.push_back( ss.str() );
    
    ss.str(s);
    ss << "      min y = " << minY();
    info.push_back( ss.str() );
    
    ss.str(s);
    ss << "      max y = " << maxY();
    info.push_back( ss.str() );
    
    ss.str(s);
    ss << "  size of z = " << z.size();
    info.push_back( ss.str() );
    
    if ( z.size() > 0 )
    {
    	ss.str(s);
    	ss << "    min z = " << minZ();
    	info.push_back( ss.str() );

		ss.str(s);
	    ss << "	   max z = " << maxZ();
	    info.push_back( ss.str() );
    }
	
	ss.str(s);ss.str(s);
	ss << " size of x0 = " << x0.size() ;
    info.push_back( ss.str() );
    
    ss.str(s);ss.str(s);ss.str(s);
    ss << "      min x0 = " << minX0();
    info.push_back( ss.str() );
    
	ss.str(s);
    ss << "      max x0 = " << maxX0();
    info.push_back( ss.str() );
    
    ss.str(s);
    ss << "   size of y0 = " << y0.size();
    info.push_back( ss.str() );

	ss.str(s);    
    ss << "      min y0 = " << minY0();
	info.push_back( ss.str() );
    
    ss.str(s);
    ss << "      max y0 = " << maxY0();
    info.push_back( ss.str() );
    
    ss.str(s);
    ss << "  size of z0 = " << z0.size();
    info.push_back( ss.str() );
    if ( z0.size() > 0 )
    {
    	ss.str(s);
    	ss << "    min z0 = " << minZ0();
    	info.push_back( ss.str() );
	    
	    ss.str(s);
	    ss << "	   max z0 = " << maxZ0();
	    info.push_back( ss.str() );
    }
    
//     ss.str(s);
//     for ( int i = 0; i < _N; ++i)
//     {
//     	ss << "particle " << i << ": (x,y) = (" << x[i] << ", " << y[i] << ")\n";
//     	ss << "            (x0,y0) = (" << x0[i] << ", " << y0[i] << ")\n";
//     	ss << "            area = " << area[i] << "\n";
//     	ss << "     incidentEdges = ";
//     	for ( int j = 0; j < incidentEdges[i].size(); ++j )
//     		ss << "(" << incidentEdges[i][j].first << ", " << incidentEdges[i][j].second << "), ";
//     	info.push_back(ss.str() );
//     	ss.str(s);
//     }
    
    return info;
}

/** @brief Releases unused memory, if too much was allocated.

	Call this method after a Particles object has been initialized, and only if no additional particles will be added.
*/
void Particles::shrinkMemory()
{
	if ( _N < _nMax )
	{
		for ( int i = _N+1; i < _nMax; ++i)
		{
			x.erase( x.cend() );
			y.erase( x.cend() );
			x0.erase( x.cend() );
			y0.erase( x.cend() );
			area.erase( x.cend() );
		}
		
		x.shrink_to_fit();
		y.shrink_to_fit();
		x0.shrink_to_fit();
		y0.shrink_to_fit();
	
		if ( _nDim == 3)
		{
			for ( int i = _N+1; i < _nMax; ++i)
			{
				z.erase( x.cend() );
				z0.erase( x.cend() );
				volume.erase( x.cend() );
			}
						
			z.shrink_to_fit();
			z0.shrink_to_fit();
			volume.shrink_to_fit();
		}
		_nMax = _N;
	}	
};

double Particles::totalArea() const 
{
	double result = 0.0;
	for ( int i = 0; i < _N; ++i )
		result += area[i];
	return result;
};

double Particles::totalVolume() const
{
	double result = 0.0;
	for ( int i = 0; i < _N; ++i)
		result += volume[i];
	return result;
};

void Particles::rescale( const double ampFactor )
{
	for ( int k = 0; k < _N; ++k )
	{
		x[k] *= ampFactor;
		y[k] *= ampFactor;
		x0[k] *= ampFactor;
		y0[k] *= ampFactor;
	}
	if ( _nDim == 3 )
	{
		for ( int k = 0; k < _N; ++k )
		{
			z[k] *= ampFactor;
			z0[k] *= ampFactor;
		}
	}
};


void Particles::sortIncidentEdgesAtParticle( const int particleIndex )
{
	if ( incidentEdges[particleIndex].size() == 0 )
	{
		OutputMessage statusMessage("no edges recorded at particle.", OutputMessage::errorPriority,
									"Particles::sortIncidentEdgesAtParticle");
		log->logMessage(statusMessage);
	}
	else if ( area[particleIndex] > 0.0 )
	{
		OutputMessage statusMessage("vertices should have zero area", OutputMessage::errorPriority,
									"Particles::sortIncidentEdgesAtParticle");
		log->logMessage(statusMessage);
	}
	else
	{
		std::sort( incidentEdges[particleIndex].begin(), incidentEdges[particleIndex].end() );
	}
};

std::ostream& operator << ( std::ostream& os, const Particles& aParticles)
{
	os << "particles N = " << aParticles.N() << std::endl;
	os << "particles nMax = " << aParticles.nMax() << std::endl;
	
	return os;
};