//
//  Edges.cpp
//  LPM
//
//  Created by Peter Bosler on 11/5/14.
//  Copyright (c) 2014 Peter Bosler. All rights reserved.
//
/** @file Edges.cpp 
	@brief Edges class implementation
	@author Peter Bosler, Sandia National Laboratories, Multiphysics Applications
*/

#include "Edges.h"
#include <vector>
#include <assert.h>
#include "GlobalConstants.h"
#include "OutputMessage.h"
#include "Logger.h"
#include "xyzVector.h"
#include "Particles.h"
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>

Edges::Edges( const int nMax, const int procRank, const int numProcs, const xyzVector::geometryKind gkin)
{
	assert( nMax > 0 );
	assert( procRank >= 0 );
	assert( numProcs >= 1 );
	
	logLevel = OutputMessage::debugPriority;
	log = Logger::Instance( logLevel, procRank, numProcs );
	
	_N = 0;
	_nMax = nMax;
	gk = gkin;
	
	_orig = std::vector<int>(nMax, -1);
	_dest = std::vector<int>(nMax, -1);
	_rightFace = std::vector<int>(nMax, -1);
	_leftFace = std::vector<int>( nMax, -1);
	_parent = std::vector<int>(nMax, -1);
	_hasChildren = std::vector<bool>(nMax, false);
	_child1 = std::vector<int>(nMax, -1);
	_child2 = std::vector<int>(nMax, -1);
};

void Edges::insertEdge( const int origIndex, const int destIndex, 
						const int leftFace, const int rightFace, Particles& particles )
{
	assert( _N < _nMax );
	
	_orig[_N] = origIndex;
	_dest[_N] = destIndex;
	_leftFace[_N] = leftFace;
	_rightFace[_N] = rightFace;
	
	recordIncidentEdgeAtParticles(_N, particles);

	_N+=1;
};

void Edges::recordIncidentEdgeAtParticles( const int edgeIndex, Particles& particles ) const
{	
	if ( _hasChildren[edgeIndex] )
		return;

	if ( particles.nDim() == 2 && gk == xyzVector::EuclideanGeometry )
	{
		// Record edge at origin vertex
		bool duplicateEdge = false;
		for ( int j =0; j < particles.incidentEdges[ _orig[edgeIndex]].size(); ++j )
		{
			if ( particles.incidentEdges[_orig[edgeIndex]][j].second == edgeIndex )
				duplicateEdge = true;
		}
		double angleVal;
		// define angle relative to positive real axis (planar geometry only)
		if ( !duplicateEdge )
		{
			angleVal = std::atan2( particles.y[_dest[edgeIndex]] - particles.y[_orig[edgeIndex]], 
										  particles.x[_dest[edgeIndex]] - particles.x[_orig[edgeIndex]]);
			particles.incidentEdges[_orig[edgeIndex]].push_back(std::pair<double,int>(angleVal,edgeIndex));
		}
		duplicateEdge = false;
		for ( int j = 0; j < particles.incidentEdges[_dest[edgeIndex]].size(); ++j)
		{
			if ( particles.incidentEdges[_dest[edgeIndex]][j].second == edgeIndex )
				duplicateEdge = true;
		}
		if ( !duplicateEdge )
		{
			angleVal = std::atan2( particles.y[_orig[edgeIndex]] - particles.y[_dest[edgeIndex]],
								   particles.x[_orig[edgeIndex]] - particles.x[_dest[edgeIndex]]);
			particles.incidentEdges[_dest[edgeIndex]].push_back(std::pair<double,int>( angleVal, edgeIndex));
		}
	}
	else if ( particles.nDim() == 3 )
	{
		// record edge at origin vertex
		bool duplicateEdge = false;
		for ( int j = 0; j < particles.incidentEdges[_orig[edgeIndex] ].size(); ++j )
		{
			if ( particles.incidentEdges[ _orig[edgeIndex] ][j].second == edgeIndex )
				duplicateEdge = true;
		}	
		double angleVal;
		// define angle relative to first edge at vertex
		if ( !duplicateEdge )
		{
			if ( particles.incidentEdges[ _orig[edgeIndex] ].empty() )
				particles.incidentEdges[ _orig[edgeIndex] ].push_back( std::pair<double, int>( 0.0, edgeIndex ) );
			else
			{
				angleVal = edgeAngleAtOrig(edgeIndex, particles);
				particles.incidentEdges[_orig[edgeIndex]].push_back( std::pair<double,int>(angleVal,edgeIndex));
			}
		}
		
		// record edge at destination vertex
		duplicateEdge = false;
		for ( int j = 0; j < particles.incidentEdges[_dest[edgeIndex]].size(); ++j)
		{
			if ( particles.incidentEdges[_dest[edgeIndex]][j].second == edgeIndex )
				duplicateEdge = true;
		}
		if ( !duplicateEdge )
		{
			if ( particles.incidentEdges[_dest[edgeIndex]].empty() )
				particles.incidentEdges[_dest[edgeIndex]].push_back( std::pair<double,int>(0.0, edgeIndex ) );
			else
			{		
				angleVal = edgeAngleAtDest(edgeIndex,particles);			
				particles.incidentEdges[_dest[edgeIndex]].push_back( std::pair<double,int>(angleVal, edgeIndex ) );
			}
		}
	}
};

double Edges::edgeAngleAtOrig( const int edgeIndex, const Particles& particles ) const
{
	GlobalConstants* consts = GlobalConstants::Instance();
	const double PI = consts->Pi();
	double angleVal;
	xyzVector v0;
	xyzVector v1;
	if ( particles.incidentEdges[ _orig[edgeIndex] ].empty() )
		angleVal = 0.0;			
	else
	{
		if ( _orig[ particles.incidentEdges[ _orig[edgeIndex]][0].second ] == _orig[edgeIndex] )
		{
			v0 = particles.physCoord( _orig[particles.incidentEdges[_orig[edgeIndex]][0].second] );
			v1 = particles.physCoord( _dest[particles.incidentEdges[_orig[edgeIndex]][0].second] );
		}
		else if ( _dest[particles.incidentEdges[_orig[edgeIndex]][0].second] == _orig[edgeIndex] )
		{
			v0 = particles.physCoord( _dest[particles.incidentEdges[_orig[edgeIndex]][0].second] );
			v1 = particles.physCoord( _orig[particles.incidentEdges[_orig[edgeIndex]][0].second] );
		}
		else
		{
			std::stringstream ss;
			std::string s;
			
			ss << "connectivity error at new edge origin ; index = " << edgeIndex;
			
			OutputMessage statusMsg(ss.str(), OutputMessage::errorPriority, "Edges::edgeAngleAtOrig");
			log->logMessage(statusMsg);
			ss.str(s);
		}
		
		xyzVector e0 = v1 - v0;
		
		xyzVector newO = particles.physCoord( _orig[edgeIndex] );
		xyzVector newD = particles.physCoord( _dest[edgeIndex] );
		xyzVector newE = newD - newO;
		
		e0.normalize();
		newE.normalize();
		xyzVector cp = e0.crossProduct( newE );		
		
		angleVal = atan4( cp.magnitude(), e0.dotProduct(newE) );
		if ( e0.dotProduct(newE) < 0.0 )
			angleVal = 2.0 * PI - angleVal;
	}
	
	return angleVal;
};

double Edges::edgeAngleAtDest( const int edgeIndex, const Particles& particles ) const
{
	GlobalConstants* consts = GlobalConstants::Instance();
	const double PI = consts->Pi();
	double angleVal;
	xyzVector v0; 
	xyzVector v1;
	if ( particles.incidentEdges[_dest[edgeIndex]].empty() )
		angleVal = 0.0;
	else
	{	
		if ( _orig[ particles.incidentEdges[ _dest[edgeIndex]][0].second ] == _dest[edgeIndex] )
		{
			v0 = particles.physCoord( _orig[particles.incidentEdges[_dest[edgeIndex]][0].second] );
			v1 = particles.physCoord( _dest[particles.incidentEdges[_dest[edgeIndex]][0].second] );
		}
		else if ( _dest[particles.incidentEdges[_dest[edgeIndex]][0].second] == _dest[edgeIndex] )
		{
			v0 = particles.physCoord( _dest[particles.incidentEdges[_dest[edgeIndex]][0].second] );
			v1 = particles.physCoord( _orig[particles.incidentEdges[_dest[edgeIndex]][0].second] );
		}
		else
		{
			std::stringstream ss;
			std::string s;
			ss << "connectivity error at new edge destination, index = " << edgeIndex << ", orig = "
				<< _orig[edgeIndex] << ", dest = " << _dest[edgeIndex] << ", edges at dest particle = ";
				for ( int i = 0; i < particles.incidentEdges[_dest[edgeIndex]].size(); ++i )
				{
					ss << particles.incidentEdges[_dest[edgeIndex]][i].second << ", ";
				}
				ss << std::endl;
			OutputMessage statusMsg( ss.str(), OutputMessage::errorPriority, 
									"Edges::edgeAngleAtDest");
			log->logMessage(statusMsg);
			ss.str(s);
		}
		
		xyzVector e0 = v1 - v0;
		xyzVector newO = particles.physCoord( _orig[edgeIndex] );
		xyzVector newD = particles.physCoord( _dest[edgeIndex] );
		xyzVector newE = newO - newD;
		
		e0.normalize();
		newE.normalize();
		xyzVector cp = e0.crossProduct( newE );
		
		angleVal = atan4( cp.magnitude(), e0.dotProduct(newE) );
		if ( e0.dotProduct(newE) < 0.0 )
			angleVal = 2.0 * PI - angleVal;
	}
	
	return angleVal;
};

void Edges::divideEdge( const int index, Particles& particles )
{
	assert( _N+2 <= _nMax );

	xyzVector midpt;
	xyzVector lagMidpt;
	if ( gk == xyzVector::EuclideanGeometry )
	{
		midpt = cartMidpoint( particles.physCoord( _orig[index] ), particles.physCoord( _dest[index] ) );
		lagMidpt = cartMidpoint( particles.lagCoord( _orig[index] ), particles.lagCoord( _dest[index] ) );
	}
	else if ( gk == xyzVector::SphericalGeometry )
	{
		midpt = sphereMidpoint(particles.physCoord( _orig[index] ), particles.physCoord( _dest[index] ) );
		lagMidpt = sphereMidpoint(particles.lagCoord( _orig[index]), particles.lagCoord(_dest[index]) );
	}
	
	const int pinsert = particles.N();
	particles.insertParticle( midpt, lagMidpt );
	
	_hasChildren[index] = true;
	_child1[index] = _N;
	_child2[index] = _N+1;
	_parent[_N] = index;
	_parent[_N+1] = index;
	
	_orig[_N] = _orig[index];
	_dest[_N] = pinsert;
	_leftFace[_N] = _leftFace[index];
	_rightFace[_N] = _rightFace[index];
	
	_orig[_N+1] = pinsert;
	_dest[_N+1] = _dest[index];
	_leftFace[_N+1] = _leftFace[index];
	_rightFace[_N+1] = _rightFace[index];
	
	replaceIncidentEdgeWithChild( index, _N, particles);
	recordIncidentEdgeAtParticles(_N, particles );
	recordIncidentEdgeAtParticles(_N+1, particles );
	_N += 2;
};

void Edges::replaceIncidentEdgeWithChild( const int parentIndex, const int child1Index, Particles& particles ) const
{
	// replace edge at origin vertex
	int edgeIndexAtVertex = -1;
	for ( int i = 0; i < particles.incidentEdges[ _orig[parentIndex] ].size(); ++i)
	{
		if ( particles.incidentEdges[_orig[parentIndex]][i].second == parentIndex )
		{
			edgeIndexAtVertex = i;
			break;
		}
	}
	if ( edgeIndexAtVertex == -1 )
	{
		OutputMessage statusMsg("parent edge not found at origin", OutputMessage::errorPriority,
			 "Edges::replaceIncidentEdgeWithChild");
		log->logMessage(statusMsg);
	}
	else
	{
		particles.incidentEdges[_orig[parentIndex]][edgeIndexAtVertex].second = child1Index;	
	}
	
	// replace edge at destination vertex
	edgeIndexAtVertex = -1;
	for ( int i = 0; i < particles.incidentEdges[_dest[parentIndex]].size(); ++i)
	{
		if ( particles.incidentEdges[_dest[parentIndex]][i].second == parentIndex )
		{
			edgeIndexAtVertex = i;
			break;
		}
	}
	if ( edgeIndexAtVertex == -1 )
	{
		OutputMessage statusMsg("parent edge not found at destination", OutputMessage::errorPriority,
			 "Edges::replaceIncidentEdgeWithChild");
		log->logMessage(statusMsg);
	}
	else
	{
		particles.incidentEdges[_dest[parentIndex]][edgeIndexAtVertex].second = child1Index+1;	
	}
};

void Edges::writeVariablesToMatlab( std::ostream& fs ) const
{
	fs << "edgeN = " << _N << ";" << std::endl;
	fs << "edgeNMax = " << _nMax << ";" <<std::endl;
	
	fs << "edgeOrig = [";
	for ( int i = 0; i < _N-1; ++i)
	{
		fs << _orig[i] << " , ";
	}
	fs << _orig[_N-1] << " ];" << std::endl << std::endl;
	
	fs << "edgeDest = [";
	for ( int i = 0; i < _N-1; ++i)
	{
		fs << _dest[i] << " , ";
	}
	fs << _dest[_N-1] << " ];" << std::endl << std::endl;
	
	fs << "edgeLeft = [";
	for ( int i = 0; i < _N-1; ++i)
	{
		fs << _leftFace[i] << " , ";
	}
	fs << _leftFace[_N-1] << " ];" << std::endl << std::endl;
	
	fs << "edgeRight = [";
	for ( int i = 0; i < _N-1; ++i)
	{
		fs << _rightFace[i] << " , ";
	}
	fs << _rightFace[_N-1] << " ];" << std::endl << std::endl;

	fs << "edgeHasChildren = [" ;
	for ( int i = 0; i < _N-1; ++i )
	{
		if ( _hasChildren[i] )
			fs << "true" << " , ";
		else
			fs << "false" << " , ";
	}
	if ( _hasChildren[_N-1] )
		fs << "true" << "];" << std::endl << std::endl;
	else
		fs << "false" << "];" << std::endl << std::endl;
		
	fs << "edgeChildren = [";
	for ( int i = 0; i < _N-1; ++i)
	{
		fs << _child1[i] << " , " << _child2[i] << " ; " ;
	}	
	fs << _child1[_N-1] << " , " << _child2[_N-1] << "];" << std::endl << std::endl;
	
	fs << "edgeParent = [";
	for ( int i = 0; i < _N-1; ++i )
	{
		fs << _parent[i] << " , ";
	}	
	fs << _parent[_N-1] << "];" << std::endl << std::endl;
};

std::ostream& operator << ( std::ostream& os, const Edges& edges)
{
	os << "Edges info : \n ";
	os << "    N = " << edges.N() << std::endl;
	os << "   nMax = " << edges.nMax() << std::endl;
	return os;
};

std::vector< std::string > Edges::getInfo() const
{
	std::vector< std::string > info;
	std::stringstream ss;
	std::string s;
	
	ss << "    N = " << _N;
	info.push_back( ss.str() );
	ss.str(s);
	
	ss << "    nMax = " << _nMax;
	info.push_back( ss.str() );
	ss.str(s);
	
	for ( int i = 0; i < _N; ++i)
	{
		ss << "edge " << i << ": orig = " << _orig[i] << ", dest = " << _dest[i] 
			<< ", left = " << _leftFace[i] << ", right = " << _rightFace[i] ;
		if ( _hasChildren[i] )
		{
			ss << ", DIVIDED, parent = " << _parent[i] 
				<< ", children = " << _child1[i] << ", " << _child2[i] ;
		}
		else 
		{
			ss << ", not divided, parent = "<< _parent[i];
		}
		
		info.push_back( ss.str() );
		ss.str(s);
	}
	
	return info;
};

void Edges::shrinkMemory()
{
	if ( _N < _nMax )
	{
		for ( int i = _N+1; i < _nMax; ++i)
		{
			_orig.erase( _orig.cend() );
			_dest.erase( _dest.cend() );
			_rightFace.erase( _rightFace.cend() );
			_leftFace.erase( _leftFace.cend() );
			_parent.erase( _parent.cend() );
			_child1.erase( _child1.cend() );
			_child2.erase( _child2.cend() );
			_hasChildren.erase( _hasChildren.cend() );
		}
		_orig.shrink_to_fit();
		_dest.shrink_to_fit();
		_rightFace.shrink_to_fit();
		_leftFace.shrink_to_fit();
		_parent.shrink_to_fit();
		_hasChildren.shrink_to_fit();
		_child1.shrink_to_fit();
		_child2.shrink_to_fit();
		_nMax = _N;
	}	
};