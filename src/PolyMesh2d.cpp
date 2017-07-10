//
//  PolyMesh2d.cpp
//  LPM
//
//  Created by Peter Bosler on 11/5/14.
//  Copyright (c) 2014 Peter Bosler. All rights reserved.
//

/** @file PolyMesh2d.cpp
	@brief PolyMesh2d class implementation
	@author Peter Bosler, Sandia National Laboratories, Multiphysics Applications
*/

#include "PolyMesh2d.h"
#include "Particles.h"
#include "Edges.h"
#include "Faces.h"
#include "OutputMessage.h"
#include "Logger.h"
#include "Field.h"
#include "xyzVector.h"
#include <iostream>
#include <vector>
#include <assert.h>
#include <cmath>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <limits>
#include <iomanip>

PolyMesh2d::PolyMesh2d( const int initNest, const meshSeed seed,
					const double maxR, const int procRank, const int numProcs )
{
	assert(initNest >= 0 );
	assert(maxR > 0.0 );
	assert(procRank >= 0 );
	assert(numProcs >= 1 );
	
	logLevel = OutputMessage::debugPriority;
	log = Logger::Instance( logLevel, procRank, numProcs );
	
	_initNest = initNest;
	_maxR = maxR;
	_seed = seed;
	
	OutputMessage statusMsg;
	
	switch ( seed )
	{
		case triHexSeed :
		{
			_faceKind = Faces::triangularFaces;
			gk = xyzVector::EuclideanGeometry;
			_nCoordDim = 2;
			break;
		}
		case quadRectSeed :
		{
			_faceKind = Faces::quadrilateralFaces;
			gk = xyzVector::EuclideanGeometry;
			_nCoordDim = 2;
			break;
		}
		case polarDiscSeed :
		{
			_nCoordDim = 2;
			statusMsg = OutputMessage("polar discs not implemented yet.", OutputMessage::errorPriority,
													  "PolyMesh2d constructor");
			log->logMessage(statusMsg);
			break;
		}
		case quadRectPeriodic :
		{	
			_nCoordDim = 2;
			statusMsg = OutputMessage("periodic boundaries not implemented yet.", OutputMessage::errorPriority,
													  "PolyMesh2d constructor");
			log->logMessage(statusMsg);
			break;
		}
		case icosTriSphereSeed :
		{
			_faceKind = Faces::triangularFaces;
			gk = xyzVector::SphericalGeometry;
			_nCoordDim = 3;
			break;
		}
		case cubedSphereSeed :
		{
			_faceKind = Faces::quadrilateralFaces;
			gk = xyzVector::SphericalGeometry;
			_nCoordDim = 3;
			break;
		}
		default :
		{
			statusMsg = OutputMessage("invalid meshSeed", OutputMessage::errorPriority, "PoilyMesh2d constructor");
			log->logMessage(statusMsg);
		}
	} // switch seed
	
	int nMaxParticles = 0;
	int nMaxFaces = 0;
	int nMaxEdges = 0;
	nMaxParticles = nVertices(initNest) + nFaces(initNest);
	std::cout << "nMaxParticles = " << nMaxParticles << std::endl;
	for ( int k = 0; k <= initNest; ++k )
	{
		nMaxFaces += nFaces(k);
		nMaxEdges += nEdges( nVertices(k), nFaces(k) );
	}
	std::cout << "nMaxFaces = " << nMaxFaces << std::endl;
	std::cout << "nMaxEdges = " << nMaxEdges << std::endl;
	particles = Particles( _nCoordDim, nMaxParticles, procRank, numProcs, gk );
	edges = Edges( nMaxEdges, procRank, numProcs, gk );
	faces = Faces( nMaxFaces, _faceKind, procRank, numProcs, gk );
	
	initializeFromSeed( _maxR );
	int startIndex = 0;
	int nOldFaces = 0;
	for ( int k = 1; k <= initNest; ++k )
	{
		nOldFaces = faces.N();
		for ( int i = startIndex; i < nOldFaces; ++i )
		{
			if ( _faceKind == Faces::triangularFaces )
			{
				faces.divideTriFace( i, particles, edges );
			}
			else if ( _faceKind == Faces::quadrilateralFaces )
			{
				faces.divideQuadFace( i, particles, edges );
			}
		}
		startIndex = nOldFaces;
	}
	
	sortEdgesAroundVertices();
};

void PolyMesh2d::shrinkMemory()
{
	particles.shrinkMemory();
	edges.shrinkMemory();
	faces.shrinkMemory();
};

xyzVector PolyMesh2d::facePosition( const int faceIndex ) const
{
	return particles.physCoord( faces.centerParticle( faceIndex ) );
};

xyzVector PolyMesh2d::faceCentroid( const int faceIndex ) const
{
	xyzVector result = xyzVector();
	std::vector<int> vertIndices = faces.vertices(faceIndex);
	std::vector< xyzVector > faceVerts(vertIndices.size() );
	for (int i = 0; i < vertIndices.size(); ++i )
	{
		faceVerts[i] = particles.physCoord(vertIndices[i]);
	}
	if ( gk == xyzVector::EuclideanGeometry )
		result = cartCentroid( faceVerts );
	else if ( gk == xyzVector::SphericalGeometry )
		result = sphereCentroid( faceVerts );
	
	return result;
};

void PolyMesh2d::setSurfaceAreaFromParticles()
{
	_surfaceArea = particles.totalArea();
};

void PolyMesh2d::setSurfaceAreaByCalculation()
{
	_surfaceArea = 0.0;
	for ( int i = 0; i < faces.N(); ++i )
	{
		if ( faces.isActive(i) )
		{
			faces.setArea(i, particles, edges );
			_surfaceArea += faces.area( i, particles );
		}
	}
};

int PolyMesh2d::locatePointInMeshFace( const xyzVector& queryPt ) const
{
	int treeStart = nearestRootFace(queryPt);
	int walkStart = locatePointTreeSearch( queryPt, treeStart );
	return locatePointWalkSearch( queryPt, walkStart );
};

int PolyMesh2d::nearestRootFace( const xyzVector& queryPt ) const
{
	int result = 0;
	double dist = distance( faceCentroid(0), queryPt, gk );
	double testDist = dist;
	int nRootFaces;
	switch ( _seed )
	{
		case triHexSeed :
		{
			nRootFaces = 6;
			break;
		}
		case quadRectSeed :
		{
			nRootFaces = 4;
			break;
		}
		case polarDiscSeed :
		{
			nRootFaces = 4;
			break;
		}
		case quadRectPeriodic :
		{
			nRootFaces = 4; 
			break;
		}
		case icosTriSphereSeed :
		{
			nRootFaces = 20;
			break;
		}
		case cubedSphereSeed :
		{
			nRootFaces = 6;
			break;
		}
	}
	for ( int i = 1; i < nRootFaces; ++i )
	{
		testDist = distance( faceCentroid(i), queryPt, gk );
		if ( testDist < dist )
		{
			result = i;
			dist = testDist;
		}
	}
	return result;
};

int PolyMesh2d::locatePointTreeSearch( const xyzVector& queryPt, const int startIndex ) const
{
	double currentMinDist = 1.0e20;
	double testDist;
	int nearestChild = -1;
	if ( faces.hasChildren( startIndex ) )
	{
		std::vector<int> childFaces = faces.children(startIndex);
		xyzVector cntd;
		for ( int i = 0; i < 4; ++i )
		{
			cntd = faceCentroid( childFaces[i] );
			testDist = distance( queryPt, cntd, gk );
			if ( testDist < currentMinDist )
			{
				currentMinDist = testDist;
				nearestChild = childFaces[i];
			}
		}
		return locatePointTreeSearch( queryPt, nearestChild );
	}
	else
	{
		return startIndex;
	}
};

bool PolyMesh2d::pointIsOutsideMesh( const xyzVector& queryPt ) const
{
	bool result = false;
	int faceInd = locatePointInMeshFace(queryPt);
	std::vector<int> faceEdges = ccwEdgesAroundFace( faceInd );
	std::vector<int> bndyEdges;
	for ( int i = 0; i < faceEdges.size(); ++i )
	{
		if ( edges.onBoundary( faceEdges[i] ) )
			bndyEdges.push_back(faceEdges[i]) ;
	}
	if ( bndyEdges.size() > 0 )
	{
		xyzVector cntd = faceCentroid( faceInd );
		double interiorDist = cartDistance( queryPt, cntd );
		for ( int i = 0; i < bndyEdges.size(); ++i )
		{
			xyzVector v1 = particles.physCoord( edges.orig(bndyEdges[i]) );
			xyzVector v2 = particles.physCoord( edges.dest(bndyEdges[i]) );
			xyzVector q = v2 - v1;
			q.normalize();
			xyzVector p = cntd - v1;
			double scalarProj = p.dotProduct(q);
			q.scale(2.0 * scalarProj);
			xyzVector reflection = q - p  + v1;
			
			double extDist = cartDistance( queryPt, reflection );
			
			if ( extDist < interiorDist ) 
				result = true;
		}
	}
	return result;
};

std::vector<int> PolyMesh2d::ccwAdjacentFaces( const int faceIndex ) const
{	
	std::vector<int> result;
	result.reserve(8);
	std::vector<int> ccwEdges = ccwEdgesAroundFace( faceIndex );
	if ( faces.hasChildren(faceIndex) )
	{
		OutputMessage statusMsg("returning ccw adjacent faces for a divided face.", OutputMessage::warningPriority,
								"PolyMesh2d::ccwAdjacentFaces");
		log->logMessage(statusMsg);
	}
	for ( int i = 0; i < ccwEdges.size(); ++i )
	{
		if ( edges.leftFace( ccwEdges[i] ) == faceIndex )
			result.push_back( edges.rightFace( ccwEdges[i] ) );
		else if ( edges.rightFace( ccwEdges[i] ) == faceIndex )
			result.push_back( edges.leftFace( ccwEdges[i] ) );
		else
		{
			std::vector<std::string> info;
			std::stringstream ss;
			std::string s;
			std::vector<int> faceChildren = faces.children(faceIndex);
			ss << "face " << faceIndex << ": parent = " << faces.parent(faceIndex) << ", children = ("
					<< faceChildren[0] << ", " << faceChildren[1] << ", " << faceChildren[2] << ", " 
					<< faceChildren[3] << "); ccwEdges = ";
				for ( int j = 0; j < ccwEdges.size(); ++j )
					ss << ccwEdges[j] << ", ";
				info.push_back( ss.str() );
				ss.str(s);
				ss << "face " << faceIndex << "; edge " << ccwEdges[i] << ": orig = " << edges.orig(ccwEdges[i])
					<< ", dest = " << edges.dest(ccwEdges[i]) << ", leftFace = " << edges.leftFace(ccwEdges[i])
					<< ", rightFace = " << edges.rightFace(ccwEdges[i]) << ", parent = " << edges.parent(ccwEdges[i])
					<< ", children = " << edges.child1(ccwEdges[i]) << ", " << edges.child2(ccwEdges[i]);
				info.push_back(ss.str());
				LongMessage statusMsg("connectivity error; edge does not connect to face.", 
					OutputMessage::errorPriority, "PolyMesh2d::ccwAdjacentFaces",info);
				log->logMessage(statusMsg);	
		}
	}
	return result;
};

std::vector<int> PolyMesh2d::ccwFacesAroundVertex( const int particleIndex) const
{	
	std::vector<int> result;
	if ( particles.area[particleIndex] > 0.0 )
	{
		OutputMessage statusMsg("vertices should have area = 0.0", OutputMessage::errorPriority,
								"PolyMesh2d::ccwFacesAroundVertex");
		log->logMessage(statusMsg);
		return result;
	}
	result.reserve(6);
	for ( int i = 0; i < particles.incidentEdges[particleIndex].size(); ++i )
	{
		if ( particleIndex == edges.orig( particles.incidentEdges[particleIndex][i].second ) )
			result.push_back( edges.leftFace( particles.incidentEdges[particleIndex][i].second ) );
		else if ( particleIndex == edges.dest( particles.incidentEdges[particleIndex][i].second ) )
			result.push_back( edges.rightFace( particles.incidentEdges[particleIndex][i].second) );
	}
	return result;	
};

void PolyMesh2d::initializeFromSeed( const double ampFactor )
{
	OutputMessage statusMsg;
	GlobalConstants* physConsts = GlobalConstants::Instance();
	
	switch ( _seed )
	{
		case PolyMesh2d::triHexSeed :
		{
			std::vector< xyzVector > coords(13);
			std::vector<int> edgeOrigs(12,-1) ;
			std::vector<int> edgeDests(12,-1) ;
			std::vector<int> edgeLefts(12,-1) ;
			std::vector<int> edgeRghts(12,-1) ;
			std::vector< std::vector<int> > triFaceVerts(6, std::vector<int>(3, -1) );
			std::vector< std::vector<int> > triFaceEdges(6, std::vector<int>(3, -1) );
			std::vector< std::vector< std::pair<double,int> > > incidentEdges(7);
	
			coords[0] = xyzVector( 0.0, 0.0 );
			coords[1] = xyzVector( 0.5, 0.866025403784438597 );
			coords[2] = xyzVector(-0.5, 0.866025403784438597 );
			coords[3] = xyzVector(-1.0, 0.0 );
			coords[4] = xyzVector(-0.5, -0.866025403784438597 );
			coords[5] = xyzVector( 0.5, -0.866025403784438597 );
			coords[6] = xyzVector( 1.0, 0.0 );
			coords[7] = xyzVector( 0.0, 0.577350269189625731  );
			coords[8] = xyzVector(-0.5, 0.288675134594812810 );
			coords[9] = xyzVector(-0.5, -0.288675134594812810 );
			coords[10]= xyzVector( 0.0, -0.577350269189625731 );
			coords[11]= xyzVector( 0.5, -0.288675134594812810 );
			coords[12]= xyzVector( 0.5, 0.288675134594812810 );
			for ( int k = 0; k < 13; ++k)
			{
				coords[k].scale(ampFactor);
				particles.insertParticle( coords[k], coords[k] );
			}
			assert( particles.N() == 13 );
			
			const double PI = physConsts->Pi();
			
			incidentEdges[0].push_back( std::pair<double,int>( -2.0*PI/3.0, 9) );
			incidentEdges[0].push_back( std::pair<double,int>( -PI/3.0, 10) );
			incidentEdges[0].push_back( std::pair<double,int>( 0.0, 11) );
			incidentEdges[0].push_back( std::pair<double,int>( PI/3.0, 0) );
			incidentEdges[0].push_back( std::pair<double,int>( 2.0 * PI/3.0, 7 ) );
			incidentEdges[0].push_back( std::pair<double,int>( PI, 8) );
			
			incidentEdges[1].push_back( std::pair<double,int>( -2.0*PI/3.0, 0) );
			incidentEdges[1].push_back( std::pair<double,int>( -PI/3.0, 6) );
			incidentEdges[1].push_back( std::pair<double,int>( PI, 1) );
			
			incidentEdges[2].push_back( std::pair<double,int>( -2.0*PI/3.0, 2) );
			incidentEdges[2].push_back( std::pair<double,int>( - PI/3.0, 7) );
			incidentEdges[2].push_back( std::pair<double,int>( 0.0, 1) );
			
			incidentEdges[3].push_back( std::pair<double,int>( -PI/3.0, 3) );
			incidentEdges[3].push_back( std::pair<double,int>( 0.0, 8) );
			incidentEdges[3].push_back( std::pair<double,int>( PI/3.0, 2) );
			
			incidentEdges[4].push_back( std::pair<double,int>( 0.0, 4) );
			incidentEdges[4].push_back( std::pair<double,int>( PI/3.0, 9) );
			incidentEdges[4].push_back( std::pair<double,int>( 2.0 * PI/3.0, 3) );
			
			incidentEdges[5].push_back( std::pair<double,int>( PI/3.0, 5) );
			incidentEdges[5].push_back( std::pair<double,int>( 2.0 * PI/3.0, 10) );
			incidentEdges[5].push_back( std::pair<double,int>( PI , 4) );
			
			incidentEdges[6].push_back( std::pair<double,int>( -2.0*PI/3.0 , 5) );
			incidentEdges[6].push_back( std::pair<double,int>( 2.0*PI/3.0 , 6) );
			incidentEdges[6].push_back( std::pair<double,int>( PI , 11) );	
					
			for ( int k = 0; k < 7; ++k )
			{
				particles.incidentEdges[k] = incidentEdges[k];
			}
			
			edgeOrigs[0] = 0;
			edgeDests[0] = 1;
			edgeLefts[0] = 0;
			edgeRghts[0] = 5;
			
			edgeOrigs[1] = 1;
			edgeDests[1] = 2;
			edgeLefts[1] = 0;
			
			edgeOrigs[2] = 2;
			edgeDests[2] = 3;
			edgeLefts[2] = 1;
			
			edgeOrigs[3] = 3;
			edgeDests[3] = 4;
			edgeLefts[3] = 2;
			
			edgeOrigs[4] = 4;
			edgeDests[4] = 5;
			edgeLefts[4] = 3;
			
			edgeOrigs[5] = 5;
			edgeDests[5] = 6;
			edgeLefts[5] = 4;

			edgeOrigs[6] = 6;
			edgeDests[6] = 1;
			edgeLefts[6] = 5;			
			
			edgeOrigs[7] = 2;
			edgeDests[7] = 0;
			edgeLefts[7] = 0;
			edgeRghts[7] = 1;
			
			edgeOrigs[8] = 3;
			edgeDests[8] = 0;
			edgeLefts[8] = 1;
			edgeRghts[8] = 2;
			
			edgeOrigs[9] = 4; 
			edgeDests[9] = 0; 
			edgeLefts[9] = 2;
			edgeRghts[9]= 3;
			
			edgeOrigs[10] = 0; 
			edgeDests[10] = 5; 
			edgeLefts[10] = 4;
			edgeRghts[10] = 3;
			
			edgeOrigs[11] = 0; 
			edgeDests[11] = 6; 			
			edgeLefts[11] = 5;
			edgeRghts[11] = 4;

			for ( int k = 0; k < 12; ++k)
			{
				edges.insertEdge( edgeOrigs[k], edgeDests[k], edgeLefts[k], edgeRghts[k], particles);
			}
			assert( edges.N() == 12 );
			
			triFaceVerts[0][0] = 0;
			triFaceVerts[0][1] = 1; 
			triFaceVerts[0][2] = 2;
			
			triFaceVerts[1][0] = 2;
			triFaceVerts[1][1] = 3;
			triFaceVerts[1][2] = 0; 
			
			triFaceVerts[2][0] = 4;
			triFaceVerts[2][1] = 0; 
			triFaceVerts[2][2] = 3;
			
			triFaceVerts[3][0] = 0;
			triFaceVerts[3][1] = 4;
			triFaceVerts[3][2] = 5; 
			
			triFaceVerts[4][0] = 5;
			triFaceVerts[4][1] = 6;
			triFaceVerts[4][2] = 0;
			
			triFaceVerts[5][0] = 1;
			triFaceVerts[5][1] = 0;
			triFaceVerts[5][2] = 6;
			
			triFaceEdges[0][0] = 0;
			triFaceEdges[0][1] = 1; 
			triFaceEdges[0][2] = 7;
			
			triFaceEdges[1][0] = 2;
			triFaceEdges[1][1] = 8;
			triFaceEdges[1][2] = 7; 
			
			triFaceEdges[2][0] = 9;
			triFaceEdges[2][1] = 8; 
			triFaceEdges[2][2] = 3;
			
			triFaceEdges[3][0] = 9;
			triFaceEdges[3][1] = 4;
			triFaceEdges[3][2] = 10; 
			
			triFaceEdges[4][0] = 5;
			triFaceEdges[4][1] = 11;
			triFaceEdges[4][2] = 10;
			
			triFaceEdges[5][0] = 0;
			triFaceEdges[5][1] = 11;
			triFaceEdges[5][2] = 6;
			
			for ( int k = 0; k < 6; ++k )
			{
				faces.insertFace( 7 + k, triFaceVerts[k], triFaceEdges[k] );
				faces.setArea( k, particles, edges);
			}
			faces.setNActive(6);
			setSurfaceAreaFromParticles();
			estimateMeshSize();
			
			assert( faces.N() == 6 );
			assert( faces.nActive() == 6 );
			break;
		}
		case PolyMesh2d::quadRectSeed :
		{
			std::vector< xyzVector > coords(13);
			std::vector<int> edgeOrigs(12,-1) ;
			std::vector<int> edgeDests(12,-1) ;
			std::vector<int> edgeLefts(12,-1) ;
			std::vector<int> edgeRghts(12,-1) ;
			std::vector< std::vector<int> > quadFaceVerts(4, std::vector<int>(4, -1) );
			std::vector< std::vector<int> > quadFaceEdges(4, std::vector<int>(4, -1) );
			std::vector< std::vector< std::pair<double,int> > > incidentEdges(9);
			
			coords[0] = xyzVector(-1.0, 1.0 );
			coords[1] = xyzVector(-1.0, 0.0 );
			coords[2] = xyzVector(-1.0,-1.0 );
			coords[3] = xyzVector( 0.0,-1.0 );
			coords[4] = xyzVector( 1.0,-1.0 );
			coords[5] = xyzVector( 1.0, 0.0 );
			coords[6] = xyzVector( 1.0, 1.0 );
			coords[7] = xyzVector( 0.0, 1.0 );
			coords[8] = xyzVector( 0.0, 0.0 );
			coords[9] = xyzVector(-0.5, 0.5 );
			coords[10]= xyzVector(-0.5,-0.5 );
			coords[11]= xyzVector( 0.5,-0.5 );
			coords[12]= xyzVector( 0.5, 0.5 );
			for ( int k = 0; k < 13; ++k )
			{
				coords[k].scale(ampFactor);
				particles.insertParticle( coords[k], coords[k] );
			}
			assert( particles.N() == 13 );
			
			incidentEdges[0].push_back(std::pair<double,int>(-physConsts->Pi()/2.0, 0) );
			incidentEdges[0].push_back(std::pair<double,int>(0.0, 7 ) );

			incidentEdges[1].push_back(std::pair<double,int>( -physConsts->Pi()/2.0, 1 ) );
			incidentEdges[1].push_back(std::pair<double,int>( 0.0, 8 ) );
			incidentEdges[1].push_back(std::pair<double,int>( physConsts->Pi()/2.0, 0 ) );
			
			incidentEdges[2].push_back(std::pair<double,int>(0.0, 2 ) );
			incidentEdges[2].push_back(std::pair<double,int>(physConsts->Pi()/2.0, 1 ) );
			
			incidentEdges[3].push_back(std::pair<double,int>(-physConsts->Pi()/2.0, 2) );
			incidentEdges[3].push_back(std::pair<double,int>( 0.0, 3) );
			incidentEdges[3].push_back(std::pair<double,int>(physConsts->Pi()/2.0, 10 ) );
			
			incidentEdges[4].push_back(std::pair<double,int>(-physConsts->Pi(), 3 ));
			incidentEdges[4].push_back(std::pair<double,int>( physConsts->Pi()/2.0, 4) );
			
			incidentEdges[5].push_back(std::pair<double,int>(-physConsts->Pi(), 9 ));
			incidentEdges[5].push_back(std::pair<double,int>(-physConsts->Pi()/2.0, 4));
			incidentEdges[5].push_back(std::pair<double,int>(physConsts->Pi()/2.0, 5 ));
			
			incidentEdges[6].push_back(std::pair<double,int>(-physConsts->Pi()/2.0, 5));
			incidentEdges[6].push_back(std::pair<double,int>(-physConsts->Pi(), 6 ));
			
			incidentEdges[7].push_back(std::pair<double,int>(-physConsts->Pi(), 7));
			incidentEdges[7].push_back(std::pair<double,int>(-physConsts->Pi()/2.0, 11 ));
			incidentEdges[7].push_back(std::pair<double,int>( 0.0, 6 ) );
			
			incidentEdges[8].push_back(std::pair<double,int>(-physConsts->Pi(),8));
			incidentEdges[8].push_back(std::pair<double,int>(-physConsts->Pi()/2.0,10));
			incidentEdges[8].push_back(std::pair<double,int>(0.0, 9 ));
			incidentEdges[8].push_back(std::pair<double,int>(physConsts->Pi()/2.0,11));
			
			for ( int k = 0; k < 9; ++k )
			{
				particles.incidentEdges[k] = incidentEdges[k];
			}
			
			for ( int k = 0; k < 8; ++k )
			{
				edgeOrigs[k] = k;
				edgeDests[k] = (k+1)%8;
			}
			edgeOrigs[8] = 1;
			edgeDests[8] = 8;
			edgeOrigs[9] = 8;
			edgeDests[9] = 5;
			edgeOrigs[10]= 3;
			edgeDests[10]= 8;
			edgeOrigs[11]= 8;
			edgeDests[11]= 7;
			edgeLefts[0] = 0;
			edgeLefts[1] = 1;
			edgeLefts[2] = 1;
			edgeLefts[3] = 2;
			edgeLefts[4] = 2;
			edgeLefts[5] = 3;
			edgeLefts[6] = 3;
			edgeLefts[7] = 0;
			edgeLefts[8] = 0;
			edgeLefts[9] = 3;
			edgeLefts[10]= 1;
			edgeLefts[11]= 0;
			edgeRghts[8] = 1;
			edgeRghts[9] = 2;
			edgeRghts[10]= 2;
			edgeRghts[11]= 3;
			
			for ( int k = 0; k < 12; ++k)
			{
				edges.insertEdge( edgeOrigs[k], edgeDests[k], edgeLefts[k], edgeRghts[k], particles );
			}
			assert( edges.N() == 12 );
			
			quadFaceVerts[0][0] = 0;
			quadFaceVerts[0][1] = 1;
			quadFaceVerts[0][2] = 8;
			quadFaceVerts[0][3] = 7;
			
			quadFaceVerts[1][0] = 1;
			quadFaceVerts[1][1] = 2;
			quadFaceVerts[1][2] = 3;
			quadFaceVerts[1][3] = 8;
			
			quadFaceVerts[2][0] = 8;
			quadFaceVerts[2][1] = 3;
			quadFaceVerts[2][2] = 4;
			quadFaceVerts[2][3] = 5;
			
			quadFaceVerts[3][0] = 7;
			quadFaceVerts[3][1] = 8;
			quadFaceVerts[3][2] = 5;
			quadFaceVerts[3][3] = 6;
			
			quadFaceEdges[0][0] = 0;
			quadFaceEdges[0][1] = 8;
			quadFaceEdges[0][2] = 11;
			quadFaceEdges[0][3] = 7;
			
			quadFaceEdges[1][0] = 1;
			quadFaceEdges[1][1] = 2;
			quadFaceEdges[1][2] = 10;
			quadFaceEdges[1][3] = 8;
			
			quadFaceEdges[2][0] = 10;
			quadFaceEdges[2][1] = 3;
			quadFaceEdges[2][2] = 4;
			quadFaceEdges[2][3] = 9;
			
			quadFaceEdges[3][0] = 11;
			quadFaceEdges[3][1] = 9;
			quadFaceEdges[3][2] = 5;
			quadFaceEdges[3][3] = 6;
	
			for ( int k = 0; k < 4; ++k )
			{
				faces.insertFace( 9 + k, quadFaceVerts[k], quadFaceEdges[k] );
				faces.setArea(k, particles, edges);
			}
			faces.setNActive(4);
			
			setSurfaceAreaFromParticles();
			estimateMeshSize();
			assert(faces.N() == 4 );
			assert(faces.nActive() == 4);
			
			break;
		}
		case PolyMesh2d::polarDiscSeed :
		{
			
			break;
		}
		case PolyMesh2d::quadRectPeriodic :
		{
			break;
		}
		case PolyMesh2d::icosTriSphereSeed :
		{
			std::vector< xyzVector > coords(31);
			std::vector<int> edgeOrigs(30);
			std::vector<int> edgeDests(30);
			std::vector<int> edgeLefts(30);
			std::vector<int> edgeRights(30);
			std::vector< std::vector<int> > triFaceVerts(20, std::vector<int>(3));
			std::vector< std::vector<int> > triFaceEdges(20, std::vector<int>(3));
			std::vector< std::vector< std::pair<double, int> > > incidentEdges(12);
			
			coords[0] = xyzVector(0.0, 0.0, 1.0);
			coords[1] = xyzVector(0.723606797749978969640917366873, 0.525731112119133606025669084848, 0.447213595499957939281834733746);
			coords[2] = xyzVector( -0.276393202250021030359082633126, 0.850650808352039932181540497063, 0.447213595499957939281834733746);
			coords[3] = xyzVector(-0.894427190999915878563669467492, 0.0, 0.447213595499957939281834733746 );
			coords[4] = xyzVector(-0.276393202250021030359082633127, -0.850650808352039932181540497063, 0.447213595499957939281834733746 );
			coords[5] = xyzVector(0.723606797749978969640917366873, -0.525731112119133606025669084848,  0.447213595499957939281834733746 );
			coords[6] = xyzVector(0.894427190999915878563669467492,  0.0, -0.447213595499957939281834733746 );
			coords[7] = xyzVector(0.276393202250021030359082633127,  0.850650808352039932181540497063, -0.447213595499957939281834733746 );
			coords[8] = xyzVector(-0.723606797749978969640917366873, 0.525731112119133606025669084848, -0.447213595499957939281834733746 );
			coords[9] = xyzVector(-0.723606797749978969640917366873, -0.525731112119133606025669084848, -0.447213595499957939281834733746 );
			coords[10] = xyzVector(0.276393202250021030359082633127, -0.850650808352039932181540497063, -0.447213595499957939281834733746 );
			coords[11] = xyzVector( 0.0, 0.0, -1.0 );
			
			const double angleInc = 2.0 * physConsts->Pi() / 5.0;
			
			incidentEdges[0].push_back( std::pair<double,int>(0.0, 0 ));
			incidentEdges[0].push_back( std::pair<double,int>(angleInc, 2 ));
			incidentEdges[0].push_back( std::pair<double,int>(2.0 * angleInc,  4 ));
			incidentEdges[0].push_back( std::pair<double,int>(3.0 * angleInc, 6 ));
			incidentEdges[0].push_back( std::pair<double,int>(4.0 * angleInc, 8 ));
			
			incidentEdges[1].push_back( std::pair<double,int>(0.0, 1 ));
			incidentEdges[1].push_back( std::pair<double,int>(angleInc, 0 ));
			incidentEdges[1].push_back( std::pair<double,int>(2.0 * angleInc,  9 ));
			incidentEdges[1].push_back( std::pair<double,int>(3.0 * angleInc, 10 ));
			incidentEdges[1].push_back( std::pair<double,int>(4.0 * angleInc, 12 ));
			
			incidentEdges[2].push_back( std::pair<double,int>(0.0, 3 ));
			incidentEdges[2].push_back( std::pair<double,int>(angleInc, 2 ));
			incidentEdges[2].push_back( std::pair<double,int>(2.0 * angleInc,  1 ));
			incidentEdges[2].push_back( std::pair<double,int>(3.0 * angleInc, 13 ));
			incidentEdges[2].push_back( std::pair<double,int>(4.0 * angleInc, 15 ));
			
			incidentEdges[3].push_back( std::pair<double,int>(0.0, 5 ));
			incidentEdges[3].push_back( std::pair<double,int>(angleInc, 4 ));
			incidentEdges[3].push_back( std::pair<double,int>(2.0 * angleInc,  3 ));
			incidentEdges[3].push_back( std::pair<double,int>(3.0 * angleInc, 16 ));
			incidentEdges[3].push_back( std::pair<double,int>(4.0 * angleInc, 18 ));
			
			incidentEdges[4].push_back( std::pair<double,int>(0.0, 7 ));
			incidentEdges[4].push_back( std::pair<double,int>(angleInc, 6 ));
			incidentEdges[4].push_back( std::pair<double,int>(2.0 * angleInc,  5 ));
			incidentEdges[4].push_back( std::pair<double,int>(3.0 * angleInc, 19 ));
			incidentEdges[4].push_back( std::pair<double,int>(4.0 * angleInc, 21 ));
			
			incidentEdges[5].push_back( std::pair<double,int>(0.0, 9 ));
			incidentEdges[5].push_back( std::pair<double,int>(angleInc, 8 ));
			incidentEdges[5].push_back( std::pair<double,int>(2.0 * angleInc,  7 ));
			incidentEdges[5].push_back( std::pair<double,int>(3.0 * angleInc, 22 ));
			incidentEdges[5].push_back( std::pair<double,int>(4.0 * angleInc, 24 ));
			
			incidentEdges[6].push_back( std::pair<double,int>(0.0, 11 ));
			incidentEdges[6].push_back( std::pair<double,int>(angleInc, 10 ));
			incidentEdges[6].push_back( std::pair<double,int>(2.0 * angleInc,  24 ));
			incidentEdges[6].push_back( std::pair<double,int>(3.0 * angleInc, 23 ));
			incidentEdges[6].push_back( std::pair<double,int>(4.0 * angleInc, 25 ));
			
			incidentEdges[7].push_back( std::pair<double,int>(0.0, 14 ));
			incidentEdges[7].push_back( std::pair<double,int>(angleInc, 13 ));
			incidentEdges[7].push_back( std::pair<double,int>(2.0 * angleInc,  12 ));
			incidentEdges[7].push_back( std::pair<double,int>(3.0 * angleInc, 11 ));
			incidentEdges[7].push_back( std::pair<double,int>(4.0 * angleInc, 26 ));
			
			incidentEdges[8].push_back( std::pair<double,int>(0.0, 17 ));
			incidentEdges[8].push_back( std::pair<double,int>(angleInc, 16 ));
			incidentEdges[8].push_back( std::pair<double,int>(2.0 * angleInc, 15 ));
			incidentEdges[8].push_back( std::pair<double,int>(3.0 * angleInc, 14 ));
			incidentEdges[8].push_back( std::pair<double,int>(4.0 * angleInc, 27 ));
			
			incidentEdges[9].push_back( std::pair<double,int>(0.0, 20 ));
			incidentEdges[9].push_back( std::pair<double,int>(angleInc, 19 ));
			incidentEdges[9].push_back( std::pair<double,int>(2.0 * angleInc, 18 ));
			incidentEdges[9].push_back( std::pair<double,int>(3.0 * angleInc, 17 ));
			incidentEdges[9].push_back( std::pair<double,int>(4.0 * angleInc, 28 ));
			
			incidentEdges[10].push_back( std::pair<double,int>(0.0, 23 ));
			incidentEdges[10].push_back( std::pair<double,int>(angleInc, 22 ));
			incidentEdges[10].push_back( std::pair<double,int>(2.0 * angleInc, 21 ));
			incidentEdges[10].push_back( std::pair<double,int>(3.0 * angleInc, 20 ));
			incidentEdges[10].push_back( std::pair<double,int>(4.0 * angleInc, 29 ));
			
			incidentEdges[11].push_back( std::pair<double,int>(0.0, 25 ));
			incidentEdges[11].push_back( std::pair<double,int>(angleInc, 29 ));
			incidentEdges[11].push_back( std::pair<double,int>(2.0 * angleInc, 28 ));
			incidentEdges[11].push_back( std::pair<double,int>(3.0 * angleInc, 27 ));
			incidentEdges[11].push_back( std::pair<double,int>(4.0 * angleInc, 26 ));
			
			for ( int i = 0; i < 12; ++i)
			{
				coords[i].scale(ampFactor);
				particles.insertParticle( coords[i], coords[i] );
				particles.incidentEdges[i] = incidentEdges[i];
			}
			
			edgeOrigs[0] = 0;
			edgeDests[0] = 1;
			edgeLefts[0] = 0;
			edgeRights[0] = 4;
			
			edgeOrigs[1] = 1;
			edgeDests[1] = 2;
			edgeLefts[1] = 0;
			edgeRights[1] = 6;
			
			edgeOrigs[2] = 2; 
			edgeDests[2] = 0;
			edgeLefts[2] = 0; 
			edgeRights[2] = 1;

			edgeOrigs[3] = 2; 
			edgeDests[3] = 3;
			edgeLefts[3] = 1;
			edgeRights[3] = 8;

			edgeOrigs[4] = 0;
			edgeDests[4] = 3;
			edgeLefts[4] = 2;
			edgeRights[4] = 1;

			edgeOrigs[5] = 3;
			edgeDests[5] = 4;
			edgeLefts[5] = 2;
			edgeRights[5] = 10;

			edgeOrigs[6] = 4;
			edgeDests[6] = 0;
			edgeLefts[6] = 2;
			edgeRights[6] = 3;

			edgeOrigs[7] = 4;
			edgeDests[7] = 5;
			edgeLefts[7] = 3;
			edgeRights[7] = 12;

			edgeOrigs[8] = 5;
			edgeDests[8] = 0;
			edgeLefts[8] = 3;
			edgeRights[8] = 4;

			edgeOrigs[9] = 5;
			edgeDests[9] = 1;
			edgeLefts[9] = 4;
			edgeRights[9] = 14;

			edgeOrigs[10] = 1;
			edgeDests[10] = 6;
			edgeLefts[10] = 5;
			edgeRights[10] = 14;

			edgeOrigs[11] = 6;
			edgeDests[11] = 7;
			edgeLefts[11] = 5; 
			edgeRights[11] = 15; 

			edgeOrigs[12] = 7;
			edgeDests[12] = 1;
			edgeLefts[12] = 5;
			edgeRights[12] = 6;

			edgeOrigs[13] = 7;
			edgeDests[13] = 2;
			edgeLefts[13] = 6;
			edgeRights[13] = 7;

			edgeOrigs[14] = 7;
			edgeDests[14] = 8;
			edgeLefts[14] = 7;
			edgeRights[14] = 16;

			edgeOrigs[15] = 8;
			edgeDests[15] = 2;
			edgeLefts[15] = 7;
			edgeRights[15] = 8;

			edgeOrigs[16] = 8;
			edgeDests[16] = 3;
			edgeLefts[16] = 8;
			edgeRights[16] = 9;

			edgeOrigs[17] = 8;
			edgeDests[17] = 9;
			edgeLefts[17] = 9;
			edgeRights[17] = 17;

			edgeOrigs[18] = 3;
			edgeDests[18] = 9;
			edgeLefts[18] = 10;
			edgeRights[18] = 9;

			edgeOrigs[19] = 9;
			edgeDests[19] = 4;
			edgeLefts[19] = 10;
			edgeRights[19] = 11; 

			edgeOrigs[20] = 9;
			edgeDests[20] = 10;
			edgeLefts[20] = 11;
			edgeRights[20] = 18;

			edgeOrigs[21] = 10;
			edgeDests[21] = 4;
			edgeLefts[21] = 11;
			edgeRights[21] = 12;

			edgeOrigs[22] = 10;
			edgeDests[22] = 5;
			edgeLefts[22] = 12;
			edgeRights[22] = 13;

			edgeOrigs[23] = 10;
			edgeDests[23] = 6;
			edgeLefts[23] = 13;
			edgeRights[23] = 19;

			edgeOrigs[24] = 6;
			edgeDests[24] = 5;
			edgeLefts[24] = 13;
			edgeRights[24] = 14;

			edgeOrigs[25] = 11;
			edgeDests[25] = 6;
			edgeLefts[25] = 19;
			edgeRights[25] = 15;

			edgeOrigs[26] = 11;
			edgeDests[26] = 7;
			edgeLefts[26] = 15;
			edgeRights[26] = 16;

			edgeOrigs[27] = 8;
			edgeDests[27] = 11;
			edgeLefts[27] = 17;
			edgeRights[27] = 16; 

			edgeOrigs[28] = 11;
			edgeDests[28] = 9;
			edgeLefts[28] = 17;
			edgeRights[28] = 18;

			edgeOrigs[29] = 10;
			edgeDests[29] = 11;
			edgeLefts[29] = 19;
			edgeRights[29] = 18;
			
			for ( int i = 0; i < 30; ++i)
			{
				edges.insertEdge( edgeOrigs[i], edgeDests[i], edgeLefts[i], edgeRights[i], particles );
			}
			assert( edges.N() == 30);		
				
			triFaceVerts[0][0] = 0;
			triFaceVerts[0][1] = 1;
			triFaceVerts[0][2] = 2;
			
			triFaceVerts[1][0] = 0;
			triFaceVerts[1][1] = 2;
			triFaceVerts[1][2] = 3;
			
			triFaceVerts[2][0] = 0;
			triFaceVerts[2][1] = 3;
			triFaceVerts[2][2] = 4;
			
			triFaceVerts[3][0] = 0;
			triFaceVerts[3][1] = 4;
			triFaceVerts[3][2] = 5;
			
			triFaceVerts[4][0] = 0;
			triFaceVerts[4][1] = 5;
			triFaceVerts[4][2] = 1;
			
			triFaceVerts[5][0] = 1;
			triFaceVerts[5][1] = 6;
			triFaceVerts[5][2] = 7;
			
			triFaceVerts[6][0] = 7;
			triFaceVerts[6][1] = 2;
			triFaceVerts[6][2] = 1;
			
			triFaceVerts[7][0] = 2;
			triFaceVerts[7][1] = 7;
			triFaceVerts[7][2] = 8;
			
			triFaceVerts[8][0] = 8;
			triFaceVerts[8][1] = 3;
			triFaceVerts[8][2] = 2;
			
			triFaceVerts[9][0] = 3;
			triFaceVerts[9][1] = 8;
			triFaceVerts[9][2] = 9;
			
			triFaceVerts[10][0] = 9;
			triFaceVerts[10][1] = 4;
			triFaceVerts[10][2] = 3;
			
			triFaceVerts[11][0] = 4;
			triFaceVerts[11][1] = 9;
			triFaceVerts[11][2] = 10;
			
			triFaceVerts[12][0] = 10;
			triFaceVerts[12][1] = 5;
			triFaceVerts[12][2] = 4;
			
			triFaceVerts[13][0] = 5;
			triFaceVerts[13][1] = 10;
			triFaceVerts[13][2] = 6;	
			
			triFaceVerts[14][0] = 6;
			triFaceVerts[14][1] = 1;
			triFaceVerts[14][2] = 5;
			
			triFaceVerts[15][0] = 11;
			triFaceVerts[15][1] = 7;
			triFaceVerts[15][2] = 6;
			
			triFaceVerts[16][0] = 11;
			triFaceVerts[16][1] = 8;
			triFaceVerts[16][2] = 7;
			
			triFaceVerts[17][0] = 11;
			triFaceVerts[17][1] = 9;
			triFaceVerts[17][2] = 8;
			
			triFaceVerts[18][0] = 11;
			triFaceVerts[18][1] = 10;
			triFaceVerts[18][2] = 9;
			
			triFaceVerts[19][0] = 11;
			triFaceVerts[19][1] = 6;
			triFaceVerts[19][2] = 10;
			
			triFaceEdges[0][0] = 0;
			triFaceEdges[0][1] = 1;
			triFaceEdges[0][2] = 2;

			triFaceEdges[1][0] = 2; 
			triFaceEdges[1][1] = 3;
			triFaceEdges[1][2] = 4;

			triFaceEdges[2][0] = 4;
			triFaceEdges[2][1] = 5;
			triFaceEdges[2][2] = 6;

			triFaceEdges[3][0] = 6;
			triFaceEdges[3][1] = 7;
			triFaceEdges[3][2] = 8;

			triFaceEdges[4][0] = 8;
			triFaceEdges[4][1] = 9;
			triFaceEdges[4][2] = 0;

			triFaceEdges[5][0] = 10;
			triFaceEdges[5][1] = 11;
			triFaceEdges[5][2] = 12;

			triFaceEdges[6][0] = 13;
			triFaceEdges[6][1] = 1;
			triFaceEdges[6][2] = 12;

			triFaceEdges[7][0] = 13;
			triFaceEdges[7][1] = 14;
			triFaceEdges[7][2] = 15;

			triFaceEdges[8][0] = 16;
			triFaceEdges[8][1] = 3;
			triFaceEdges[8][2] = 15;

			triFaceEdges[9][0] = 16;
			triFaceEdges[9][1] = 17;
			triFaceEdges[9][2] = 18;

			triFaceEdges[10][0] = 19;
			triFaceEdges[10][1] = 5;
			triFaceEdges[10][2] = 18;

			triFaceEdges[11][0] = 19;
			triFaceEdges[11][1] = 20;
			triFaceEdges[11][2] = 21;

			triFaceEdges[12][0] = 22;
			triFaceEdges[12][1] = 7;
			triFaceEdges[12][2] = 21;

			triFaceEdges[13][0] = 22;
			triFaceEdges[13][1] = 23;
			triFaceEdges[13][2] = 24;

			triFaceEdges[14][0] = 10;
			triFaceEdges[14][1] = 9;
			triFaceEdges[14][2] = 24;

			triFaceEdges[15][0] = 26;
			triFaceEdges[15][1] = 11;
			triFaceEdges[15][2] = 25;

			triFaceEdges[16][0] = 27;
			triFaceEdges[16][1] = 14;
			triFaceEdges[16][2] = 26;

			triFaceEdges[17][0] = 28;
			triFaceEdges[17][1] = 17;
			triFaceEdges[17][2] = 27;

			triFaceEdges[18][0] = 29;
			triFaceEdges[18][1] = 20;
			triFaceEdges[18][2] = 28;

			triFaceEdges[19][0] = 25;
			triFaceEdges[19][1] = 23;
			triFaceEdges[19][2] = 29;
			
			for ( int i = 0; i < 20; ++i)
			{
				std::vector< xyzVector> faceVertCoords(3);
				xyzVector faceCenter;
				for ( int j = 0; j < 3; ++j )
				{
					faceVertCoords[j] = coords[ triFaceVerts[i][j] ];
					faceCenter = sphereCentroid( faceVertCoords );
				}
				particles.insertParticle( faceCenter, faceCenter );
			}
			
			for ( int i = 0; i < 20; ++i )
			{
				faces.insertFace( 12 +i, triFaceVerts[i], triFaceEdges[i] );
				faces.setArea( i, particles, edges );
			}
			
			faces.setNActive(20);
			setSurfaceAreaFromParticles();
			
			assert(faces.N() == 20);



			break;
		}
		case PolyMesh2d::cubedSphereSeed :
		{
			const double a = 1.0 / std::sqrt(3.0);
			std::vector< xyzVector > coords(14);
			std::vector<int> edgeOrigs(12);
			std::vector<int> edgeDests(12);
			std::vector<int> edgeLefts(12);
			std::vector<int> edgeRights(12);
			std::vector< std::vector<int> > quadFaceVerts(6, std::vector<int>(4) );
			std::vector< std::vector<int> > quadFaceEdges(6, std::vector<int>(4) );
			std::vector< std::vector< std::pair<double, int> > > incidentEdges(8);
			
			coords[0] = xyzVector( a, -a,  a );
			coords[1] = xyzVector( a, -a, -a );
			coords[2] = xyzVector( a,  a, -a );
			coords[3] = xyzVector( a,  a,  a );
			coords[4] = xyzVector(-a,  a, -a );
			coords[5] = xyzVector(-a,  a,  a );
			coords[6] = xyzVector(-a, -a, -a );
			coords[7] = xyzVector(-a, -a,  a );
			coords[8] = xyzVector( 1.0, 0.0, 0.0 );
			coords[9] = xyzVector( 0.0, 1.0, 0.0 );
			coords[10] = xyzVector( -1.0, 0.0, 0.0 );
			coords[11] = xyzVector( 0.0, -1.0, 0.0 );
			coords[12] = xyzVector( 0.0, 0.0, 1.0 );
			coords[13] = xyzVector( 0.0, 0.0, -1.0 );
	
			for (int i = 0; i < 14; ++i )
			{
				coords[i].scale(ampFactor);
				particles.insertParticle(coords[i], coords[i]);
			}
			
			const double PI = physConsts->Pi();
			
			incidentEdges[0].push_back( std::pair<double,int>( -2.0*PI/3.0, 0 ) );
			incidentEdges[0].push_back( std::pair<double,int>( 0.0, 3 ) );
			incidentEdges[0].push_back( std::pair<double,int>( 2.0*PI/3.0, 11 ) );
			
			incidentEdges[1].push_back( std::pair<double,int>( -2.0*PI/3.0, 10 ) );
			incidentEdges[1].push_back( std::pair<double,int>( 0.0, 1 ) );
			incidentEdges[1].push_back( std::pair<double,int>( 2.0*PI/3.0, 0 ) );
			
			incidentEdges[2].push_back( std::pair<double,int>( -2.0*PI/3.0, 4 ) );
			incidentEdges[2].push_back( std::pair<double,int>( 0.0, 2 ) );
			incidentEdges[2].push_back( std::pair<double,int>( 2.0*PI/3.0, 1 ) );
			
			incidentEdges[3].push_back( std::pair<double,int>( -2.0*PI/3.0, 2 ) );
			incidentEdges[3].push_back( std::pair<double,int>( 0.0, 6 ) );
			incidentEdges[3].push_back( std::pair<double,int>( 2.0*PI/3.0, 3 ) );
			
			incidentEdges[4].push_back( std::pair<double,int>( -2.0*PI/3.0, 7 ) );
			incidentEdges[4].push_back( std::pair<double,int>( 0.0, 5 ) );
			incidentEdges[4].push_back( std::pair<double,int>( 2.0*PI/3.0, 4 ) );
			
			incidentEdges[5].push_back( std::pair<double,int>( -2.0*PI/3.0, 5 ) );
			incidentEdges[5].push_back( std::pair<double,int>( 0.0, 9 ) );
			incidentEdges[5].push_back( std::pair<double,int>( 2.0*PI/3.0, 6 ) );
			
			incidentEdges[6].push_back( std::pair<double,int>( -2.0*PI/3.0, 10 ) );
			incidentEdges[6].push_back( std::pair<double,int>( 0.0, 8 ) );
			incidentEdges[6].push_back( std::pair<double,int>( 2.0*PI/3.0, 7 ) );
			
			incidentEdges[7].push_back( std::pair<double,int>( -2.0*PI/3.0, 8 ) );
			incidentEdges[7].push_back( std::pair<double,int>( 0.0, 11 ) );
			incidentEdges[7].push_back( std::pair<double,int>( 2.0*PI/3.0, 9 ) );
			
			for ( int i = 0; i < 8; ++i)
			{
				particles.incidentEdges[i] = incidentEdges[i];
			}
			
			edgeOrigs[0] = 0;
			edgeDests[0] = 1;
			edgeLefts[0] = 0;
			edgeRights[0] = 3;
			
			edgeOrigs[1] = 1;
			edgeDests[1] = 2;
			edgeLefts[1] = 0;
			edgeRights[1] = 5;
			
			edgeOrigs[2] = 2;
			edgeDests[2] = 3;
			edgeLefts[2] = 0;
			edgeRights[2] = 1;
			
			edgeOrigs[3] = 3;
			edgeDests[3] = 0;
			edgeLefts[3] = 0;
			edgeRights[3] = 4;
			
			edgeOrigs[4] = 2;
			edgeDests[4] = 4;
			edgeLefts[4] = 1;
			edgeRights[4] = 5;
			
			edgeOrigs[5] = 4;
			edgeDests[5] = 5;
			edgeLefts[5] = 1;
			edgeRights[5] = 2;
			
			edgeOrigs[6] = 5;
			edgeDests[6] = 3;
			edgeLefts[6] = 1;
			edgeRights[6] = 4;
			
			edgeOrigs[7] = 4;
			edgeDests[7] = 6;
			edgeLefts[7] = 2;
			edgeRights[7] = 5;
			
			edgeOrigs[8] = 6;
			edgeDests[8] = 7;
			edgeLefts[8] = 2;
			edgeRights[8] = 3;
			
			edgeOrigs[9] = 7;
			edgeDests[9] = 5;
			edgeLefts[9] = 2;
			edgeRights[9] = 4;
			
			edgeOrigs[10] = 6;
			edgeDests[10] = 1;
			edgeLefts[10] = 3;
			edgeRights[10] = 5;
			
			edgeOrigs[11] = 0;
			edgeDests[11] = 7;
			edgeLefts[11] = 3;
			edgeRights[11] = 4;
			
			for (int i = 0; i < 12; ++i )
			{
				edges.insertEdge( edgeOrigs[i], edgeDests[i], edgeLefts[i], edgeRights[i], particles );
			}	
			assert( edges.N() == 12);
			
			quadFaceVerts[0][0] = 0;
			quadFaceVerts[0][1] = 1;
			quadFaceVerts[0][2] = 2;
			quadFaceVerts[0][3] = 3; 
			
			quadFaceVerts[1][0] = 3;
			quadFaceVerts[1][1] = 2;
			quadFaceVerts[1][2] = 4;
			quadFaceVerts[1][3] = 5;
			
			quadFaceVerts[2][0] = 5;
			quadFaceVerts[2][1] = 4;
			quadFaceVerts[2][2] = 6;
			quadFaceVerts[2][3] = 7;
			
			quadFaceVerts[3][0] = 7;
			quadFaceVerts[3][1] = 6;
			quadFaceVerts[3][2] = 1;
			quadFaceVerts[3][3] = 0;
			
			quadFaceVerts[4][0] = 7;
			quadFaceVerts[4][1] = 0;
			quadFaceVerts[4][2] = 3;
			quadFaceVerts[4][3] = 5;
			
			quadFaceVerts[5][0] = 1;
			quadFaceVerts[5][1] = 6;
			quadFaceVerts[5][2] = 4;
			quadFaceVerts[5][3] = 2;
			
			quadFaceEdges[0][0] = 0;
			quadFaceEdges[0][1] = 1;
			quadFaceEdges[0][2] = 2;
			quadFaceEdges[0][3] = 3;
			
			quadFaceEdges[1][0] = 2;
			quadFaceEdges[1][1] = 4;
			quadFaceEdges[1][2] = 5;
			quadFaceEdges[1][3] = 6;
			
			quadFaceEdges[2][0] = 5;
			quadFaceEdges[2][1] = 7;
			quadFaceEdges[2][2] = 8;
			quadFaceEdges[2][3] = 9;
			
			quadFaceEdges[3][0] = 8;
			quadFaceEdges[3][1] = 10;
			quadFaceEdges[3][2] = 0;
			quadFaceEdges[3][3] = 11;
			
			quadFaceEdges[4][0] = 11;
			quadFaceEdges[4][1] = 3;
			quadFaceEdges[4][2] = 6;
			quadFaceEdges[4][3] = 9;
			
			quadFaceEdges[5][0] = 10;
			quadFaceEdges[5][1] = 7;
			quadFaceEdges[5][2] = 4;
			quadFaceEdges[5][3] = 1;
			
			for ( int i = 0; i < 6; ++i )
			{
				faces.insertFace( 8 + i, quadFaceVerts[i], quadFaceEdges[i] );
				faces.setArea(i, particles, edges);
			}
			faces.setNActive(6);
			setSurfaceAreaFromParticles();
			
			assert( faces.N() == 6 );	
			break;
		}
		default :
		{
			statusMsg = OutputMessage( "invalid mesh seed.", OutputMessage::errorPriority,
													  "PolyMesh2d::initializeFromSeed");
			log->logMessage(statusMsg);
			break;
		}
	}
};

int PolyMesh2d::locatePointWalkSearch( const xyzVector& queryPt, const int startIndex ) const
{
	if ( faces.hasChildren(startIndex) )
	{
		OutputMessage statusMsg("expected a low-level face as input", OutputMessage::errorPriority, 
								"PolyMesh2d::locatePointWalkSearch" );
		log->logMessage(statusMsg);
		return -1;
	}
	xyzVector cntd = faceCentroid(startIndex);
	double currentMin = distance( cntd, queryPt, gk );
	int currentFace = startIndex;
	
	std::vector<int> adjFaces = ccwAdjacentFaces( currentFace );
	for ( int i = 0; i < adjFaces.size(); ++i)
	{
		if ( adjFaces[i] >= 0 )
		{
			cntd = faceCentroid( adjFaces[i] );
			double testDist = distance( cntd, queryPt, gk );
			if ( testDist < currentMin )
			{
				currentMin = testDist;
				currentFace = adjFaces[i];
			}
		}
	}
	if ( currentFace == startIndex )
		return currentFace;
	else
		return locatePointWalkSearch( queryPt, currentFace );
};

int PolyMesh2d::nVertices( const int initNest ) const
{
	int result = 0;
	OutputMessage statusMsg;
	
	switch ( _seed )
	{
		case triHexSeed :
		{
			// int steps = std::pow(2,initNest);
			// for ( int k = 0; k <= steps; ++k)
			for (int k = std::pow(2,initNest) + 1; k <= std::pow(2,initNest+1); ++k)
			{
				result +=  k;
				std::cout << "DEBUG k " << k <<": nVerts = " << result << std::endl;
			}
			result = 2 * result + std::pow(2,initNest+1) + 1;
			std::cout << "DEBUG : nVerts = " << result << std::endl;
			break;
		}
		case quadRectSeed :
		{
			result = 3;
			for ( int i = 1; i <= initNest; ++i )
			{
				result += std::pow(2,i);
			}
			result *= result;
			break;
		}
		case polarDiscSeed :
		{
			break;
		}
		case quadRectPeriodic :
		{
			result = 3;
			for ( int i = 1; i <= initNest; ++i )
			{
				result += std::pow(2,i);
			}
			result *= result;
			break;
		}
		case icosTriSphereSeed :
		{
			result = 2 + 10 * std::pow(4,initNest);
			break;
		}
		case cubedSphereSeed :
		{
			result = 2 + 6 * std::pow(4,initNest);
			break;
		}
	}
	return result;
};

int PolyMesh2d::nFaces( const int initNest ) const
{
	OutputMessage statusMsg;
	int result = 0;
	switch ( _seed )
	{
		case triHexSeed :
		{
			result = 6 * std::pow(4, initNest );
			break;
		}
		case quadRectSeed :
		{
			result = 4 * std::pow(4, initNest );
			break;
		}
		case polarDiscSeed :
		{
			break;
		}
		case quadRectPeriodic :
		{
			break;
		}
		case icosTriSphereSeed :
		{
			result = 20 * std::pow(4,initNest);
			break;
		}
		case cubedSphereSeed :
		{
			result = 6 * std::pow(4,initNest);
			break;
		}
	}
	return result;
};

int PolyMesh2d::nEdges( const int nVert, const int nFaces ) const
{
	int result = 0;
	switch ( gk )
	{
		case ( xyzVector::EuclideanGeometry ) :
		{
			result = nFaces + nVert - 1;
			break;
		}
		case ( xyzVector::SphericalGeometry ) :
		{
			result = nFaces + nVert - 2;
			break;
		}
	}
	return result;
};

void PolyMesh2d::outputToMatlab( const std::string filename ) const
{
	std::ofstream file( filename.c_str() );
	if ( !file )
	{
		OutputMessage statusMsg("cannot open output file.", OutputMessage::errorPriority, "PolyMesh2d::outputToMatlab");
		log->logMessage(statusMsg);
		return;
	}
	particles.writeVariablesToMatlab( file );
	edges.writeVariablesToMatlab( file );
	faces.writeVariablesToMatlab( file );
};

void PolyMesh2d::sortEdgesAroundVertices()
{
	std::vector< std::pair<double,int> > nullVec;
	for ( int i = 0; i < particles.N(); ++i)
	{
		if ( particles.area[i] > 0.0 )
			particles.incidentEdges[i] = nullVec;
		else
			particles.sortIncidentEdgesAtParticle(i);
	}
};

void PolyMesh2d::outputToVTK( const std::string filename, const std::vector<Field> fields ) const
{
	std::ofstream file(filename.c_str() );
	if ( !file )
	{
		OutputMessage statusMsg("cannot open output file.", OutputMessage::errorPriority, "PolyMesh2d::outputToVTK");
		log->logMessage(statusMsg);
		return;
	}
	particles.writePointsToVTK( file );
	faces.writePolygonsToVTK( file );
	particles.writeLagCoordToVTK( file );
	if ( !fields.empty() )
	{
		for (int i = 0; i < fields.size(); ++i)
		{
			fields[i].writePointDataToVTK( file );
		}
	}
};

void PolyMesh2d::writeActiveParticlesToCSV(const std::string filename) const{
    std::ofstream file(filename.c_str());
	if ( !file )
	{
		OutputMessage statusMsg("cannot open output file.", OutputMessage::errorPriority, 
		    "PolyMesh2d::writeActiveParticlesToCSV");
		log->logMessage(statusMsg);
		return;
	}
	file << "x,   y,   z,   boundary_id,   area\n";
	file << faces.nActive() << std::endl;
	const std::string sepStr = ",";
	for (int i = 0; i < particles.N(); ++i){
	    if ( particles.area[i] > 0.0 ){
            file << std::setprecision(std::numeric_limits<long double>::digits10 + 1) <<
            particles.x[i] << sepStr << std::setprecision(std::numeric_limits<long double>::digits10 + 1) <<
            particles.y[i] << sepStr << std::setprecision(std::numeric_limits<long double>::digits10 + 1) <<
            particles.z[i] << sepStr << "0" << sepStr << std::setprecision(std::numeric_limits<long double>::digits10 + 1) <<
            particles.area[i] << std::endl;
//             std::cout << "i = " << i << ", nactive = " << faces.nActive() << ", particles.N = " << particles.N() << std::endl;
	    }
	}
}

void PolyMesh2d::writeActiveParticlesToVTK(const std::string filename, const std::string title) const {
    std::ofstream file(filename.c_str());
	if ( !file )
	{
		OutputMessage statusMsg("cannot open output file.", OutputMessage::errorPriority, 
		    "PolyMesh2d::writeActiveParticlesToVTK");
		log->logMessage(statusMsg);
		return;
	}
	file << "# vtk DataFile Version 2.0 " << std::endl;
	file << title << std::endl;
	file << "ASCII" << std::endl;
	file << "DATASET POLYDATA" << std::endl;
	
	file << "POINTS " << faces.nActive() << " double " << std::endl;
	const std::string sepStr = "  ";
	for (int i = 0; i < particles.N(); ++i){
	    if ( particles.area[i] > 0.0 ){
            file << std::setprecision(std::numeric_limits<long double>::digits10 + 1) <<
            particles.x[i] << sepStr << std::setprecision(std::numeric_limits<long double>::digits10 + 1) <<
            particles.y[i] << sepStr << std::setprecision(std::numeric_limits<long double>::digits10 + 1) <<
            particles.z[i] << std::endl;
        }    
     }
     file << "POINT_DATA " << faces.nActive() << std::endl;
     file << "SCALARS area double 1" << std::endl;
     file << "LOOKUP_TABLE default" << std::endl;
     for (int i = 0; i < particles.N(); ++i) {
        if (particles.area[i] > 0.0) {
            file << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << particles.area[i] << std::endl;
//             std::cout << "i = " << i << ", nactive = " << faces.nActive() << ", particles.N = " << particles.N() << std::endl;
	    }
	}
}

std::vector< std::string > PolyMesh2d::getInfo() const 
{
	std::vector< std::string > info;
	std::stringstream ss;
	std::string s;
	
	ss << " nFaces = " << faces.N();
	info.push_back( ss.str() );
	ss.str(s);
	
	ss << " nActive Faces = " << faces.nActive();
	info.push_back( ss.str() );
	ss.str(s);
	
	ss << " nEdges = " << edges.N();
	info.push_back( ss.str() );
	ss.str(s);
	
	ss << " nParticles = " << particles.N();
	info.push_back( ss.str() );
	ss.str(s);
	
	ss << " mesh surf area = " << _surfaceArea;
	info.push_back( ss.str() );
	ss.str(s);

	ss << " mesh size = " << _meshSpacing;
	info.push_back( ss.str() );
	ss.str(s);
	
	ss << " PARTICLES INFO : " ;
	info.push_back( ss.str() );
	ss.str(s);
	
	std::vector< std::string > particlesInfo = particles.getInfo();
	for ( int i = 0; i < particlesInfo.size(); ++i )
		info.push_back( particlesInfo[i] );

	ss << " EDGES INFO : " ;
	info.push_back( ss.str() );
	ss.str(s);
	std::vector< std::string > edgesInfo = edges.getInfo();
	for ( int i = 0; i < edgesInfo.size(); ++i)
		info.push_back( edgesInfo[i] );

	ss << " FACES INFO : ";
	info.push_back( ss.str() );
	ss.str();
	std::vector< std::string > facesInfo = faces.getInfo();
	for ( int i = 0; i < facesInfo.size(); ++i)
		info.push_back( facesInfo[i] );

	return info;	
};

std::ostream& operator<<( std::ostream& os, const PolyMesh2d& aMesh)
{
	os << "planar mesh nFaces = " << aMesh.nFaces() << std::endl;
	os << "planar mesh nActive = " << aMesh.nActiveFaces() << std::endl;
	os << "planar mesh nEdges = " << aMesh.nEdges() << std::endl;
	os << "planar mesh nParticles = " << aMesh.nParticles() << std::endl;
	os << "planar mesh surface area = " << aMesh.surfaceArea() << std::endl;
	
	return os;	
};