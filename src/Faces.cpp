//
//  Faces.cpp
//  LPM
//
//  Created by Peter Bosler on 11/5/14.
//  Copyright (c) 2014 Peter Bosler. All rights reserved.
//

/** @file Faces.cpp
	@brief Faces class implementation
	@author Peter Bosler, Sandia National Laboratories, Multiphysics Applications
*/

#include "Faces.h"
#include <vector>
#include <assert.h>
#include "xyzVector.h"
#include "OutputMessage.h"
#include "Logger.h"
#include <sstream>
#include <string>
#include "xyzVector.h"

Faces::Faces( const int nMax, const facesKind kind, const int procRank, const int numProcs, 
			  const xyzVector::geometryKind gkin)
{
	assert( nMax > 0 );
	assert( procRank >= 0 );
	assert( numProcs >= 1 );
	
	logLevel = OutputMessage::debugPriority;
	log = Logger::Instance( logLevel, procRank, numProcs);
	
	_N = 0;
	_nMax = nMax;
	_nActive = 0;
	_kind = kind;
	gk = gkin;
	
	int nVertsPerFace;
	if ( kind == Faces::triangularFaces )
		nVertsPerFace = 3;
	else if ( kind == Faces::quadrilateralFaces )
		nVertsPerFace = 4;
	else
	{
		OutputMessage statusMsg("invalid facesKind", OutputMessage::errorPriority, "Faces constructor");
		log->logMessage(statusMsg);
	}
	
	_centerParticle = std::vector<int>(nMax, -1);
	_vertices = std::vector< std::vector<int> >(nMax, std::vector<int>(nVertsPerFace,-1) );
	_edges = std::vector< std::vector<int> >(nMax, std::vector<int>(nVertsPerFace,-1)  );
	_hasChildren = std::vector<bool>(nMax, false);
	_children = std::vector< std::vector<int> >( nMax, std::vector<int>(4,-1) );
	_parent = std::vector<int>(nMax, -1);
};


void Faces::setArea( const int index, Particles& particles, const Edges& edges ) const 
{
	double area = 0.0;

	std::vector<xyzVector> tri(3);
	tri[1] = particles.physCoord( _centerParticle[index] );
	
	const int nverts = _vertices[index].size();
	
	for ( int i = 0; i < nverts; ++i )
	{
		tri[0] = particles.physCoord( _vertices[index][i] );
		tri[2] = particles.physCoord( _vertices[index][ (i+1)%nverts ] );

		if ( gk == xyzVector::EuclideanGeometry )
			area += cartTriArea( tri[0], tri[1], tri[2] );
		else if ( gk == xyzVector::SphericalGeometry )
			area += sphereTriArea( tri[0], tri[1], tri[2] );
	}
	
	particles.area[_centerParticle[index]] = area;	
};

/** @brief Inserts a face to the end of the Faces data structure.
	
	@param cntrParticle
	@param vertices
	@param edges
*/
void Faces::insertFace( const int cntrParticle, std::vector<int> vertices, std::vector<int> edges )
{
	assert( _N+1 <= _nMax );
	assert( vertices.size() == _vertices[_N].size() );
	assert( edges.size() == _edges[_N].size() );
	
	_centerParticle[_N] = cntrParticle;
	_vertices[_N] = vertices;
	_edges[_N] = edges;	
	
	_N++;
};

/** @brief Divides a triangular parent face to create four triangular child faces.
	
	The center particle of the parent face becomes the center particle of child 4.
	
	@param index
	@param particles
	@param edges
*/
void Faces::divideTriFace( const int index, Particles& particles, Edges& edges )
{
	assert( _N+4 <= _nMax );
	assert( _edges[index].size() == 3 );

	OutputMessage statusMsg;
	std::string s;
	std::stringstream ss;
	
	std::vector< std::vector<int> > newFaceVerts(4, std::vector<int>(3, -1) );
	std::vector< std::vector<int> > newFaceEdges(4, std::vector<int>(3, -1) );
	
	//
	// connect child faces to vertices of parent
	//
	for ( int i = 0; i < 3; ++i )
		newFaceVerts[i][i] = _vertices[index][i];
	
	//
	// loop over parent edges
	//
	int k = 0;
	int parentEdge, childEdge1, childEdge2;
	for ( int i = 0; i < 3; ++i )
	{
		parentEdge = _edges[index][i];
				
		if ( edges.hasChildren( parentEdge ) ) // parent edge has already been divided
		{	
			childEdge1 = edges.child1( parentEdge );
			childEdge2 = edges.child2( parentEdge );

			// connect child faces to child edges
			if ( positiveEdge( index, parentEdge, edges ) )
			{
				// parent edge has positive orientation
				newFaceEdges[k][i] = childEdge1;
				edges.setLeftFace( childEdge1, _N + k );
				newFaceEdges[(++k)%3][i] = childEdge2;
				edges.setLeftFace( childEdge2, _N + k%3 );
			}
			else
			{
				// parent edge has negative orientation
				newFaceEdges[k][i] = childEdge2;
				edges.setRightFace( childEdge2, _N + k );
				newFaceEdges[(++k)%3][i] = childEdge1;
				edges.setRightFace( childEdge1, _N + k%3 );
			}
		}
		else  // parent edge needs to be divided
		{
			childEdge1 = edges.N();
			childEdge2 = edges.N() + 1;							
			edges.divideEdge( parentEdge, particles );
			
			if ( positiveEdge( index, parentEdge, edges ) )
			{
				newFaceEdges[k][i] = childEdge1;
				edges.setLeftFace( childEdge1, _N + k );
				newFaceEdges[(++k)%3][i] = childEdge2;
				edges.setLeftFace( childEdge2, _N + k%3 );
			}
			else
			{
				newFaceEdges[k][i] = childEdge2;
				edges.setRightFace( childEdge2, _N + k);
				newFaceEdges[(++k)%3][i] = childEdge1;
				edges.setRightFace( childEdge1, _N + k%3 );
			}			
		}
		// connect child faces to particle at midpoint of parent edge
		newFaceVerts[i][(i+1)%3] = edges.dest( childEdge1 );
		newFaceVerts[(i+1)%3][i] = edges.dest( childEdge1 );
		if ( i == 2 )
			assert(k == 3);
	}
	//
	// DEBUG : check vertex connectivity
	//
	if ( logLevel == OutputMessage::debugPriority )
	{
		for ( int i = 0; i < 3; ++i )
		{
			for ( int j = 0; j < 3; ++j )
			{
				if ( newFaceVerts[i][j] < 0 )
				{
					ss << "exterior connectivity error at parent face "<< index << " child " << i << " vertex " << j ;
					statusMsg = OutputMessage(ss.str(), OutputMessage::errorPriority,"Faces::divideTriFace");
					log->logMessage(statusMsg);
					ss.str(s);
				}
			}			
		}
	}

	//
	// new interior edges
	//
	newFaceVerts[3][0] = newFaceVerts[2][1];
	newFaceVerts[3][1] = newFaceVerts[2][0];
	newFaceVerts[3][2] = newFaceVerts[0][1];
	k = edges.N();
	
	edges.insertEdge( newFaceVerts[3][0], newFaceVerts[3][1], _N+3, _N+2, particles);
	newFaceEdges[2][0] = k;
	newFaceEdges[3][0] = k++;
	
	edges.insertEdge( newFaceVerts[3][1], newFaceVerts[3][2], _N+3, _N, particles );
	newFaceEdges[0][1] = k;
	newFaceEdges[3][1] = k++;
	
	edges.insertEdge( newFaceVerts[3][2], newFaceVerts[3][0], _N+3, _N+1, particles);
	newFaceEdges[1][2] = k;
	newFaceEdges[3][2] = k;

	//
	// DEBUG : check vertex connectivity
	//
	if ( logLevel == OutputMessage::debugPriority )
	{
		for ( int i = 0; i < 4; ++i )
		{
			for ( int j = 0; j < 3; ++j )
			{
				if ( newFaceVerts[i][j] < 0 )
				{
					ss << "interior connectivity error at parent face "<< index << " child " << i << " vertex " << j ;
					statusMsg = OutputMessage(ss.str(), OutputMessage::errorPriority,"Faces::divideTriFace");
					log->logMessage(statusMsg);
					ss.str(s);
				}
			}			
		}
		//
		// DEBUG : check edge connectivity
		//
		for ( int i = 0; i < 4; ++i )
		{
			for ( int j = 0; j < 3; ++j )
			{
				if ( newFaceEdges[i][j] < 0 )
				{
					std::stringstream ss;
					ss << "interior connectivity error at parent face "<< index << " child " << i << " edge " << j ;
					statusMsg = OutputMessage(ss.str(), OutputMessage::errorPriority,"Faces::divideTriFace");
					log->logMessage(statusMsg);
					ss.str(s);
				}
			}			
		}
	}
	
	//
	// child face centers
	//
	//	note: for triangular faces, center particle of parent face become center particle of child 4
	//
	std::vector<xyzVector> newFaceCenters(4);
	std::vector<xyzVector> newFaceLagCenters(4);
	std::vector<xyzVector> tri(3);
	std::vector<xyzVector> lagTri(3);
	for ( int i = 0; i < 4; ++i )
	{
		tri[0] = particles.physCoord( newFaceVerts[i][0] );
		tri[1] = particles.physCoord( newFaceVerts[i][1] );
		tri[2] = particles.physCoord( newFaceVerts[i][2] );
		lagTri[0] = particles.lagCoord( newFaceVerts[i][0] );
		lagTri[1] = particles.lagCoord( newFaceVerts[i][1] );
		lagTri[2] = particles.lagCoord( newFaceVerts[i][2] );
		if ( gk == xyzVector::EuclideanGeometry )
		{
			newFaceCenters[i] = cartCentroid( tri );
			newFaceLagCenters[i] = cartCentroid( lagTri );
		}
		else if ( gk == xyzVector::SphericalGeometry )
		{
			newFaceCenters[i] = sphereCentroid( tri );
			newFaceLagCenters[i] = sphereCentroid( lagTri );
		}
	}

	k = particles.N();	
	for ( int i = 0; i < 3; ++i )
	{
		particles.insertParticle( newFaceCenters[i], newFaceLagCenters[i] );
		_centerParticle[_N+i] = k++;
		_vertices[_N+i] = newFaceVerts[i];
		_edges[_N+i] = newFaceEdges[i];
		_children[index][i] = _N + i;
		_parent[_N+i] = index;	
	}	
	particles.setPhysicalCoords( _centerParticle[index], newFaceCenters[3] );
	particles.setLagrangianCoords( _centerParticle[index], newFaceLagCenters[3] );
	_centerParticle[_N+3] = _centerParticle[index];
	_vertices[_N+3] = newFaceVerts[3];
	_edges[_N+3] = newFaceEdges[3];
	_children[index][3] = _N+3;
	_parent[_N+3] = index;	
	
	//
	// set child face areas
	//
	for ( int i = 0; i < 4; ++i )
		setArea( _N+i, particles, edges );
	
	_hasChildren[index] = true;
	_N += 4;
	_nActive += 3;
}

int Faces::countDivided() const
{
	int result = 0;
	for ( int i = 0; i < _nMax; ++i )
	{
		if ( _hasChildren[i] )
			++result;
	}
	return result;
};

void Faces::divideQuadFace( const int index, Particles& particles, Edges& edges )
{
	assert( _N+4 <= _nMax );
	assert( _edges[index].size() == 4 );
	
	OutputMessage statusMsg;
	std::string s;
	std::stringstream ss;
	
	std::vector< std::vector<int> > newFaceVerts(4, std::vector<int>(4, -1) );
	std::vector< std::vector<int> > newFaceEdges(4, std::vector<int>(4, -1) );
	//
	// connect parent vertices to child faces
	//
	for ( int i = 0; i < 4; ++i )
		newFaceVerts[i][i] = _vertices[index][i];
	//
	// loop over parent edges
	//
	int k = 0;
	int parentEdge, childEdge1, childEdge2;
	for ( int i = 0; i < 4; ++i)
	{
		parentEdge = _edges[index][i];
		
		if ( edges.hasChildren(parentEdge) ) // parent edge has already been divided
		{
			childEdge1 = edges.child1( parentEdge );
			childEdge2 = edges.child2( parentEdge );
			
			// connect child faces to child Edges
			if ( positiveEdge( index, parentEdge, edges) ) 
			{
				// parent edge has positive orientation
				newFaceEdges[k][i] = childEdge1;
				edges.setLeftFace( childEdge1, _N + k );
				newFaceEdges[(++k)%4][i] = childEdge2;
				edges.setLeftFace( childEdge2, _N + k%4 );
			}
			else
			{
				// parent edge has negative orientation
				newFaceEdges[k][i] = childEdge2;
				edges.setRightFace( childEdge2, _N + k );
				newFaceEdges[(++k)%4][i] = childEdge1;
				edges.setRightFace( childEdge1, _N + k%4 );
			}
			
		}
		else // parent edge needs to be divided
		{
			childEdge1 = edges.N();
			childEdge2 = edges.N() + 1;
			edges.divideEdge( parentEdge, particles );
			// connect child faces to child Edges
			if ( positiveEdge( index, parentEdge, edges) ) 
			{
				// parent edge has positive orientation
				newFaceEdges[k][i] = childEdge1;
				edges.setLeftFace( childEdge1, _N + k );
				newFaceEdges[(++k)%4][i] = childEdge2;
				edges.setLeftFace( childEdge2, _N + k%4 );
			}
			else
			{
				// parent edge has negative orientation
				newFaceEdges[k][i] = childEdge2;
				edges.setRightFace( childEdge2, _N + k );
				newFaceEdges[(++k)%4][i] = childEdge1;
				edges.setRightFace( childEdge1, _N + k%4 );
			}
		}
		// connect child faces to midpoints of parent edges
		newFaceVerts[i][(i+1)%4] = edges.dest( childEdge1 );
		newFaceVerts[(i+1)%4][i] = edges.dest( childEdge1 );
	}
	// connect center particle to child faces
	for ( int i = 0; i < 4; ++i )
		newFaceVerts[i][ (i+2)%4 ]  = _centerParticle[index];
	//
	// DEBUG : check vertex connectivity
	//
	if ( logLevel == OutputMessage::debugPriority )
	{
		for ( int i = 0; i < 4; ++i )
		{
			for ( int j = 0; j < 4; ++j )
			{
				if ( newFaceVerts[i][j] < 0 )
				{
					ss << "vertex connectivity error at parent face "<< index << " child " << i << " vertex " << j ;
					statusMsg = OutputMessage(ss.str(), OutputMessage::errorPriority,"Faces::divideTriFace");
					log->logMessage(statusMsg);
					ss.str(s);
				}
			}
		}
	}
	
	//
	// new interior edges
	//
	k = edges.N();
	edges.insertEdge( newFaceVerts[0][1], newFaceVerts[0][2], _N, _N+1, particles );
	newFaceEdges[0][1] = k;
	newFaceEdges[1][3] = k++;
	
	edges.insertEdge( newFaceVerts[3][1], newFaceVerts[3][2], _N+3, _N+2, particles );
	newFaceEdges[3][1] = k;
	newFaceEdges[2][3] = k++;
	
	edges.insertEdge( newFaceVerts[2][1], newFaceVerts[2][0], _N+1, _N+2, particles );
	newFaceEdges[1][2] = k;
	newFaceEdges[2][0] = k++;
	
	edges.insertEdge( newFaceVerts[0][2], newFaceVerts[0][3], _N, _N+3, particles );
	newFaceEdges[0][2] = k;
	newFaceEdges[3][0] = k;

	//
	// DEBUG : check edge connectivity
	//
	if ( logLevel == OutputMessage::debugPriority ) 
	{
		for ( int i = 0; i < 4; ++i )
		{
			for (int j = 0; j < 4; ++j)
			{
				if ( newFaceEdges[i][j] < 0 ) 
				{
					ss << "edge connectivity error at parent face "<< index << " child " << i << " edge " << j ;
					statusMsg = OutputMessage(ss.str(), OutputMessage::errorPriority,"Faces::divideTriFace");
					log->logMessage(statusMsg);
					ss.str(s);
				}	
			}
		}
	}
	
	//
	// child face centers
	//
	//	note : for quadrilateral panels, the center particle of the parent becomes a vertex particle for each child
	//
	std::vector<xyzVector> newFaceCenters(4);
	std::vector<xyzVector> newFaceLagCenters(4);
	std::vector<xyzVector> quad(4);
	std::vector<xyzVector> lagQuad(4);
	
	for ( int i = 0; i < 4; ++i)
	{
		for ( int j = 0; j < 4; ++j)
		{
			quad[j] = particles.physCoord( newFaceVerts[i][j] );
			lagQuad[j] = particles.lagCoord( newFaceVerts[i][j] );
		}
		if ( gk == xyzVector::EuclideanGeometry )
		{
			newFaceCenters[i] = cartCentroid(quad);
			newFaceLagCenters[i] = cartCentroid(lagQuad);
		}
		else if ( gk == xyzVector::SphericalGeometry )
		{
			newFaceCenters[i] = sphereCentroid(quad);
			newFaceLagCenters[i] = sphereCentroid(lagQuad);
		}
	}
	
	k = particles.N();
	for ( int i = 0; i < 4; ++i )
	{
		particles.insertParticle( newFaceCenters[i], newFaceLagCenters[i] );
		_centerParticle[_N+i] = k++;
		_vertices[_N+i] = newFaceVerts[i];
		_edges[_N+i] = newFaceEdges[i];
		_children[index][i] = _N+i;
		_parent[_N+i] = index;
	}
	//
	// set area
	//
	for ( int i = 0; i < 4; ++i)
		setArea(_N+i, particles, edges );
	setAreaToZero(index, particles);
	
	_hasChildren[index] = true;
	_N += 4;
	_nActive +=3;
};

void Faces::shrinkMemory()
{
	if ( _N < _nMax )
	{
		for ( int i = _N+1; i < _nMax; ++i )
		{
			_centerParticle.erase( _centerParticle.cend() );
			_vertices.erase( _vertices.cend() );
			_hasChildren.erase( _hasChildren.cend() );
			_children.erase( _children.cend() );
			_parent.erase( _parent.cend() );
		}
		
		_centerParticle.shrink_to_fit();
		_vertices.shrink_to_fit();
		_hasChildren.shrink_to_fit();
		_children.shrink_to_fit();
		_parent.shrink_to_fit();	
			
		_nMax = _N;
	}
};

void Faces::writePolygonsToVTK( std::ostream& fs ) const
{
	int nCells;
	int cellListSize;
	int nVerts;
	if ( _kind == triangularFaces )
	{
		nVerts = 3;
		nCells = nVerts * _nActive;
		cellListSize = 4 * nCells;
	}
	else if ( _kind == quadrilateralFaces )
	{
		nVerts = 4;
		nCells = nVerts * _nActive;
		cellListSize = 4 * nCells;
	}
	else
	{
		OutputMessage statusMsg("Polygon type not supported.", OutputMessage::errorPriority, "Faces::writePolygonsToVTK");
		log->logMessage(statusMsg);
		return;
	}
		
	fs << "POLYGONS " << nCells << "     " << cellListSize << std::endl;
	for ( int i = 0; i < _N; ++i )
	{
		if ( !_hasChildren[i] )
		{
			for ( int j = 0; j < nVerts; ++j )
			{
				fs << 3 << "     " << _vertices[i][j] << "     " << _vertices[i][(j+1)%nVerts]<< "     "  << _centerParticle[i] << std::endl;
			}
		}
	}
	
};

void Faces::writeVariablesToMatlab( std::ostream& fs ) const
{
	int nVerts;
	if ( _kind == triangularFaces )
	{
		nVerts = 3;
		fs << "faceKind = 'triangular';" << std::endl;
	}
	else if ( _kind == quadrilateralFaces )
	{
		nVerts = 4;
		fs << "faceKind = 'quadrilateral';" << std::endl;
	}
	
	fs << "faceN = " << _N << ";" << std::endl;
	fs << "faceNActive = " << _nActive << " ; " << std::endl;
	fs << "faceNMax = " << _nMax << ";" << std::endl;
	
		
	fs << "faceCenterParticle = [";
	for ( int i = 0; i < _N-1; ++i)
	{
		fs << _centerParticle[i] << " , " ;
	}
	fs << _centerParticle[_N-1] << " ];" << std::endl << std::endl;
	
	fs << "faceVertices = [";
	for ( int i = 0; i < _N-1; ++i)
	{
		for ( int j = 0; j < nVerts-1; ++j )
		{
			fs << _vertices[i][j] << " , ";
		}	
		fs << _vertices[i][nVerts - 1] << " ; ";
	}
	for ( int j = 0; j < nVerts-1; ++j )
	{
		fs << _vertices[_N-1][j] << " , ";
	}
	fs << _vertices[_N-1][nVerts - 1] << " ]; " << std::endl << std::endl;
	
	fs << "faceEdges = [";
	for ( int i = 0; i < _N-1; ++i)
	{
		for ( int j = 0; j < nVerts-1; ++j )
		{
			fs << _edges[i][j] << " , ";
		}	
		fs << _edges[i][nVerts - 1] << " ; ";
	}
	for ( int j = 0; j < nVerts-1; ++j )
	{
		fs << _edges[_N-1][j] << " , ";
	}
	fs << _edges[_N-1][nVerts - 1] << " ]; " << std::endl << std::endl;
	
	fs << "faceHasChildren = [";
	for ( int i = 0; i < _N-1; ++i)
	{
		if ( _hasChildren[i] )
			fs << "true" << " , ";
		else
			fs << "false" << " , ";
	}
	if ( _hasChildren[_N-1] )
		fs << "true" << " ]; " << std::endl << std::endl;
	else
		fs << "false" << " ]; " << std::endl << std::endl;
		
	fs << "faceChildren = [";
	for ( int i = 0; i < _N-1; ++i )
	{
		for ( int j = 0; j < 3; ++j )
		{
			fs << _children[i][j] << " , ";
		}
		fs << _children[i][3] << " ; ";
	}
	for ( int j = 0; j < 3; ++j )
	{
		fs << _children[_N-1][j] << " , ";
	}
	fs << _children[_N-1][3] << " ]; " << std::endl << std::endl;
	
	fs << "faceParent = [";
	for ( int i = 0; i < _N-1; ++i )
	{
		fs << _parent[i] << " , ";
	}
	fs << _parent[_N-1] << "];" << std::endl << std::endl;
};

std::vector< std::string > Faces::getInfo() const
{
	std::vector< std::string > info;
	std::stringstream ss;
	std::string s;
	
	ss << "   N = " << _N;
	info.push_back( ss.str() );
	ss.str(s);
	
	ss << "   nActive = " << _nActive;
	info.push_back(ss.str() );
	ss.str(s);
	
	ss << "    nMax = " << _nMax;
	info.push_back( ss.str() );
	ss.str(s);
	
	for ( int i = 0; i < _N; ++i )
	{
		ss << "face " << i << ": centerParticle = " << _centerParticle[i] << "; vertices = ";
		for ( int j = 0; j < _vertices[i].size(); ++j )
			ss << _vertices[i][j] << ", ";
		ss << "; edges = ";
		for ( int j = 0; j < _edges[i].size(); ++j)
			ss << _edges[i][j] << ", ";
		if ( _hasChildren[i] ) 
		{
			ss << "; DIVIDED, parent = " << _parent[i] << ", children = ";
			for ( int j = 0; j < 4; ++j)
				ss << _children[i][j] << ", ";
		}
		else 
		{
			ss << "; not divided, parent = " << _parent[i] ;
		}
		info.push_back( ss.str() );
		ss.str(s);
	}
	
	return info;
};

void Faces::setAreaToZero( const int index, Particles& particles ) const 
{
	particles.area[ _centerParticle[index] ] = 0.0;
};

bool Faces::positiveEdge( const int faceIndex, const int edgeIndex,  const Edges& edges ) const
{
	return (edges.leftFace(edgeIndex) == faceIndex);
};