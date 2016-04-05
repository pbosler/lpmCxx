#ifndef __LPM_FACES_HPP__
#define __LPM_FACES_HPP__

#include <iostream>
#include <cmath>
#include <vector>
#include <sstream>
#include "OutputMessage.h"
#include "Logger.h"
#include "GlobalConstants.h"
#include "LpmXyzVector.hpp"
#include "LpmParticles.hpp"
#include "LpmEdges.hpp"

using LpmXyzVector::XyzVector;
/*
Base class for LPM Faces; will be subclassed for Trianguler, Quadrilateral, and Voronoi faces.
*/
template<typename scalarType> class LpmFaces {
	public :
	
		virtual ~LpmFaces() {};
		
		virtual void divide( const int index, LpmParticles<scalarType>& particles, LpmEdges<scalarType>& edges ) = 0;
		virtual void computeFaceArea( const int index, LpmParticles<scalarType>& particles, 
									  const LpmEdges<scalarType>& edges) const {
			scalarType area = 0.0;
			const int indexB = _centerParticle[index];
			const int nVerts = _vertices[index].size();
			for ( int i = 0; i < nVerts; ++i ) {
				const int indexA = _vertices[index][i];
				const int indexC = _vertices[index][ (i+1)%nVerts ];
				area += particles.triArea(indexA, indexB, indexC);
			}
			particles.setArea( _centerParticle[index], area );
		}
		
		scalarType faceArea( const int index, const LpmParticles<scalarType>& particles ) const { 
			return particles.area( _centerParticle[index] ); }
		
		void setFaceArea( const int index, const scalarType newArea, LpmParticles<scalarType>& particles ) const { 
			particles.setArea( _centerParticle[index], newArea ); }
			
		int centerParticle( const size_t index ) const { return _centerParticle[index]; }
		std::vector<int> faceVertices( const size_t index ) const { return _vertices[index]; }
		std::vector<int> faceEdges( const size_t index ) const { return _edges[index]; }
		
		virtual void insert( const int cntrParticle, const std::vector<int>& vertIndices, 
							 const std::vector<int>& edgeIndices ) {
			_centerParticle.push_back(cntrParticle);
			_vertices.push_back(vertIndices);
			_edges.push_back(edgeIndices);
			_hasChildren.push_back( false );
			_children.push_back( std::vector<int>(4,-1) );
			_parent.push_back(-1);
		}; 
		
		void setVertices( const size_t index, const std::vector<int>& vertIndices ) { _vertices[index] = vertIndices; }
		
		void setEdges( const size_t index, const std::vector<int>& edgeIndices ) { _edges[index] = edgeIndices; }
		
		size_t size() const { return _centerParticle.size(); }
		int nActive() const { return _nActive; }
		bool hasChildren( const size_t index ) const { return _hasChildren[index]; }
		bool isDivided( const size_t index ) const { return _hasChildren[index]; }
		bool isActive( const size_t index ) const { return !_hasChildren[index]; }
		
		virtual int nVertsPerFace( const size_t index ) = 0;
		virtual int nEdgesPerFace( const size_t index ) = 0;
		
	protected :
		LpmFaces( const int nMax, const int nMaxVerticesPerFace, const int nMaxEdgesPerFace, const int pRank = 0, 
				  const int nProcs = 1 ) : _nMax(nMax), _nActive(0), _nMaxVerts(nMaxVerticesPerFace),
				  _nMaxEdges(nMaxEdgesPerFace) {
				  log = Logger::Instance( OutputMessage::debugPriority, pRank, nProcs );
				  
				  _centerParticle.reserve(nMax);
				  _hasChildren.reserve(nMax);
				  _parent.reserve(nMax);			  
				  }
	
		std::vector<int> _centerParticle;
		std::vector< std::vector<int> > _vertices;
		std::vector< std::vector<int> > _edges;
		std::vector<bool> _hasChildren;
		std::vector< std::vector<int> > _children;
		std::vector<int> _parent;
		int _nActive;
		int _nMax;
		int _nMaxVerts;
		int _nMaxEdges;
		
		Logger* log;
};

#endif
