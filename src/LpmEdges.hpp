#ifndef __LPM_EDGES_HPP__
#define __LPM_EDGES_HPP__

#include <iostream>
#include <cmath>
#include <vector>
#include <sstream>
#include "OutputMessage.h"
#include "Logger.h"
#include "GlobalConstants.h"
#include "LpmXyzVector.hpp"
#include "LpmParticles.hpp"
#include "LpmCoords.hpp"

using LpmXyzVector::XyzVector;

template <typename scalarType> class LpmEdges {
	public :
		LpmEdges( const size_t nMax, const int nPts = 2, const int pRank = 0, const int nProcs = 1 ) :
			_nMax(nMax), _nPtsPerEdge(nPts), _nActive(0) {
			log = Logger::Instance( OutputMessage::debugPriority, pRank, nProcs);
			_orig.reserve(nMax);
			_dest.reserve(nMax);
			_leftFace.reserve(nMax);
			_rightFace.reserve(nMax);
			_parent.reserve(nMax);
			_child1.reserve(nMax);
			_child2.reserve(nMax);
		};
		virtual ~LpmEdges() {};
		
		size_t size() const { return _orig.size(); }
		int nActive() const { return _nActive; }
		int nMax() const { return _nMax; }
		int nDivided() const { return size() - _nActive; }
		
		int orig( const size_t index ) const { return _orig[index]; }
		int dest( const size_t index ) const { return _dest[index]; }
		int leftFace( const size_t index ) const { return _leftFace[index]; }
		int rightFace(const size_t index ) const { return _rightFace[index]; }
		
		std::vector<int> getLeafEdgesFromParent( const size_t index ) const { 
			return replaceEdgeIndicesWithChildren( std::vector<int>(1, index) );
		};
		
		std::vector<int> getAllVerticesAlongEdge( const size_t index ) const {
			std::vector<int> leafEdges = getLeafeEdgesFromParent( index );
			std::vector<int> result;
			for ( int i = 0; i < leafEdges.size(); ++i )
				result.push_back( _orig[ leafEdges[i] ] );
			result.push_back( _dest[leafEdges.back()] );
			return result;
		}
		
		scalarType length( const size_t index, const LpmParticles<scalarType>& particles ) const { 
			return particles.distance( _orig[index], _dest[index]); }
		
		XyzVector<scalarType> midpoint( const size_t index, const LpmParticles<scalarType>& particles ) const {
			return particles.midpoint( _orig[index], _dest[index] ); }
			
		XyzVector<scalarType> lagMidpoint( const size_t index, const LpmParticles<scalarType>& particles ) const {
			return particles.lagMidpoint( _orig[index], _dest[index] ); }
			
		XyzVector<scalarType> edgeVector( const size_t index, const LpmParticles<scalarType>& particles ) const {
			return particles.physCoordVec( _dest[index] ) - particles.physCoordVec( _orig[index] );
		};
		
		void setLeftFace( const size_t index, const size_t lFaceIndex ) {
			_leftFace[index] = lFaceIndex;
		};
		void setRightFace( const size_t index, const size_t rFaceIndex ) {
			_rightFace[index] = rFaceIndex;
		};
		void setNActive( const nA ) { _nActive = nA; }
		
		bool isDivided(const size_t index) const { return ( _child1[index] >= 0 ); }
		bool isActive(const size_t index) const { return !isDivided(index); }
		bool hasChildren(const size_t index) const { return isDivided(index); }
		
		int parent( const size_t index ) const { return _parent(index); }
		int child1( const size_t index ) const { return _child1(index); }
		int child2( const size_t index ) const { return _child2(index); }
		
		bool onBoundary( const size_t index) const { return ( _leftFace[index] < 0  || _rightFace[index] < 0 ); }
		
		virtual void insert(const int origIndex, const int destIndex, const int leftIndex, const int rightIndex) {
			_orig.push_back(origIndex);
			_dest.push_back(destIndex);
			_leftFace.push_back(leftIndex);
			_rightFace.push_back(rightIndex);
			_parent.push_back(-1);
			_child1.push_back(-1);
			_child2.push_back(-1);
		};
		
		virtual void divide( const int index, LpmParticles<scalarType>& particles ) {
			const XyzVector<scalarType> physMidpt = midpoint(index, particles);
			const XyzVector<scalarType> lagMidpt = lagMidpoint(index, particles);
			
			const int particleInsertPoint = particles.size();
			const int edgeInsertPoint = size();
			
			particles.insert(physMidpt, lagMidpt);
			
			_child1[index] = edgeInsertPoint;
			_child2[index] = edgeInsertPoint + 1;
			
			// child edge 1
			insert( _orig[index], particleInsertPoint, _left[index], _right[index] );
			_parent[edgeInsertPoint] = index;
			
			// child edge 2
			insert( particleInsertPoint, _dest[index], _left[index], _right[index] );
			_parent[edgeInsertPoint+1] = index;
			
			_nActive += 1;
		};
		

	
	protected :
		std::vector<int> _orig;
		std::vector<int> _dest;
		std::vector<int> _leftFace;
		std::vector<int> _rightFace;
		std::vector<int> _parent;
		std::vector<int> _child1;
		std::vector<int> _child2;
		
		std::vector<int> replaceEdgeIndicesWithChildren( std::vector<int> edgeList ) const {
			std::vector<int> result;
			bool keepGoing = false;
			for ( int i = 0; i < edgeList.size(); ++i ) {
				if ( hasChildren( edgeList[i] ) ) {
					result.push_back( _child1[edgeList[i]] );
					result.push_back( _child2[edgeList[i]] );
					if (hasChildren( _child1[edgeList[i]] ) || hasChildren( _child2[edgeList[i]] ) )
						keepGoing = true;
				}
				else {
					result.push_back( edgeList[i] );
				}
			}
			if ( !keepGoing ) 
				return result;
			else
				return replaceEdgeIndicesWithChildren( result );
		}
		
		int _nPtsPerEdge;
		int _nActive;
		int _nMax;
		
		Logger* log;
	
};



#endif
