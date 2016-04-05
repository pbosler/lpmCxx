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
		virtual ~LpmEdges() {};
		
		size_t size() const { return _orig.size(); }
		int nActive() const { return _nActive; }
		int nMax() const { return _nMax; }
		int nDivided() const { return size() - _nActive; }
		
		int orig( const size_t index ) const { return _orig[index]; }
		int dest( const size_t index ) const { return _dest[index]; }
		int leftFace( const size_t index ) const { return _leftFace[index]; }
		int rightFace(const size_t index ) const { return _rightFace[index]; }
		
		void setLeftFace( const size_t index, const size_t lFaceIndex ) {
			_leftFace[index] = lFaceIndex;
		};
		void setRightFace( const size_t index, const size_t rFaceIndex ) {
			_rightFace[index] = rFaceIndex;
		};
		
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
		};
		virtual void divide() {};
		
		virtual scalarType edgeLenght( const size_t index ) const = 0;
	
	protected :
		LpmEdges( const size_t nMax, const int pRank = 0, const int nProcs = 1 ) :
			_nMax(nMax), _nPtsPerEdge(2), _nActive(0) {
			log = Logger::Instance( OutputMessage::debugPriority, pRank, nProcs);
			_orig.reserve(nMax);
			_dest.reserve(nMax);
			_leftFace.reserve(nMax);
			_rightFace.reserve(nMax);
			_parent.reserve(nMax);
			_child1.reserve(nMax);
			_child2.reserve(nMax);
		};
	
		std::vector<int> _orig;
		std::vector<int> _dest;
		std::vector<int> _leftFace;
		std::vector<int> _rightFace;
		std::vector<int> _parent;
		std::vector<int> _child1;
		std::vector<int> _child2;
		
		int _nPtsPerEdge;
		int _nActive;
		int _nMax;
		
		Logger* log;
	
};



#endif
