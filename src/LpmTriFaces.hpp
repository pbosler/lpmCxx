#ifndef __LPM_TRI_FACES_HPP__
#define __LPM_TRI_FACES_HPP__

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
#include "LpmFaces.hpp"

using LpmXyzVector::XyzVector;

template <typename scalarType> class LpmTriFaces : public LpmFaces<scalarType> {

	int nVertsPerFace( const size_t index ) { return 3; }
	int nEdgesPerFace( const size_t index ) { return 3; }

	void divide( const int index, LpmParticles<scalarType>& particles, LpmEdges<scalarType>& edges ) {
		std::vector< std::vector<int> > newFaceVerts(4, std::vector<int>(3,-1) );
		std::vector< std::vector<int> > newFacesEdges(4, std::vector<int>(3,-1) );
		
		/*
			connect child faces to parent vertices
		*/
		for ( int i = 0; i < 3; ++i )
			newFaceVerts[i][i] = LpmFaces<scalarType>::_vertices[index][i];
			
		/*
			Loop over parent edges
		*/
	};

};



#endif
