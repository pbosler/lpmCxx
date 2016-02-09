//
//  RefinementType.cpp
//  LPM
//
//  Created by Peter Bosler on 2/7/16.
//  Copyright (c) 2014 Peter Bosler. All rights reserved.
//

#include "RefinementType.h"
#include "Field.h"
#include "PolyMesh2d.h"
#include "xyzVector.h"

bool ScalarIntegralRefinement::FlagFunction( const int faceIndex ) const
{
	const int pIndex = _refineMesh->faceParticleIndex( faceIndex );
	
	return ( std::abs(_refineField->scalar[pIndex]) * _refineMesh->particleArea(pIndex) > _tol );
}

bool ScalarVariationRefinement::FlagFunction( const int faceIndex ) const
{
	const int pIndex = _refineMesh->faceParticleIndex( faceIndex );
	
	const std::vector<int> faceVerts = _refineMesh->ccwVerticesAroundFace( faceIndex );
	
	double maxScalar = _refineField->scalar[pIndex];
	double minScalar = maxScalar;
	
	for (int i = 0; i < faceVerts.size(); ++i )
	{
		if ( _refineField->scalar[faceVerts[i]] > maxScalar)
			maxScalar = _refineField->scalar[faceVerts[i]];
		if ( _refineField->scalar[faceVerts[i]] < minScalar )
			minScalar = _refineField->scalar[faceVerts[i]];
	}
	
	return (( maxScalar - minScalar ) > _tol );
}

bool FlowMapVariationRefinement::FlagFunction( const int faceIndex ) const
{
	const int pIndex = _refineMesh->faceParticleIndex( faceIndex );
	
	xyzVector maxx0 = _refineMesh->particleLagCoord( pIndex );
	xyzVector minx0(maxx0);
	
	const std::vector<int> faceVerts = _refineMesh->ccwVerticesAroundFace( faceIndex );
	
	for ( int i = 0; i < faceVerts.size(); ++i)
	{
		const xyzVector testVec = _refineMesh->particleLagCoord( faceVerts[i] );
		if ( testVec.x > maxx0.x )
			maxx0.x = testVec.x;
		if ( testVec.y > maxx0.y )
			maxx0.y = testVec.y;
		if ( testVec.z > maxx0.z )
			maxx0.z = testVec.z;
		if ( testVec.x < minx0.x )
			minx0.x = testVec.x;
		if ( testVec.y < minx0.y )
			minx0.y = testVec.y;
		if ( testVec.z < minx0.z )
			minx0.z = testVec.z;
	}
	return ( ( maxx0.x - minx0.x + maxx0.y - minx0.y + maxx0.z - minx0.z ) > _tol );
}
