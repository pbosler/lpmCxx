//
//  RefinementType.h
//  LPM
//
//  Created by Peter Bosler on 2/7/16.
//  Copyright (c) 2014 Peter Bosler. All rights reserved.
//

#ifndef __LPM__RefinementType__
#define __LPM__RefinementType__

#include <vector>
#include <iostream>
#include "Field.h"
#include "PolyMesh2d.h"
#include "OutputMessage.h"
#include "Logger.h"


class RefinementType
{
	public:
		RefinementType( const PolyMesh2d* meshPtr, const Field* fieldPtr, const double tol ) : 
			_refineMesh( meshPtr ), _refineField( fieldPtr ), _tol(tol) {};
	
		virtual bool FlagFunction( const int faceIndex ) const = 0;
	
	protected:
		const PolyMesh2d* _refineMesh;
		const Field* _refineField;
		const double _tol;
};

class ScalarIntegralRefinement : public RefinementType
{
	public :
		virtual bool FlagFunction( const int faceIndex ) const;
};

class ScalarVariationRefinement : public RefinementType
{
	public :
		virtual bool FlagFunction(const int faceIndex ) const;
};

class FlowMapVariationRefinement : public RefinementType
{
	public :
		virtual bool FlagFunction( const int faceIndex ) const;
};

#endif
