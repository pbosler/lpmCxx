//
//  Refinement2d.h
//  LPM
//
//  Created by Peter Bosler on 2/7/16.
//  Copyright (c) 2014 Peter Bosler. All rights reserved.
//

#ifndef __LPM__Refinement2d__
#define __LPM__Refinement2d__

#include <vector>
#include <iostream>
#include "Field.h"
#include "PolyMesh2d.h"
#include "OutputMessage.h"
#include "Logger.h"
#include "RefinementType.h"

class Refinement2d
{
	public:
		Refinement2d( const PolyMesh2d& aMesh );
			
	protected:
		std::vector<bool> _refineFlag;
};

#endif
