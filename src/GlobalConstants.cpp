//
//  GlobalConstants.cpp
//  LPM
//
//  Created by Peter Bosler on 10/31/14.
//  Copyright (c) 2014 Peter Bosler. All rights reserved.
//

#include "GlobalConstants.h"
/**	@brief GlobalConstants implementation
	@file
	@author Peter Bosler <pabosle@sandia.gov>
 */
GlobalConstants* GlobalConstants::_instance = 0;

/**
 @brief Returns a pointer to the single GlobalConstants object on each process. 
 */
GlobalConstants* GlobalConstants::Instance()
{
    if ( _instance == 0 )
        _instance = new GlobalConstants();
    return _instance;
}

/** @brief Private constructor.
*/
GlobalConstants::GlobalConstants()
{
    _pi = 3.1415926535897932384626433832795027975;
    _rad2deg = 180.0/_pi;
    _deg2rad = _pi/180.0;
    _gravity = 9.80616;
    _oneDay = 86140.0;
    _earthRadius = 6371220.0;
    _Omega = 2.0 * _pi / _oneDay;
    _zeroTol = 1.0e-14;
    _domainRadius = 1.0;
}