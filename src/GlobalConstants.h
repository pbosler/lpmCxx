//
//  GlobalConstants.h
//  LPM
//
//  Created by Peter Bosler on 10/31/14.
//  Copyright (c) 2014 Peter Bosler. All rights reserved.
//

#ifndef __LPM__GlobalConstants__
#define __LPM__GlobalConstants__

/**	@file GlobalConstants.h
	@brief GlobalConstants header file
	@author Peter Bosler <pabosle@sandia.gov>
*/

/** 
    @class GlobalConstants
    @brief A singleton object for storing and accessing global physical constants.
    All clients must access this class via the GlobalConstants::Instance() member function.
 */
class GlobalConstants
{
public:
    static GlobalConstants* Instance();
    
    /// Pi
    inline double Pi() const { return _pi;}
    
    /// conversion factor, radians to degrees
    inline double Rad2Deg() const { return _rad2deg;}
    
    /// conversion factor, degrees to radians
    inline double Deg2Rad() const { return _deg2rad;}
    
    /// returns the constant g
    inline double Gravity() const { return _gravity;}
    
    /// returns the length of one day, in seconds
    inline double OneDayInSeconds() const { return _oneDay;}
    
    /// returns the radius of Earth
    inline double EarthRadius() const {return _earthRadius;}
    
    inline double DomainRadius() const { return _domainRadius; }
    
    /// returns the Earth's angular velocity
    inline double Omega() const { return _Omega; }
    
    /// returns the zero tolerance used for comparing real numbers
    inline double ZeroTol() const { return _zeroTol; }
    
    /** @brief Set the constant g.  Default is 9.8106 m s^{-2}.
        @param newGrav
        */
    inline void SetGravity( const double newGrav ){ _gravity = newGrav; }
    
    
    /** @brief Set the period of the sphere's rotation.  Default value is one sidereal day, in seconds.
    This function updates the angular velocity of the sphere @ref _Omega to match the input value.
        @param newDay
        */
    inline void SetSiderealDay( const double newDay)
    {
        _oneDay = newDay;
        _Omega = 2.0 * _pi / _oneDay;
    }
    
    /** @brief Set the angular velocity of the sphere.  Does not affect the day length.
    	Default value is appropriate for an Earth-sized sphere, with units of 1/seconds. 
        @param newOmega
        */
    inline void SetOmega( const double newOmega ){ _Omega = newOmega; }
    
    /** @brief Set the radius of the sphere.  Default is Earth's mean radius in meters.
    	@param newRad
    	*/
    inline void SetEarthRadius( ){ _domainRadius = _earthRadius; }
    
    inline void SetDomainRadius( const double radius ){ _domainRadius = radius; }
    
    /** @brief set the tolerance for comparing equality of two real numbers
    	@param newTol
    */
    inline void SetZeroTol( const double newTol ){ _zeroTol = newTol; }
    
protected:
    GlobalConstants();
    
private:

    static GlobalConstants* _instance;
   
     double _pi;
     double _rad2deg;
     double _deg2rad;
     double _gravity;
     double _Omega;
     double _earthRadius;
     double _oneDay;
     double _zeroTol;
     double _domainRadius;
};

#endif /* defined(__LPM__GlobalConstants__) */
