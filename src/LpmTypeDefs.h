#ifndef _LPM_TYPE_DEFS_H_
#define _LPM_TYPE_DEFS_H_

#ifdef USE_NANOFLANN
#include "nanoflann.hpp"
#endif

namespace Lpm {

    /// Real number type
    typedef double scalar_type;
    
    /// Memory index type
    typedef long index_type;
    
    /// Pi
    static const scalar_type PI = 3.1415926535897932384626433832795027975;
    
    /// Radians to degrees conversion factor
    static const scalar_type rad2deg = 180.0 / PI;
    
    /// Gravitational acceleration
    static const scalar_type g = 9.80616;
    
    /// Mean sea level radius of the Earth (meters)
    static const scalar_type earthRadiusMeters = 6371220.0;
    
    /// One sidereal day, in units of seconds
    static const scalar_type oneDaySeconds = 24.0 * 3600.0;
    
    /// Rotational rate of Earth about its z-axis
    static const scalar_type earthOmega = 2.0 * PI / oneDaySeconds;
    
    /// Floating point zero
    static const scalar_type ZERO_TOL = 1.0e-13;
    
    /// Kinds of geometry used by Lpm classes
    enum GeometryType {PLANAR_GEOMETRY, SPHERICAL_SURFACE_GEOMETRY, CARTESIAN_3D_GEOMETRY};
    
    /// Kinds of boundary conditions
    enum BoundaryCondition {FREE, DIRICHLET, NEUMANN, PERIODIC};

}
#endif
