#ifndef _LPM_TYPE_DEFS_H_
#define _LPM_TYPE_DEFS_H_

#include "LpmConfig.h"
#include <string>

namespace Lpm {

    /// Real number type
    typedef double scalar_type;
    
    /// Memory index type
    typedef long index_type;
    
    /// Pi
    static constexpr scalar_type PI = 3.1415926535897932384626433832795027975;
    
    /// Radians to degrees conversion factor
    static constexpr scalar_type rad2deg = 180.0 / PI;
    
    /// Gravitational acceleration
    static constexpr scalar_type g = 9.80616;
    
    /// Mean sea level radius of the Earth (meters)
    static constexpr scalar_type earthRadiusMeters = 6371220.0;
    
    /// One sidereal day, in units of seconds
    static constexpr scalar_type oneDaySeconds = 24.0 * 3600.0;
    
    /// Rotational rate of Earth about its z-axis
    static constexpr scalar_type earthOmega = 2.0 * PI / oneDaySeconds;
    
    /// Floating point zero
    static constexpr scalar_type ZERO_TOL = 1.0e-13;
    
    /// Kinds of geometry used by Lpm classes
    enum GeometryType {ONED_FREE, ONED_PERIODIC, PLANAR_GEOMETRY, SPHERICAL_SURFACE_GEOMETRY, CARTESIAN_3D_GEOMETRY};
    
    /// Basic face types
    enum FaceType {TRI, QUAD, QUAD_CUBIC, VORONOI};
    
    inline std::string facekindString(const FaceType& fkind) {
        switch (fkind) {
            case (TRI) : {
                return "TRI";
                break;
            }
            case (QUAD) : {
                return "QUAD";
                break;
            }
            case (QUAD_CUBIC) : {
                return "QUAD_CUBIC";
                break;
            }
            case (VORONOI) : {
                return "VORONOI";
                break;
            }
        }
    }
    
    inline std::string geometryString(const GeometryType& geom) {
        std::string result;
        switch (geom) {
            case (ONED_FREE) : {
                result = "ONED_FREE";
                break;
            }
            case (ONED_PERIODIC) : {
                result = "ONED_PERIODIC";
                break;
            }
            case (PLANAR_GEOMETRY) : {
                result = "PLANAR_GEOMETRY";
                break;
            }
            case (SPHERICAL_SURFACE_GEOMETRY) : {
                result = "SPHERICAL_SURFACE_GEOMETRY";
                break;
            }
            case (CARTESIAN_3D_GEOMETRY) : {
                result = "CARTESIAN_3D_GEOMETRY";
                break;
            }
        }
        return result;
    }
    
    inline std::string geom_weight_name(const GeometryType& geom)  {
        switch(geom) {
            case (ONED_FREE) : {
                return "length";
                break;
            }
            case (ONED_PERIODIC) : {
                return "length";
                break;
            }
            case (PLANAR_GEOMETRY) : {
                return "area";
                break;
            }
            case (SPHERICAL_SURFACE_GEOMETRY) : {
                return "area";
                break;
            }
            case (CARTESIAN_3D_GEOMETRY) : {
                return "volume";
                break;
            }
    }
  }
    
    /// Kinds of boundary conditions
    enum BoundaryCondition {FREE, DIRICHLET, NEUMANN, PERIODIC};

}
#endif
