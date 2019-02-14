#include "LpmTypeDefs.h"
#include <sstream>

namespace Lpm {

std::string facekindString(const FaceType& fkind) {
    std::string result;
    switch (fkind) {
            case (TRI) : {
                result = "TRI";
                break;
            }
            case (QUAD) : {
                result = "QUAD";
                break;
            }
            case (QUAD_CUBIC) : {
                result = "QUAD_CUBIC";
                break;
            }
            case (VORONOI) : {
                result = "VORONOI";
                break;
            }
        }
    return result;
}

std::string geometryString(const GeometryType& geom) {
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

std::string geom_weight_name(const GeometryType& geom) {
    std::string result;
    switch(geom) {
        case (ONED_FREE) : {
            result = "length";
            break;
        }
        case (ONED_PERIODIC) : {
            result = "length";
            break;
        }
        case (PLANAR_GEOMETRY) : {
            result = "area";
            break;
        }
        case (SPHERICAL_SURFACE_GEOMETRY) : {
            result = "area";
            break;
        }
        case (CARTESIAN_3D_GEOMETRY) : {
            result = "volume";
            break;
        }
    }
    return result;
}

}
