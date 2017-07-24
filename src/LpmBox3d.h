#ifndef _LPM_BOX_3D_H_
#define _LPM_BOX_3D_H_

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmXyzVector.h"

namespace Lpm {

struct Box3d {
    Box3d() : xmin(0.0), xmax(0.0), ymin(0.0), ymax(0.0), zmin(0.0), zmax(0.0) {};
    Box3d(const scalar_type x0, const scalar_type x1, const scalar_type y0, const scalar_type y1, const scalar_type z0,
        const scalar_type z1) : xmin(x0), xmax(x1), ymin(y0), ymax(y1), zmin(z0), zmax(z1) {};
    
    inline scalar_type volume() const {return (xmax - xmin) * (ymax - ymin) * (zmax - zmin);}
    inline scalar_type area2d() const {return (xmax - xmin) * (ymax - ymin);}
    
    inline XyzVector centroid() const {return XyzVector(0.5 * (xmax + xmin), 0.5 * (ymax + ymin), 0.5 * (zmax + zmin));}
    
    inline bool containsPoint(const XyzVector& vec) const {return (xmin <= vec.x && vec.x < xmax) &&
                                                                  (ymin <= vec.y && vec.y < ymax) &&
                                                                  (zmin <= vec.z && vec.z < zmax);}
    
    scalar_type longestEdge() const;
    scalar_type shortestEdge() const;
    scalar_type aspectRatio() const;
    
    scalar_type edgeLength(const int dim) const;
    
    std::vector<Box3d> bisectAll() const;
    std::vector<Box3d> bisectAlongDims(const bool* dims) const;
    
    std::string infoString() const;
    
    scalar_type radius() const;
    
    /// minimum x-coordinate of this box
    scalar_type xmin;
    /// maximum x-coordinate of this box
    scalar_type xmax;
    /// minimum y-coordinate of this box
    scalar_type ymin;
    /// maximum y-coordinate of this box
    scalar_type ymax;
    /// minimum z-coordinate of this box
    scalar_type zmin;
    /// maximum z-coordinate of this box
    scalar_type zmax;
};

}

#endif
