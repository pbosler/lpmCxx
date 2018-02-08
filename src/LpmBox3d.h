#ifndef _LPM_BOX_3D_H_
#define _LPM_BOX_3D_H_

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmXyzVector.h"

#define BOX_PADDING_FACTOR 0.00001

namespace Lpm {

struct Box3d {
    Box3d() : xmin(0.0), xmax(0.0), ymin(0.0), ymax(0.0), zmin(0.0), zmax(0.0) {};
    Box3d(const scalar_type x0, const scalar_type x1, const scalar_type y0, const scalar_type y1, const scalar_type z0,
        const scalar_type z1) : xmin(x0), xmax(x1), ymin(y0), ymax(y1), zmin(z0), zmax(z1) {};
    Box3d(const scalar_type maxr, bool pad=true) : xmin(-maxr), xmax(maxr), ymin(-maxr), ymax(maxr), zmin(-maxr), zmax(maxr) {
        if (pad) {
            xmin -= BOX_PADDING_FACTOR;
            xmax += BOX_PADDING_FACTOR;
            ymin -= BOX_PADDING_FACTOR;
            ymax += BOX_PADDING_FACTOR;
            zmin -= BOX_PADDING_FACTOR;
            zmax += BOX_PADDING_FACTOR;
        }
    };
    
    inline scalar_type volume() const {return (xmax - xmin) * (ymax - ymin) * (zmax - zmin);}
    inline scalar_type area2d() const {return (xmax - xmin) * (ymax - ymin);}
    
    inline XyzVector centroid() const {return XyzVector(0.5 * (xmax + xmin), 0.5 * (ymax + ymin), 0.5 * (zmax + zmin));}
    
    inline bool containsPoint(const XyzVector& vec) const {return (xmin <= vec.x && vec.x < xmax) &&
                                                                  (ymin <= vec.y && vec.y < ymax) &&
                                                                  (zmin <= vec.z && vec.z < zmax);}
    
    inline XyzVector mins() const {return XyzVector(xmin, ymin, zmin);}
    inline XyzVector maxs() const {return XyzVector(xmax, ymax, zmax);}
    
    std::vector<XyzVector> faceCentroids() const;
    
    // returns the point inside the box closest to a query point
    XyzVector closestPointInBox(const XyzVector& query = XyzVector()) const;
    // returns the point inside a box farthest from a query point
    XyzVector farthestPointInBox(const XyzVector& query = XyzVector()) const;
    
    bool intersectsSphere(const XyzVector& center = XyzVector(), const scalar_type sphRadius = 1.0) const;
    bool containsSphere(const XyzVector& center = XyzVector(), const scalar_type sphRadius = 1.0) const;
    
    std::vector<XyzVector> corners() const;
    
    scalar_type longestEdge() const;
    scalar_type shortestEdge() const;
    scalar_type aspectRatio() const;
    
    scalar_type edgeLength(const int dim) const;
    
    std::vector<Box3d> bisectAll() const;
    std::vector<Box3d> bisectAlongDims(const bool* dims) const;
    std::vector<Box3d> bisectAlongDims(const std::vector<bool>& dims) const;
    
    std::string infoString() const;
    
    scalar_type maxRadius() const;
    scalar_type minRadius() const;
    
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
