#include "LpmBox3d.h"
#include <sstream>
#include <limits>

namespace Lpm {

scalar_type Box3d::longestEdge() const {
    scalar_type result = xmax - xmin;
    if (ymax-ymin > result)
        result = ymax - ymin;
    if (zmax-zmin > result)
        result = zmax - zmin;
    return result;
}

scalar_type Box3d::shortestEdge() const {
    const scalar_type dx = xmax - xmin;
    const scalar_type dy = ymax - ymin;
    const scalar_type dz = zmax - zmin;
    scalar_type result = dx;
    if (dy < dx)
        result = dy;
    if (dz > 0.0 && dz < result)
        result = dz;
    return result;
}

scalar_type Box3d::aspectRatio() const {
    return longestEdge() / shortestEdge();
}

scalar_type Box3d::maxRadius() const {
    const XyzVector cntd = centroid();
    const std::vector<XyzVector> crnrs = corners();
    scalar_type result = 0.0;
    for (int i = 0; i < 8; ++i) {
        const scalar_type testDist = distance(cntd, crnrs[i]);
        if (testDist > result)
            result = testDist;
    }
    return result;
}

scalar_type Box3d::minRadius() const {
    const XyzVector cntd = centroid();
    const std::vector<XyzVector> face_ctds = faceCentroids();
    scalar_type result = std::numeric_limits<scalar_type>::max();
    for (int i=0; i<6; ++i) {
        const scalar_type testDist = distance(cntd, face_ctds[i]);
        if (testDist < result)
            result = testDist;
    }
    return result;
}


std::vector<XyzVector> Box3d::corners() const {
    std::vector<XyzVector> result(8);
    result[0] = XyzVector(xmin, ymin, zmin);
    result[1] = XyzVector(xmin, ymax, zmin);
    result[2] = XyzVector(xmin, ymin, zmax);
    result[3] = XyzVector(xmin, ymax, zmax);
    result[4] = XyzVector(xmax, ymin, zmin);
    result[5] = XyzVector(xmax, ymax, zmin);
    result[6] = XyzVector(xmax, ymin, zmax);
    result[7] = XyzVector(xmax, ymax, zmax);
    return result;
}

std::vector<XyzVector> Box3d::faceCentroids() const {
    std::vector<XyzVector> result(6);
    result[0] = XyzVector(0.5*(xmin + xmax), ymin, 0.5*(zmin+zmax));
    result[1] = XyzVector(xmax, 0.5*(ymin + ymax), 0.5*(zmin+zmax));
    result[2] = XyzVector(0.5*(xmin + xmax), ymax, 0.5*(zmin+zmax));
    result[3] = XyzVector(xmin, 0.5*(ymin + ymax), 0.5*(zmin+zmax));
    result[4] = XyzVector(0.5*(xmin+xmax), 0.5*(ymin+ymax), zmin);
    result[5] = XyzVector(0.5*(xmin+xmax), 0.5*(ymin+ymax), zmax);
    return result;
}

XyzVector Box3d::closestPointInBox(const XyzVector& query) const {
/*  returns the closest point inside the box to query point.
*/
    XyzVector rvec;
    if (containsPoint(query)) {
        rvec = query;
    }
    else {
        const std::vector<scalar_type> ll = (mins() - query).getStdVector();
        const std::vector<scalar_type> uu = (maxs() - query).getStdVector();
        std::vector<scalar_type> result(3, 0.0);
        for (int i=0; i<3; ++i) {
            if (ll[i] == uu[i]) {
                result[i] = ll[i];
            }
            else {
                if (ll[i] >= 0.0) {
                    result[i] = ll[i];
                }
                else if (uu[i] <= 0.0) {
                    result[i] = uu[i];
                }
                else {
                    result[i] = 0.0;
                }
            }
        }
        rvec = XyzVector(result) + query;
    }
    return rvec;
}

XyzVector Box3d::farthestPointInBox(const XyzVector& query) const {
    const std::vector<XyzVector> crnrs = corners();
    scalar_type dist = 0.0;
    int corner_index = 0;
    for (int i=0; i<8; ++i) {
        const scalar_type testDist = distance(crnrs[i], query);
        if (testDist > dist) {
            dist = testDist;
            corner_index = i;
        }
    }
    return crnrs[corner_index];
}

bool Box3d::intersectsSphere(const XyzVector& center, const scalar_type sphRadius) const {
    const XyzVector closest = closestPointInBox(center);
    const XyzVector farthest = farthestPointInBox(center);
    return (closest-center).magnitude() <= sphRadius && (farthest-center).magnitude() >= sphRadius;
}

bool Box3d::containsSphere(const XyzVector& center, const scalar_type sphRadius) const {
    bool result = false;
    if (containsPoint(center)) {
        const XyzVector closest = closestPointInBox(center);
        const XyzVector farthest = farthestPointInBox(center);
        
    }
}

std::vector<Box3d> Box3d::bisectAll() const {
    std::vector<Box3d> result;
    result.push_back(Box3d(xmin, 0.5 * (xmin + xmax), ymin, 0.5 * (ymin + ymax), zmin, 0.5 * (zmin + zmax)));
    result.push_back(Box3d(0.5 * (xmin + xmax), xmax, ymin, 0.5 * (ymin + ymax), zmin, 0.5 * (zmin + zmax)));
    result.push_back(Box3d(xmin, 0.5 * (xmin + xmax), 0.5 * (ymin + ymax), ymax, zmin, 0.5 * (zmin + zmax)));
    result.push_back(Box3d(0.5 * (xmin + xmax), xmax, 0.5 * (ymin + ymax), ymax, zmin, 0.5 * (zmin + zmax)));
    result.push_back(Box3d(xmin, 0.5 * (xmin + xmax), ymin, 0.5 * (ymin + ymax), 0.5 * (zmin + zmax), zmax));
    result.push_back(Box3d(0.5 * (xmin + xmax), xmax, ymin, 0.5 * (ymin + ymax), 0.5 * (zmin + zmax), zmax));
    result.push_back(Box3d(xmin, 0.5 * (xmin + xmax), 0.5 * (ymin + ymax), ymax, 0.5 * (zmin + zmax), zmax));
    result.push_back(Box3d(0.5 * (xmin + xmax), xmax, 0.5 * (ymin + ymax), ymax, 0.5 * (zmin + zmax), zmax));
    return result;
}

std::vector<Box3d> Box3d::bisectAlongDims(const bool* dims) const {
    std::vector<Box3d> result;
    int dimCount = 0;
    for (int i = 0; i < 3; ++i) {
        if (dims[i])
            dimCount += 1;
    }
    if (dimCount == 1) {
        if (dims[0]) {
            result.push_back(Box3d(xmin, 0.5 * (xmin + xmax), ymin, ymax, zmin, zmax));
            result.push_back(Box3d(0.5 * (xmin + xmax), xmax, ymin, ymax, zmin, zmax));
        }
        else if (dims[1]) {
            result.push_back(Box3d(xmin, xmax, ymin, 0.5 * (ymin + ymax), zmin, zmax));
            result.push_back(Box3d(xmin, xmax, 0.5 * (ymin + ymax), ymax, zmin, zmax));
        }
        else if (dims[2]) {
            result.push_back(Box3d(xmin, xmax, ymin, ymax, zmin, 0.5 * (zmin + zmax)));
            result.push_back(Box3d(xmin, xmax, ymin, ymax, 0.5 * (zmin + zmax), zmax));
        }
    }
    else if (dimCount == 2) {
        if (dims[0] && dims[1]) {
            result.push_back(Box3d(xmin, 0.5 * (xmin + xmax), ymin, 0.5 * (ymin + ymax), zmin, zmax));
            result.push_back(Box3d(0.5 * (xmin + xmax), xmax, ymin, 0.5 * (ymin + ymax), zmin, zmax));
            result.push_back(Box3d(xmin, 0.5 * (xmin + xmax), 0.5 * (ymin + ymax), ymax, zmin, zmax));
            result.push_back(Box3d(0.5 * (xmin + xmax), xmax, 0.5 * (ymin + ymax), ymax, zmin, zmax));
        }
        else if (dims[0] && dims[2]) {
            result.push_back(Box3d(xmin, 0.5 * (xmin + xmax), ymin, ymax, zmin, 0.5 * (zmin + zmax)));
            result.push_back(Box3d(0.5 * (xmin + xmax), xmax, ymin, ymax, zmin, 0.5 * (zmin + zmax)));
            result.push_back(Box3d(xmin, 0.5 * (xmin + xmax), ymin, ymax, 0.5 * (zmin + zmax), zmax));
            result.push_back(Box3d(0.5 * (xmin + xmax), xmax, ymin, ymax, 0.5 * (zmin + zmax), zmax));
        }
        else if (dims[1] && dims[2]) {
            result.push_back(Box3d(xmin, xmax, ymin, 0.5 * (ymin + ymax), zmin, 0.5 * (zmin + zmax)));
            result.push_back(Box3d(xmin, xmax, 0.5 * (ymin + ymax), ymax, zmin, 0.5 * (zmin + zmax)));
            result.push_back(Box3d(xmin, xmax, ymin, 0.5 * (ymin + ymax), 0.5 * (zmin + zmax), zmax));
            result.push_back(Box3d(xmin, xmax, 0.5 * (ymin + ymax), ymax, 0.5 * (zmin + zmax), zmax));
        }
    }
    else {
        result = bisectAll();
    }
    return result;
}

scalar_type Box3d::edgeLength(const int dim) const {
    scalar_type result = 0.0;
    switch (dim) {
        case 0 : {
            result = xmax - xmin;
            break;
        }
        case 1 : {
            result = ymax - ymin;
            break;
        }
        case 2 : {
            result = zmax - zmin;
            break;
        }
    }
    return result;
}

std::string Box3d::infoString() const {
    std::stringstream ss;
    ss << "x in [" << xmin << ", " << xmax << "], y in [" << ymin << ", " << ymax << "], z in [" << zmin << ", " << zmax << "]" << std::endl;
    return ss.str();
}

}