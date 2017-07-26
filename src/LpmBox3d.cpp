#include "LpmBox3d.h"
#include <sstream>

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

scalar_type Box3d::radius() const {
    const XyzVector cntd = centroid();
    std::vector<XyzVector> corners;
    corners.push_back(XyzVector(xmin, ymin, zmin));
    corners.push_back(XyzVector(xmin, ymax, zmin));
    corners.push_back(XyzVector(xmin, ymin, zmax));
    corners.push_back(XyzVector(xmin, ymax, zmax));
    corners.push_back(XyzVector(xmax, ymin, zmin));
    corners.push_back(XyzVector(xmax, ymax, zmin));
    corners.push_back(XyzVector(xmax, ymin, zmax));
    corners.push_back(XyzVector(xmax, ymax, zmax));
    scalar_type result = 0.0;
    for (int i = 0; i < 8; ++i) {
        const scalar_type testDist = distance(cntd, corners[i]);
        if (testDist > result)
            result = testDist;
    }
    return result;
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