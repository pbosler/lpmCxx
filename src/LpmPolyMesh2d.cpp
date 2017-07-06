#include "LpmPolyMesh2D.h"
#include "LpmEuclideanCoords.h"
#include "LpmSphericalCoords.h"
#include "LpmOutputMessage.h"
#include "LpmTriFaces.h"
#include "LpmQuadFaces.h"
#include <typeinfo>
#include <limits>
#include <iostream>
#include <fstream>
#include <exception>
#include <sstream>
#include <mpi.h>

namespace Lpm {

std::unique_ptr<Logger> PolyMesh2d::log(new Logger(OutputMessage::debugPriority, "PolyMesh2d_log"));

PolyMesh2d::PolyMesh2d(MeshSeed& seed, const int maxRecursionLevel, const bool isLagrangian, 
    const scalar_type domainRadius, const int prank) : lagrangian(isLagrangian), nRootFaces(seed.nRootFaces()) {
    std::stringstream ss;
    ss << "building" << (lagrangian ? " Lagrangian " : " ") << "PolyMesh2d with seed: " << seed.idString() << 
        ", recursion level " << maxRecursionLevel << ", and radius " << domainRadius;
    OutputMessage statusMsg(ss.str(), OutputMessage::tracePriority, "PolyMesh2d::PolyMesh2d");
    log->logMessage(statusMsg);
    log->startSection();
    const index_type nMaxVerts = seed.nVertices(maxRecursionLevel);
    index_type nMaxFaces = 0;
    index_type nMaxEdges = 0;
    for (int k = 0; k <= maxRecursionLevel; ++k) {
        nMaxFaces += seed.nFaces(k);
        nMaxEdges += seed.nEdges(seed.nVertices(k), seed.nFaces(k));
    }
    
    if (typeid(seed) == typeid(IcosTriSphereSeed) || typeid(seed) == typeid(CubedSphereSeed)) {
        coords = std::shared_ptr<Coords>(new SphericalCoords(nMaxVerts, domainRadius));
    }
    else if (typeid(seed) == typeid(TriHexSeed) || typeid(seed) == typeid(QuadRectSeed)) {
        coords = std::shared_ptr<Coords>(new EuclideanCoords(nMaxVerts, PLANAR_GEOMETRY));
    }
    coords->setLogProc(prank);
    
    edges = std::shared_ptr<Edges>(new Edges(nMaxEdges, coords));
    edges->setLogProc(prank);
    
    if (typeid(seed) == typeid(IcosTriSphereSeed) || typeid(seed) == typeid(TriHexSeed)) {
        faces = std::shared_ptr<Faces>(new TriFaces(nMaxFaces, edges, coords));
    }
    else if  (typeid(seed) == typeid(QuadRectSeed) || typeid(seed) == typeid(CubedSphereSeed)) {
        faces = std::shared_ptr<Faces>(new QuadFaces(nMaxFaces, edges, coords));
    }
    faces->setLogProc(prank);
    
    seed.initMeshFromSeed(coords, edges, faces, domainRadius);
//     ss.str(std::string());
// //    ss << "seed.infoString() = " << seed.infoString();
//     statusMsg.resetMsgString(ss.str());
//     log->logMessage(statusMsg);
    
    for (int k = 0; k < maxRecursionLevel; ++k) {
        const index_type nfaces = faces->n();
        for (index_type i = 0; i < nfaces; ++i) {
            if (!faces->hasChildren(i)) {
                faces->divide(i);
            }
        }
    }
    
    if (lagrangian) {
        if (typeid(seed) == typeid(IcosTriSphereSeed) || typeid(seed) == typeid(CubedSphereSeed)) {
            lagCoords = std::shared_ptr<Coords>(new SphericalCoords(*(dynamic_cast<SphericalCoords*>(coords.get()))));    
        }
        else if (typeid(seed) == typeid(TriHexSeed) || typeid(seed) == typeid(QuadRectSeed)) {
            lagCoords = std::shared_ptr<Coords>(new EuclideanCoords(*(dynamic_cast<EuclideanCoords*>(coords.get()))));
        }
        lagCoords->setLogProc(prank);
        
        edges->makeLagrangian(lagCoords);
        faces->makeLagrangian(lagCoords);
    }
    
    ss.str(std::string());
    ss << "** Mesh created **" <<std::endl;
    ss << "\tfaces->nMax() = " << faces->nMax() << ", faces->n() = " << faces->n() << ", faces->nLeaves() = " << faces->nLeaves() << std::endl;
    ss << "\tedges->nMax() = " << edges->nMax() << ", edges->n() = " << edges->n() << ", edges->nLeaves() = " << edges->nLeaves() << std::endl;    
    ss << "\tsurface area = " << faces->surfaceArea();
    statusMsg.resetMsgString(ss.str());
    log->logMessage(statusMsg);
    log->endSection();
}

index_type PolyMesh2d::locateFaceContainingPoint(const XyzVector& queryPt) const {
    const index_type rootFace = nearestRootFace(queryPt);
    index_type treeResult = treeSearch(queryPt, rootFace);
    return walkSearch(queryPt, treeResult);
}

index_type PolyMesh2d::walkSearch(const XyzVector& queryPt, index_type startIndex) const {
    if (faces->hasChildren(startIndex)) {
        OutputMessage errMsg("walk search requires a leaf face for input", OutputMessage::errorPriority, 
            "PolyMesh2d::walkSearch");
        log->logMessage(errMsg);
        return -1;
    }
    else {
        scalar_type dist = coords->distance(queryPt, faces->centroid(startIndex));
        index_type nearestFace = startIndex;
        const std::vector<index_type> adjFaces = faces->ccwAdjacentFaces(startIndex);
        for (index_type i = 0; i < adjFaces.size(); ++i) {
            if (adjFaces[i] >= 0) {
                const scalar_type test_dist = coords->distance(queryPt, faces->centroid(adjFaces[i]));
                if (test_dist < dist) {
                    dist = test_dist;
                    nearestFace = adjFaces[i];
                }
            }
        }
        if (nearestFace == startIndex) 
            return startIndex;
        else
            return walkSearch(queryPt, nearestFace);
    }
}

index_type PolyMesh2d::treeSearch(const XyzVector& queryPt, index_type startIndex) const {
    index_type result = -1;
    scalar_type dist = std::numeric_limits<scalar_type>::max();
    index_type nearestChild = -1;
    if (faces->hasChildren(startIndex)) {
        const std::vector<index_type> childFaces = faces->children(startIndex);
        for (int i = 0; i < 4; ++i) {
            const scalar_type test_dist = coords->distance(queryPt, faces->centroid(childFaces[i]));
            if (test_dist < dist) {
                dist = test_dist;
                nearestChild = childFaces[i];
            }
        }
        return treeSearch(queryPt, nearestChild);
    }
    else {
        return startIndex;
    }
}

index_type PolyMesh2d::nearestRootFace(const XyzVector& queryPt) const {
    index_type result = -1;
    scalar_type dist = std::numeric_limits<scalar_type>::max();
    for (int i = 0; i < nRootFaces; ++i) {
        const scalar_type test_dist = coords->distance(queryPt, faces->centroid(i));
        if ( test_dist < dist ) {
            dist = test_dist;
            result = i;
        }
    }
    return result;
}

std::shared_ptr<Field> PolyMesh2d::createFaceAreaField() const {
    std::shared_ptr<Field> result(new Field(faces->n(), 1, "area", "length^2"));
    for (index_type i = 0; i < faces->n(); ++i) {
        result->insert(faces->area(i));
    }
    return result;
}

void PolyMesh2d::writeToVTKFile(const std::string& fname, const std::string& desc) const {
    VtkWriter writer;
    std::ofstream fs(fname);
    if (!fs.is_open()) {
        OutputMessage errMsg("cannot open .vtk file", OutputMessage::errorPriority, "PolyMesh2d::writeToVTKFile");
        log->logMessage(errMsg);
        throw std::ios_base::failure("file write error");
    }
    writer.writeVTKHeader(fs, desc);
    writer.writeCoordsToVTKPoints(fs, coords);
    writer.writeFacesToVTKPolygons(fs, faces);
    writer.writeVTKPointDataHeader(fs, coords->n());
    if (lagrangian) {
        writer.writeCoordsToVTKPointData(fs, lagCoords, lagrangian);
    }
    writer.writeVTKCellDataHeader(fs, faces->nLeaves());
    writer.writeFaceAreaToVTKCellData(fs, faces);
}

}
