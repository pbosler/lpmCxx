#include "LpmPolyMesh2d.h"
#include "LpmEuclideanCoords.h"
#include "LpmSphericalCoords.h"
#include "LpmOutputMessage.h"
#include "LpmTriFaces.h"
#include "LpmQuadFaces.h"
#include "LpmParticles.h"
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
    log->setProcRank(prank);
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
        _coords = std::shared_ptr<Coords>(new SphericalCoords(nMaxVerts, domainRadius));
    }
    else if (typeid(seed) == typeid(TriHexSeed) || typeid(seed) == typeid(QuadRectSeed)) {
        _coords = std::shared_ptr<Coords>(new EuclideanCoords(nMaxVerts, PLANAR_GEOMETRY));
    }
    _coords->setLogProc(prank);
    
    _edges = std::shared_ptr<Edges>(new Edges(nMaxEdges, _coords));
    _edges->setLogProc(prank);
    
    if (typeid(seed) == typeid(IcosTriSphereSeed) || typeid(seed) == typeid(TriHexSeed)) {
        _faces = std::shared_ptr<Faces>(new TriFaces(nMaxFaces, _edges, _coords));
    }
    else if  (typeid(seed) == typeid(QuadRectSeed) || typeid(seed) == typeid(CubedSphereSeed)) {
        _faces = std::shared_ptr<Faces>(new QuadFaces(nMaxFaces, _edges, _coords));
    }
    _faces->setLogProc(prank);
    
    seed.initMeshFromSeed(_coords, _edges, _faces, domainRadius);
//     ss.str(std::string());
// //    ss << "seed.infoString() = " << seed.infoString();
//     statusMsg.resetMsgString(ss.str());
//     log->logMessage(statusMsg);
    
    for (int k = 0; k < maxRecursionLevel; ++k) {
        const index_type n_faces = _faces->n();
        for (index_type i = 0; i < n_faces; ++i) {
            if (!_faces->hasChildren(i)) {
                _faces->divide(i);
            }
        }
    }
    
    if (lagrangian) {
        if (typeid(seed) == typeid(IcosTriSphereSeed) || typeid(seed) == typeid(CubedSphereSeed)) {
            _lagCoords = std::shared_ptr<Coords>(new SphericalCoords(*(dynamic_cast<SphericalCoords*>(_coords.get()))));    
        }
        else if (typeid(seed) == typeid(TriHexSeed) || typeid(seed) == typeid(QuadRectSeed)) {
            _lagCoords = std::shared_ptr<Coords>(new EuclideanCoords(*(dynamic_cast<EuclideanCoords*>(_coords.get()))));
        }
        _lagCoords->setLogProc(prank);
        
        _edges->makeLagrangian(_lagCoords);
        _faces->makeLagrangian(_lagCoords);
    }
    
    ss.str(std::string());
    ss << "** Mesh created **" <<std::endl;
    ss << "\t_faces->nMax() = " << _faces->nMax() << ", _faces->n() = " << _faces->n() << ", _faces->nLeaves() = " << _faces->nLeaves() << std::endl;
    ss << "\t_edges->nMax() = " << _edges->nMax() << ", _edges->n() = " << _edges->n() << ", _edges->nLeaves() = " << _edges->nLeaves() << std::endl;    
    ss << "\tsurface area = " << _faces->surfaceArea();
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
    if (_faces->hasChildren(startIndex)) {
        OutputMessage errMsg("walk search requires a leaf face for input", OutputMessage::errorPriority, 
            "PolyMesh2d::walkSearch");
        log->logMessage(errMsg);
        return -1;
    }
    else {
        scalar_type dist = _coords->distance(queryPt, _faces->centroid(startIndex));
        index_type nearestFace = startIndex;
        const std::vector<index_type> adjFaces = _faces->ccwAdjacentFaces(startIndex);
        for (index_type i = 0; i < adjFaces.size(); ++i) {
            if (adjFaces[i] >= 0) {
                const scalar_type test_dist = _coords->distance(queryPt, _faces->centroid(adjFaces[i]));
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
    if (_faces->hasChildren(startIndex)) {
        const std::vector<index_type> childFaces = _faces->children(startIndex);
        for (int i = 0; i < 4; ++i) {
            const scalar_type test_dist = _coords->distance(queryPt, _faces->centroid(childFaces[i]));
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
        const scalar_type test_dist = _coords->distance(queryPt, _faces->centroid(i));
        if ( test_dist < dist ) {
            dist = test_dist;
            result = i;
        }
    }
    return result;
}

std::shared_ptr<Field> PolyMesh2d::createFaceAreaField() const {
    std::shared_ptr<Field> result(new Field(_faces->n(), 1, "area", "length^2"));
    for (index_type i = 0; i < _faces->n(); ++i) {
        result->insert(_faces->area(i));
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
    writer.writeCoordsToVTKPoints(fs, _coords);
    writer.writeFacesToVTKPolygons(fs, _faces);
    writer.writeVTKPointDataHeader(fs, _coords->n());
    if (lagrangian) {
        writer.writeCoordsToVTKPointData(fs, _lagCoords, lagrangian);
    }
    writer.writeVTKCellDataHeader(fs, _faces->nLeaves());
    writer.writeFaceAreaToVTKCellData(fs, _faces);
}

}
