#include "LpmMeshedParticles.h"
#include "LpmOutputMessage.h"
#include "LpmVtkFileIO.h"
#include <fstream>
#include <exception>
#include <sstream>

namespace Lpm {

std::unique_ptr<Logger> MeshedParticles::log(new Logger(OutputMessage::debugPriority, "MeshedParticles_log"));

MeshedParticles::MeshedParticles(MeshSeed& seed, const int maxRecursionLevel, const bool isLagrangian, 
    const scalar_type domainRadius, const int prank) {
    mesh2d = std::shared_ptr<PolyMesh2d>(new PolyMesh2d(seed, maxRecursionLevel, isLagrangian, domainRadius, prank));
    _coords = mesh2d->getPhysCoords();
    _lagCoords = mesh2d->getLagCoords();
    std::shared_ptr<Field> faceAreaPtr = mesh2d->createFaceAreaField();
    faceAreaPtr->rename("faceArea");
    _faceFieldMap.emplace("faceArea", faceAreaPtr);
    log->setProcRank(prank);
};

void MeshedParticles::createVertexField(const std::string& fieldName, const std::string& fieldUnits, const int fieldDim) {
    _vertexFieldMap.emplace(fieldName, std::shared_ptr<Field>(new Field(mesh2d->nVertices(), fieldDim, fieldName, fieldUnits)));
    _vertexFieldMap.at(fieldName)->initializeToConstant(mesh2d->nVertices(), 0.0);
}

void MeshedParticles::createEdgeField(const std::string& fieldName, const std::string& fieldUnits, const int fieldDim) {
    _edgeFieldMap.emplace(fieldName, std::shared_ptr<Field>(new Field(mesh2d->nEdges(), fieldDim, fieldName, fieldUnits)));
    _edgeFieldMap.at(fieldName)->initializeToConstant(mesh2d->nEdges(), 0.0);
}

void MeshedParticles::createFaceField(const std::string& fieldName, const std::string& fieldUnits, const int fieldDim) {
    _faceFieldMap.emplace(fieldName, std::shared_ptr<Field>(new Field(mesh2d->nFaces(), fieldDim, fieldName, fieldUnits)));
    _faceFieldMap.at(fieldName)->initializeToConstant(mesh2d->nFaces(), 0.0);
}

void MeshedParticles::initializeVertexFieldWithFunction(const std::string& fieldName, const AnalyticFunction* fn){
    std::shared_ptr<Field> fptr = getVertexFieldPtr(fieldName);
    fptr->clear();
    if (fptr->nDim() == 1) {
        fptr->initializeToScalarFunction(_coords.get(), fn);
    }
    else {
        fptr->initializeToVectorFunction(_coords.get(), fn);
    }
}

void MeshedParticles::initializeEdgeFieldWithFunction(const std::string& fieldName, const AnalyticFunction* fn) {
    std::shared_ptr<Field> fptr = getEdgeFieldPtr(fieldName);
    fptr->clear();
    if (fptr->nDim() == 1) {
        for (index_type i = 0; i < mesh2d->nEdges(); ++i) {
            fptr->insert((mesh2d->_edges->hasChildren(i) ? 0.0 : fn->evaluateScalar(mesh2d->_edges->midpoint(i))));
        }
    }
    else {
        for (index_type i = 0; i < mesh2d->nEdges(); ++i) {
           fptr->insert((mesh2d->_edges->hasChildren(i) ? XyzVector(0.0, 0.0, 0.0) : fn->evaluateVector(mesh2d->_edges->midpoint(i))));
        }
    }
}

void MeshedParticles::initializeFaceFieldWithFunction(const std::string& fieldName, const AnalyticFunction* fn) {
    std::shared_ptr<Field> fptr = getFaceFieldPtr(fieldName);
    fptr->clear();
    if (fptr->nDim() == 1) {
        for (index_type i = 0; i < mesh2d->nFaces(); ++i) {
            fptr->insert((mesh2d->_faces->hasChildren(i) ? 0.0 : fn->evaluateScalar(mesh2d->_faces->centroid(i))));
        }
    }
    else {
        for (index_type i = 0; i < mesh2d->nFaces(); ++i) {
            fptr->insert((mesh2d->_faces->hasChildren(i) ? XyzVector(0.0,0.0,0.0) : fn->evaluateVector(mesh2d->_faces->centroid(i))));
        }
    }
}


std::shared_ptr<Field> MeshedParticles::getVertexFieldPtr(const std::string& fieldName){
    std::shared_ptr<Field> result;
    if (!fieldName.empty()) {
        try {
            result = _vertexFieldMap.at(fieldName);
        }
        catch (std::out_of_range& oor) {
            std::stringstream ss;
            ss << "vertex field '" << fieldName << "' not registered, caught std::out_of_range() exception.";
        
            OutputMessage warnMsg(ss.str(), OutputMessage::warningPriority, "MeshedParticles::getVertexFieldPtr");
            log->logMessage(warnMsg);
        }
    }
    return result;
}

std::shared_ptr<Field> MeshedParticles::getEdgeFieldPtr(const std::string& fieldName){
    std::shared_ptr<Field> result;
    if (!fieldName.empty()) {
        try {
            result = _edgeFieldMap.at(fieldName);
        }
        catch (std::out_of_range& oor) {
            std::stringstream ss;
            ss << "edge field '" << fieldName << "' not registered, caught std::out_of_range() exception.";
        
            OutputMessage warnMsg(ss.str(), OutputMessage::warningPriority, "MeshedParticles::getEdgeFieldPtr");
            log->logMessage(warnMsg);
        }
    }
    return result;
}

std::shared_ptr<Field> MeshedParticles::getFaceFieldPtr(const std::string& fieldName){
    std::shared_ptr<Field> result;
    if (!fieldName.empty()) {
        try {
            result = _faceFieldMap.at(fieldName);
        }
        catch (std::out_of_range& oor) {
            std::stringstream ss;
            ss << "face field '" << fieldName << "' not registered, caught std::out_of_range() exception.";
        
            OutputMessage warnMsg(ss.str(), OutputMessage::warningPriority, "MeshedParticles::getFaceFieldPtr");
            log->logMessage(warnMsg);
        }
    }
    return result;
}

void MeshedParticles::writeToVtkFile(const std::string& fname, const std::string& desc) const {
    VtkWriter writer;
    std::ofstream fs(fname);
    if (!fs.is_open()) {
        OutputMessage errMsg("cannot open .vtk file", OutputMessage::errorPriority, "MeshedParticles::writeToVtkFile");
        log->logMessage(errMsg);
        throw std::ios_base::failure("file write error");
    }
    writer.writeVTKHeader(fs, desc);
    writer.writeCoordsToVTKPoints(fs, _coords);
    writer.writeFacesToVTKPolygons(fs, mesh2d->_faces);
    writer.writeVTKPointDataHeader(fs, _coords->n());
    if (mesh2d->isLagrangian()) {
        writer.writeCoordsToVTKPointData(fs, _lagCoords, true);
    }
    for (auto& elem : _vertexFieldMap) {
        writer.writeFieldToVTKData(fs, elem.second);
    }
    writer.writeVTKCellDataHeader(fs, mesh2d->nLeafFaces());
    for (auto& elem : _faceFieldMap) {
        writer.writeFaceFieldToVTKCellData(fs, mesh2d->_faces, elem.second);
    }
}

}
