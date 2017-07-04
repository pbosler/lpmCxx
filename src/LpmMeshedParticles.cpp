#include "LpmMeshedParticles.h"
#include "LpmOutputMessage.h"
#include "LpmVtkFileIO.h"
#include <fstream>
#include <exception>
#include <sstream>

namespace Lpm {

MeshedParticles::MeshedParticles(MeshSeed& seed, const int maxRecursionLevel, const bool isLagrangian, const scalar_type domainRadius) {
    mesh2d = std::shared_ptr<PolyMesh2d>(new PolyMesh2d(seed, maxRecursionLevel, isLagrangian, domainRadius));
    _coords = mesh2d->getPhysCoords();
    _lagCoords = mesh2d->getLagCoords();
    std::shared_ptr<Field> faceAreaPtr = mesh2d->createFaceAreaField();
    faceAreaPtr->rename("faceArea");
    _faceFieldMap.emplace("faceArea", faceAreaPtr);
};

void MeshedParticles::createVertexField(const std::string& fieldName, const std::string& fieldUnits, const int fieldDim) {
    _vertexFieldMap.emplace(fieldName, std::shared_ptr<Field>(new Field(mesh2d->coords->n(), fieldDim, fieldName, fieldUnits)));
    _vertexFieldMap.at(fieldName)->initializeToConstant(mesh2d->coords->n(), 0.0);
}

void MeshedParticles::createEdgeField(const std::string& fieldName, const std::string& fieldUnits, const int fieldDim) {
    _edgeFieldMap.emplace(fieldName, std::shared_ptr<Field>(new Field(mesh2d->edges->n(), fieldDim, fieldName, fieldUnits)));
    _edgeFieldMap.at(fieldName)->initializeToConstant(mesh2d->edges->n(), 0.0);
}

void MeshedParticles::createFaceField(const std::string& fieldName, const std::string& fieldUnits, const int fieldDim) {
    _faceFieldMap.emplace(fieldName, std::shared_ptr<Field>(new Field(mesh2d->faces->n(), fieldDim, fieldName, fieldUnits)));
    _faceFieldMap.at(fieldName)->initializeToConstant(mesh2d->faces->n(), 0.0);
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
        for (index_type i = 0; i < mesh2d->edges->n(); ++i) {
            fptr->insert((mesh2d->edges->hasChildren(i) ? 0.0 : fn->evaluateScalar(mesh2d->edges->midpoint(i))));
        }
    }
    else {
        for (index_type i = 0; i < mesh2d->edges->n(); ++i) {
           fptr->insert((mesh2d->edges->hasChildren(i) ? XyzVector(0.0, 0.0, 0.0) : fn->evaluateVector(mesh2d->edges->midpoint(i))));
        }
    }
}

void MeshedParticles::initializeFaceFieldWithFunction(const std::string& fieldName, const AnalyticFunction* fn) {
    std::shared_ptr<Field> fptr = getFaceFieldPtr(fieldName);
    fptr->clear();
    if (fptr->nDim() == 1) {
        for (index_type i = 0; i < mesh2d->faces->n(); ++i) {
            fptr->insert((mesh2d->faces->hasChildren(i) ? 0.0 : fn->evaluateScalar(mesh2d->faces->centroid(i))));
        }
    }
    else {
        for (index_type i = 0; i < mesh2d->faces->n(); ++i) {
            fptr->insert((mesh2d->faces->hasChildren(i) ? XyzVector(0.0,0.0,0.0) : fn->evaluateVector(mesh2d->faces->centroid(i))));
        }
    }
}


std::shared_ptr<Field> MeshedParticles::getVertexFieldPtr(const std::string& fieldName){
    std::shared_ptr<Field> result;
    try {
        result = _vertexFieldMap.at(fieldName);
    }
    catch (std::out_of_range& oor) {
        std::stringstream ss;
        ss << "vertex field '" << fieldName << "' not registered, caught std::out_of_range() exception.";
        
        OutputMessage warnMsg(ss.str(), OutputMessage::warningPriority, "MeshedParticles::getVertexFieldPtr");
        log->logMessage(warnMsg);
    }
    return result;
}

std::shared_ptr<Field> MeshedParticles::getEdgeFieldPtr(const std::string& fieldName){
    std::shared_ptr<Field> result;
    try {
        result = _edgeFieldMap.at(fieldName);
    }
    catch (std::out_of_range& oor) {
        std::stringstream ss;
        ss << "edge field '" << fieldName << "' not registered, caught std::out_of_range() exception.";
        
        OutputMessage warnMsg(ss.str(), OutputMessage::warningPriority, "MeshedParticles::getEdgeFieldPtr");
        log->logMessage(warnMsg);
    }
    return result;
}

std::shared_ptr<Field> MeshedParticles::getFaceFieldPtr(const std::string& fieldName){
    std::shared_ptr<Field> result;
    try {
        result = _faceFieldMap.at(fieldName);
    }
    catch (std::out_of_range& oor) {
        std::stringstream ss;
        ss << "face field '" << fieldName << "' not registered, caught std::out_of_range() exception.";
        
        OutputMessage warnMsg(ss.str(), OutputMessage::warningPriority, "MeshedParticles::getFaceFieldPtr");
        log->logMessage(warnMsg);
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
    writer.writeFacesToVTKPolygons(fs, mesh2d->faces);
    writer.writeVTKPointDataHeader(fs, _coords->n());
    if (mesh2d->isLagrangian()) {
        writer.writeCoordsToVTKPointData(fs, _lagCoords, true);
    }
    for (auto& elem : _vertexFieldMap) {
        writer.writeFieldToVTKData(fs, elem.second);
    }
    writer.writeVTKCellDataHeader(fs, mesh2d->nLeafFaces());
    for (auto& elem : _faceFieldMap) {
        writer.writeFaceFieldToVTKCellData(fs, mesh2d->faces, elem.second);
    }
}

}
