#include "LpmParticles.h"
#include "LpmOutputMessage.h"
#include "LpmEuclideanCoords.h"
#include "LpmSphericalCoords.h"
#include <sstream>
#include <exception>

namespace Lpm {

std::unique_ptr<Logger> Particles::log(new Logger(OutputMessage::debugPriority, "Particles_log"));

Particles::Particles(const std::shared_ptr<Coords> crds, const std::shared_ptr<Coords> lagCrds, const int prank) : 
_coords(crds), _lagCoords(lagCrds) {log->setProcRank(prank);};

Particles::Particles(const std::shared_ptr<Coords> crds, const std::shared_ptr<Coords> lagCrds, 
    const std::vector<std::shared_ptr<Field>> fields, const int prank) : _coords(crds), _lagCoords(lagCrds) {
    for (index_type i = 0; i < fields.size(); ++i) {
        _fieldMap.emplace(fields[i]->name(), fields[i]);
    }    
    log->setProcRank(prank);
}

Particles::Particles(const index_type nMax, const std::vector<std::string>& fnames, const std::vector<int>& fdims, 
    const std::vector<std::string>& funits, const GeometryType gkind, const bool lagrangian, 
    const scalar_type domainRadius, const int prank) {
    
    if (gkind == PLANAR_GEOMETRY || gkind == CARTESIAN_3D_GEOMETRY) {
        _coords = std::shared_ptr<Coords>(new EuclideanCoords(nMax, gkind));
        if (lagrangian)
            _lagCoords = std::shared_ptr<Coords>(new EuclideanCoords(nMax, gkind));
    }
    else if (gkind == SPHERICAL_SURFACE_GEOMETRY) {
        _coords = std::shared_ptr<Coords>(new SphericalCoords(nMax, domainRadius));
        if (lagrangian) 
            _lagCoords = std::shared_ptr<Coords>(new SphericalCoords(nMax, domainRadius));
    }
    
    for (index_type i = 0; i < fnames.size(); ++i) {
        _fieldMap.emplace(fnames[i], std::shared_ptr<Field>(new Field(nMax, fdims[i], fnames[i], funits[i])));
    }
    log->setProcRank(prank);
}

std::shared_ptr<Field> Particles::getFieldPtr(const std::string& fieldname) {
    std::shared_ptr<Field> result;
    try {
        result = _fieldMap.at(fieldname);
    }
    catch ( std::out_of_range& oor) {
        std::stringstream ss;
        ss << "field '" << fieldname << "' not registered.";
        
        OutputMessage warnMsg(ss.str(), OutputMessage::warningPriority, "Particles::getFieldPtr");
        log->logMessage(warnMsg);
    }
    return result;
}

void Particles::insert(const XyzVector& newCoord) {
    _coords->insert(newCoord);
    if (_lagCoords)
        _lagCoords->insert(newCoord);
}

void Particles::replaceCoordinate(const index_type ind, const XyzVector& newCoord, const bool lagrangian) {
    if (lagrangian) {
        _coords->replace(ind, newCoord);
    }
    else {
        _lagCoords->replace(ind, newCoord);
    }
}

void Particles::createField(const std::string& fieldName, const std::string& fieldUnits, const int fieldDim) {
    _fieldMap.emplace(fieldName, std::shared_ptr<Field>(new Field(_coords->nMax(), fieldDim, fieldName, fieldUnits)));
}

void Particles::registerField(const std::shared_ptr<Field>& field){
    _fieldMap.emplace(field->name(), field);
}

void Particles::initializeFieldWithFunction(const std::string& fieldName, const AnalyticFunction* fn) {
    std::shared_ptr<Field> field_ptr = _fieldMap.at(fieldName);
    if (field_ptr->nDim() == 1) {
        field_ptr->initializeToScalarFunction(_coords.get(), fn);
    }
    else {
        field_ptr->initializeToVectorFunction(_coords.get(), fn);
    }
}

}
