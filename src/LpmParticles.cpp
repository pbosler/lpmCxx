#include "LpmParticles.h"
#include "LpmEuclideanCoords.h"
#include "LpmSphericalCoords.h"

namespace Lpm {

Particles::Particles(const std::shared_ptr<Coords> crds, const std::shared_ptr<Coords> lagCrds) : 
_coords(crds), _lagCoords(lagCrds) {};

Particles::Particles(const std::shared_ptr<Coords> crds, const std::shared_ptr<Coords> lagCrds, 
    const std::vector<std::shared_ptr<Field>> fields) : _coords(crds), _lagCoords(lagCrds) {
    for (index_type i = 0; i < fields.size(); ++i) {
        _fieldMap.emplace(fields[i]->name(), fields[i]);
    }    
}

Particles::Particles(const index_type nMax, const std::vector<std::string>& fnames, const std::vector<int>& fdims, 
    const std::vector<std::string>& funits, const GeometryType gkind, const bool lagrangian, 
    const scalar_type domainRadius) {
    
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
}

Particles::Particles(const std::shared_ptr<PolyMesh2d> mesh) {
    _coords = mesh->getPhysCoords();
    _lagCoords = mesh->getLagCoords();
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
