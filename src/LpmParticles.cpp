#include "LpmParticles.h"
#include "LpmEuclideanCoords.h"
#include "LpmSphericalCoords.h"

namespace Lpm {

Particles::Particles(const std::shared_ptr<Coords> crds, const bool lagrangian) : _coords(crds), 
    (lagrangian ? _lagCoords(crds) : _lagCoords(NULL)) {
};

Particles::Particles(const std::shared_ptr<Coords> crds, const bool lagrangian, 
    const std::vector<std::shared_ptr<Field>> fields) : _coords(crds,
    (lagrangian ? _lagCoords(crds) : _lagCoords(NULL)) {
    for (index_type i = 0; i < fields.size(); ++i) {
        _fieldMap.emplace(fields[i]->name(), fields[i]);
    }    
}

Particles::Particles(const index_type nMax, const std::vector<std::string>& fnames, const std::vector<int>& fdims, 
    const std::vector<std::string>& funits, const GeometryType gkind, const bool lagrangian = false, 
    const scalar_type domainRadius = 1.0) {
    
    if (gkind == PLANAR_GEOMETRY) {
        _coords = std::shared_ptr<Coords>(new EuclideanCoords(nMax));
        if (lagrangian)
            _lagCoords = std::shared_ptr<Coords>(new EuclideanCoords(nMax));
    }
    else if (gkind == CARTESIAN_3D_GEOMETRY) {
        _coords = std::shared_ptr<Coords>(new EuclideanCoords(nMax, true));
        if (lagrangian)
            _lagCoords = std::shared_ptr<Coords>(new EuclideanCoords(nMax, true));
    }
    else if (gkind == SPHERICAL_SURFACE_GEOMETRY) {
        _coords = std::shared_ptr<Coords>(new SphericalCoords(nMax, domainRadius));
        if (lagrangian) 
            _lagCoords = std::shared_ptr<Coords>(new SphericalCoords(nMax, domainRadius));
    }
    
    for (index_type i = 0; i < fnames.size(); ++i) {
        _fieldMap.emplace(fnames[i], std::shared_ptr<Field>(new Field(nMax, fdims[i], fnames[i], funits[i]));
    }
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

}
