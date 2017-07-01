#ifndef _LPM_PARTICLES_H_
#define _LPM_PARTICLES_H_

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmCoords.h"
#include "LpmFields.h"
#include "LpmXyzVector.h"
#include <memory>
#include <map>
#include <vector>
#include <string>

namespace Lpm {

class Particles {
    public:
        Particles(const std::shared_ptr<Coords> crds, const bool lagrangian = false);
        Particles(const std::shared_ptr<Coords> crds, const std::vector<std::shared_ptr<Field>> fields, const bool lagrangian = false);
        Particles(const index_type nMax, const std::vector<std::string>& fnames, const std::vector<int>& fdims, 
            const std::vector<std::string>& funits, const GeometryType gkind);
    
        void insert(const XyzVector& newCoord);
        void replaceCoordinate(const index_type ind, const XyzVector& newCoord, const bool lagrangian = false);
        void setScalarFieldValues(const index_type pIndex, const std::vector<std::string>& fnames, const std::vector<scalar_type>& vals);
        void setVectorFieldValues(const index_type pIndex, const std::vector<std::string>& fnames, const std::vector<XyzVector>& vecVals);
        
        inline scalar_type scalarFieldValue(const std::string& fname, const index_type pIndex) const {
            return _fieldMap.at(fname)->getScalar(pIndex);}
        inline XyzVector vector2dFieldValue(const std::string& fname, const index_type pIndex) const {
            return _fieldMap.at(fname)->get2dVector(pIndex);
        }
        inline XyzVector vector3dFieldValue(const std::string& fname, const index_type pIndex) const {
            return _fieldMap.at(fname)->get3dVector(pIndex);
        }
        
        inline index_type n() const {return _coords->n();}
        inline index_type nMax() const {return _coords->nMax();}
        
        void createField(const std::string fieldName, const std::string fieldUnits, const int fieldDim);
        void registerField(const std::shared_ptr<Field>);
    
        
    protected:
        
    
        std::shared_ptr<Coords> _coords;
        std::shared_ptr<Coords> _lagCoords;
        
        std::map<std::string, std::shared_ptr<Field>> _fieldMap;
};

}

#endif
