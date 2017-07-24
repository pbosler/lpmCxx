#ifndef _LPM_PARTICLES_H_
#define _LPM_PARTICLES_H_

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmCoords.h"
#include "LpmField.h"
#include "LpmXyzVector.h"
#include "LpmPolyMesh2d.h"
#include <memory>
#include <map>
#include <vector>
#include <string>

namespace Lpm {

class Particles {
    public:
        Particles(const std::shared_ptr<Coords> crds, const std::shared_ptr<Coords> lagCrds = NULL, const int prank = 0);
        Particles(const std::shared_ptr<Coords> crds, const std::shared_ptr<Coords> lagCrds = NULL, 
            const std::vector<std::shared_ptr<Field>> fields = std::vector<std::shared_ptr<Field>>(), const int prank = 0);
        Particles(const index_type nMax, const std::vector<std::string>& fnames, const std::vector<int>& fdims, 
            const std::vector<std::string>& funits, const GeometryType gkind, const bool lagrangian = false, 
            const scalar_type domainRadius = 1.0, const int prank = 0);
        
        Particles(MeshSeed& seed, const int maxRecursionLevel, const scalar_type domainRadius = 1.0, 
            const bool lagrangian = false, const int prank = 0);

        std::shared_ptr<Field> getFieldPtr(const std::string& fieldname);
        
        std::shared_ptr<Coords> getCoordPtr(const bool lagrangian = false);
        
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
        
        inline void setLogProc(const int rank) {log->setProcRank(rank);}
        
        inline index_type n() const {return _coords->n();}
        inline index_type nMax() const {return _coords->nMax();}
        
        void createField(const std::string& fieldName, const std::string& fieldUnits, const int fieldDim);
        void registerField(const std::shared_ptr<Field>& field);
        void initializeFieldWithFunction(const std::string& fieldname, const AnalyticFunction* fn);
    
        inline GeometryType geometry() const {return _coords->geometry();}
        inline std::vector<BoundaryCondition> boundaryConditions() const {return _boundaryConditions;}
        
    protected:
        Particles(){};
        std::shared_ptr<Coords> _coords;
        std::shared_ptr<Coords> _lagCoords;
        
        std::vector<BoundaryCondition> _boundaryConditions;
        
        std::map<std::string, std::shared_ptr<Field>> _fieldMap;
        
        static std::unique_ptr<Logger> log;
};

}

#endif
