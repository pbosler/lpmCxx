#ifndef _LPM_MESHED_PARTICLES_H_
#define _LPM_MESHED_PARTICLES_H_

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmParticles.h"
#include "LpmPolyMesh2d.h"

namespace Lpm {

class MeshedParticles : public Particles {
    public:
        typedef std::map<std::string, std::shared_ptr<Field>> field_map_type;
        
        inline index_type nFaces() const {return mesh2d->nFaces();}
        inline index_type nLeafFaces() const {return mesh2d->nLeafFaces();}
        
        inline index_type nVertices() const {return mesh2d->nVertices();}
        inline index_type nEdges() const {return mesh2d->nEdges();}
        inline index_type nLeafEdges() const {return mesh2d->nLeafEdges();}
        
        MeshedParticles(MeshSeed& seed, const int maxRecursionLevel, const bool isLagrangian = false, 
            const scalar_type domainRadius = 1.0, const int prank = 0);
    
        void createVertexField(const std::string& fieldName, const std::string& fieldUnits, const int fieldDim);
        void createEdgeField(const std::string& fieldName, const std::string& fieldUnits, const int fieldDim); 
        void createFaceField(const std::string& fieldName, const std::string& fieldUnits, const int fieldDim);
        
        std::shared_ptr<Field> getVertexFieldPtr(const std::string& fieldName);
        std::shared_ptr<Field> getEdgeFieldPtr(const std::string& fieldName);
        std::shared_ptr<Field> getFaceFieldPtr(const std::string& fieldName);
        
        void initializeVertexFieldWithFunction(const std::string& fieldName, const AnalyticFunction* fn);
        void initializeEdgeFieldWithFunction(const std::string& fieldName, const AnalyticFunction* fn);
        void initializeFaceFieldWithFunction(const std::string& fieldName, const AnalyticFunction* fn);
        
        void writeToVtkFile(const std::string& fname, const std::string& desc = "") const;
        
        inline std::shared_ptr<PolyMesh2d> meshPtr() {return mesh2d;}
        
        inline void setLogProc(const int rank) {log->setProcRank(rank);}
        
    protected:
        std::shared_ptr<PolyMesh2d> mesh2d;
        
        field_map_type _vertexFieldMap;
        field_map_type _edgeFieldMap;
        field_map_type _faceFieldMap;
    
        static std::unique_ptr<Logger> log;
};

}

#endif
