#ifndef _LPM_DIRECT_SUM_H_
#define _LPM_DIRECT_SUM_H_

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmXyzVector.h"
#include "LpmAnalyticFunctions.h"
#include "LpmLogger.h"
#include "LpmParticles.h"
#include "LpmMeshedParticles.h"
#include "LpmMPIReplicatedData.h"
#include "LpmScalarKernel.h"
#include <memory>

namespace Lpm {

class DirectSum {
    public: 
        DirectSum(std::shared_ptr<MeshedParticles> pmesh, const std::shared_ptr<ScalarKernel> kernel, 
            const std::string tgtFieldName, const std::string srcFieldName);
        DirectSum(std::shared_ptr<Particles> particles, const std::shared_ptr<ScalarKernel> kernel, 
            const std::string tgtFieldName, const std::string srcFieldName, const std::string srcWeightName = std::string());
            
        void meshSolve(const MPIReplicatedData& mpiVerts, const MPIReplicatedData& mpiFaces) const;
        
        void meshfreeSolve(const MPIReplicatedData& mpiParticles) const;
        
        void meshBroadcast(const MPIReplicatedData& mpiVerts, const MPIReplicatedData& mpiFaces) const;
        void meshfreeBroadcast(const MPIReplicatedData& mpiParticles) const;
    
    protected:
        std::weak_ptr<ScalarKernel> _kernel;
        
        bool useWeights;
        std::weak_ptr<Field> _srcWeights;
    
        // mesh solve variables
        std::weak_ptr<PolyMesh2d> _mesh;
        std::weak_ptr<Field> _vertTgt;
        std::weak_ptr<Field> _faceSrc;
        std::weak_ptr<Field> _faceTgt;
        
        // meshfree solve variables
        std::weak_ptr<Coords> _srcLocs;
        std::weak_ptr<Field> _srcVals;
        std::weak_ptr<Field> _tgtVals;
};

}

#endif
