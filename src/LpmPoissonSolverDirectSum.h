#ifndef _LPM_POISSON_SOLVER_DIRECT_SUM_H_
#define _LPM_POISSON_SOLVER_DIRECT_SUM_H_

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmCoords.h"
#include "LpmFaces.h"
#include "LpmMeshedParticles.h"
#include "LpmMPIReplicatedData.h"
#include "LpmPolyMesh2d.h"
#include "LpmScalarKernel.h"
#include <memory>


namespace Lpm {

class PoissonSolverDirectSum {
    public:
        PoissonSolverDirectSum(const std::shared_ptr<PolyMesh2d>& mesh, const std::shared_ptr<Field>& src, 
            const std::shared_ptr<Field>& vertexPotential, const std::shared_ptr<Field>& facePotential, 
            const ScalarKernel* kernel) : 
            _faces(mesh->getFaces()), _verts(mesh->getPhysCoords()), _src(src), _vertpot(vertexPotential), 
            _facepot(facePotential), _kernel(kernel) {};

        void solve(const MPIReplicatedData& mpiVerts, const MPIReplicatedData& mpiFaces) const;
        
        void broadcastSolution(const MPIReplicatedData& mpiVerts, const MPIReplicatedData& mpiFaces) const;
            
    protected:
        std::shared_ptr<Coords> _verts;
        std::shared_ptr<Faces> _faces;
        std::shared_ptr<Field> _src;
        std::shared_ptr<Field> _facepot;
        std::shared_ptr<Field> _vertpot;
        const ScalarKernel* _kernel;
};

}

#endif
