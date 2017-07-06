#include "LpmPoissonSolverDirectSum.h"
#include "LpmMPIReplicatedData.h"
#include "LpmTimer.h"
#include "LpmScalarKernel.h"
#include "LpmXyzVector.h"
#include <mpi.h>

namespace Lpm {

void PoissonSolverDirectSum::solve(const MPIReplicatedData& mpiVerts, const MPIReplicatedData& mpiFaces) const {
    const int myRank = mpiFaces.getRank();
    const int commSize = mpiFaces.getSize();
    
    //
    //  loop over vertices
    //
    for (index_type i = mpiVerts.startIndex(myRank); i <= mpiVerts.endIndex(myRank); ++i) {
        scalar_type potVal = 0.0;
        const XyzVector tgtLoc = _verts->getVec(i);
        for (index_type j = 0; j < _faces->n(); ++j) {
            if (!_faces->hasChildren(j)) {
                const XyzVector srcLoc = _faces->centroid(j);
                const scalar_type kernelWeight = _src->getScalar(j) * _faces->area(j);
                const scalar_type kernelVal = _kernel->evaluate(tgtLoc, srcLoc);
                potVal += kernelVal * kernelWeight;
            }
        }
        _vertpot->replace(i, potVal);
    }
    
    //
    //  loop over faces
    //
    for (index_type i = mpiFaces.startIndex(myRank); i <= mpiFaces.endIndex(myRank); ++i) {
        if (!_faces->hasChildren(i)) {
            scalar_type potVal = 0.0;
            const XyzVector tgtLoc = _faces->centroid(i);
            for (index_type j = 0; j < i; ++j) {
                if (!_faces->hasChildren(j)) {
                    const XyzVector srcLoc = _faces->centroid(j);
                    const scalar_type kernelWeight = _src->getScalar(j) * _faces->area(j);
                    const scalar_type kernelVal = _kernel->evaluate(tgtLoc, srcLoc);
                    potVal += kernelVal * kernelWeight;
                }
            }
            for (index_type j = i+1; j < _faces->n(); ++j) {
                if (!_faces->hasChildren(j)) {
                    const XyzVector srcLoc = _faces->centroid(j);
                    const scalar_type kernelWeight = _src->getScalar(j) * _faces->area(j);
                    const scalar_type kernelVal = _kernel->evaluate(tgtLoc, srcLoc);
                    potVal += kernelVal * kernelWeight;
            
                }        
            }
            _facepot->replace(i, potVal);
        }
    }
}

void PoissonSolverDirectSum::broadcastSolution(const MPIReplicatedData& mpiVerts, const MPIReplicatedData& mpiFaces) const {
    //
    //  broadcast vertex solution
    //
    for (int i = 0; i < mpiVerts.getSize(); ++i) {
        scalar_type* bufStart = _vertpot->getPtrToData();
        MPI_Bcast(bufStart + mpiVerts.startIndex(i), mpiVerts.msgSize(i), MPI_DOUBLE, i, MPI_COMM_WORLD);
    }
    //
    //  broadcast face solution
    //
    for (int i = 0; i < mpiFaces.getSize(); ++i) {
        scalar_type* bufStart = _facepot->getPtrToData();
        MPI_Bcast(bufStart + mpiFaces.startIndex(i), mpiFaces.msgSize(i), MPI_DOUBLE, i, MPI_COMM_WORLD);
    }
}

}
