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
    
    std::shared_ptr<Coords> vert_ptr = _verts.lock();
    std::shared_ptr<Faces> face_ptr = _faces.lock();
    std::shared_ptr<Field> src_ptr = _src.lock();
    std::shared_ptr<Field> face_pot_ptr = _facepot.lock();
    std::shared_ptr<Field> vert_pot_ptr = _vertpot.lock();
    
    //
    //  loop over vertices
    //
    for (index_type i = mpiVerts.startIndex(myRank); i <= mpiVerts.endIndex(myRank); ++i) {
        scalar_type potVal = 0.0;
        const XyzVector tgtLoc = vert_ptr->getVec(i);
        for (index_type j = 0; j < face_ptr->n(); ++j) {
            if (!face_ptr->hasChildren(j)) {
                const XyzVector srcLoc = face_ptr->centroid(j);
                const scalar_type kernelWeight = src_ptr->getScalar(j) * face_ptr->area(j);
                const scalar_type kernelVal = _kernel->evaluate(tgtLoc, srcLoc);
                potVal += kernelVal * kernelWeight;
            }
        }
        vert_pot_ptr->replace(i, potVal);
    }
    
    //
    //  loop over faces
    //
    for (index_type i = mpiFaces.startIndex(myRank); i <= mpiFaces.endIndex(myRank); ++i) {
        if (!face_ptr->hasChildren(i)) {
            scalar_type potVal = 0.0;
            const XyzVector tgtLoc = face_ptr->centroid(i);
            for (index_type j = 0; j < i; ++j) {
                if (!face_ptr->hasChildren(j)) {
                    const XyzVector srcLoc = face_ptr->centroid(j);
                    const scalar_type kernelWeight = src_ptr->getScalar(j) * face_ptr->area(j);
                    const scalar_type kernelVal = _kernel->evaluate(tgtLoc, srcLoc);
                    potVal += kernelVal * kernelWeight;
                }
            }
            for (index_type j = i+1; j < face_ptr->n(); ++j) {
                if (!face_ptr->hasChildren(j)) {
                    const XyzVector srcLoc = face_ptr->centroid(j);
                    const scalar_type kernelWeight = src_ptr->getScalar(j) * face_ptr->area(j);
                    const scalar_type kernelVal = _kernel->evaluate(tgtLoc, srcLoc);
                    potVal += kernelVal * kernelWeight;
            
                }        
            }
            face_pot_ptr->replace(i, potVal);
        }
    }
}

void PoissonSolverDirectSum::broadcastSolution(const MPIReplicatedData& mpiVerts, const MPIReplicatedData& mpiFaces) const {
    //
    //  broadcast vertex solution
    //
    for (int i = 0; i < mpiVerts.getSize(); ++i) {
        scalar_type* bufStart = _vertpot.lock()->getPtrToData();
        MPI_Bcast(bufStart + mpiVerts.startIndex(i), mpiVerts.msgSize(i), MPI_DOUBLE, i, MPI_COMM_WORLD);
    }
    //
    //  broadcast face solution
    //
    for (int i = 0; i < mpiFaces.getSize(); ++i) {
        scalar_type* bufStart = _facepot.lock()->getPtrToData();
        MPI_Bcast(bufStart + mpiFaces.startIndex(i), mpiFaces.msgSize(i), MPI_DOUBLE, i, MPI_COMM_WORLD);
    }
}

}
