#include "LpmDirectSum.h"
#include <mpi.h>

namespace Lpm {

DirectSum::DirectSum(std::shared_ptr<MeshedParticles> pmesh, const std::shared_ptr<ScalarKernel> kernel, 
    const std::string tgtFieldName, const std::string srcFieldName) : _mesh(pmesh->meshPtr()),
     _vertTgt(pmesh->getVertexFieldPtr(tgtFieldName)), _faceTgt(pmesh->getFaceFieldPtr(tgtFieldName)), 
     _faceSrc(pmesh->getFaceFieldPtr(srcFieldName)), useWeights(false), _kernel(kernel) {};
    
DirectSum::DirectSum(std::shared_ptr<Particles> particles, const std::shared_ptr<ScalarKernel> kernel, 
    const std::string tgtFieldName, const std::string srcFieldName, const std::string srcWeightName) : 
    _srcLocs(particles->getCoordPtr(false)), _srcVals(particles->getFieldPtr(srcFieldName)),
    _srcWeights(particles->getFieldPtr(srcWeightName)), _tgtVals(particles->getFieldPtr(tgtFieldName)), 
    _kernel(kernel) { useWeights = !_srcWeights.expired();}

void DirectSum::meshSolve(const MPIReplicatedData& mpiVerts, const MPIReplicatedData& mpiFaces) const {
    const int myRank = mpiVerts.getRank();
    const int commSize = mpiVerts.getSize();
    
    std::shared_ptr<Coords> vert_ptr = _mesh.lock()->getPhysCoords();
    std::shared_ptr<Faces> face_ptr = _mesh.lock()->getFaces();
    std::shared_ptr<Field> src_ptr = _faceSrc.lock();
    std::shared_ptr<Field> vert_tgt_ptr = _vertTgt.lock();
    std::shared_ptr<Field> face_tgt_ptr = _faceTgt.lock();
    std::shared_ptr<ScalarKernel> kernel = _kernel.lock();

    for (index_type i = mpiVerts.startIndex(myRank); i <= mpiVerts.endIndex(myRank); ++i) {
        scalar_type tgtVal = 0.0;
        const XyzVector tgtLoc = vert_ptr->getVec(i);
        for (index_type j = 0; j < face_ptr->n(); ++j) {
            if (!face_ptr->hasChildren(j)) {
                const XyzVector srcLoc = face_ptr->centroid(j);
                const scalar_type strength = src_ptr->getScalar(j) * face_ptr->area(j);
                const scalar_type kval = kernel->evaluate(tgtLoc, srcLoc);
                tgtVal += kval * strength;
            }
        }
        vert_tgt_ptr->replace(i, tgtVal);
    }
    
    if (kernel->isSingular()) {
        for (index_type i = mpiFaces.startIndex(myRank); i <= mpiFaces.endIndex(myRank); ++i) {
            if (!face_ptr->hasChildren(i)) {
                scalar_type tgtVal = 0.0;
                const XyzVector tgtLoc = face_ptr->centroid(i);
                for (index_type j = 0; j < i; ++j) {
                    if (!face_ptr->hasChildren(j)) {
                        const XyzVector srcLoc = face_ptr->centroid(j);
                        const scalar_type strength = src_ptr->getScalar(j) * face_ptr->area(j);
                        const scalar_type kval = kernel->evaluate(tgtLoc, srcLoc);
                        tgtVal += kval * strength;
                    }
                }
                for (index_type j = i+1; j < face_ptr->n(); ++j) {
                    if (!face_ptr->hasChildren(j)) {
                        const XyzVector srcLoc = face_ptr->centroid(j);
                        const scalar_type strength = src_ptr->getScalar(j) * face_ptr->area(j);
                        const scalar_type kval = kernel->evaluate(tgtLoc, srcLoc);
                        tgtVal += kval * strength;
                    }
                }
                face_tgt_ptr->replace(i, tgtVal);
            }
        }
    }
    else {
        for (index_type i = mpiFaces.startIndex(myRank); i <= mpiFaces.endIndex(myRank); ++i) {
            if (!face_ptr->hasChildren(i)) {
                scalar_type tgtVal = 0.0;
                const XyzVector tgtLoc = face_ptr->centroid(i);
                for (index_type j = 0; j < face_ptr->n(); ++j) {
                    if (!face_ptr->hasChildren(j)) {
                        const XyzVector srcLoc = face_ptr->centroid(j);
                        const scalar_type strength = src_ptr->getScalar(j) * face_ptr->area(j);
                        const scalar_type kval = kernel->evaluate(tgtLoc, srcLoc);
                        tgtVal += kval * strength;
                    }
                }
                face_tgt_ptr->replace(i, tgtVal);
            }
        }
    }

}

void DirectSum::meshfreeSolve(const MPIReplicatedData& mpiParticles) const {
    const int myRank = mpiParticles.getRank();
    const int commSize = mpiParticles.getSize();
    
    std::shared_ptr<Coords> crd_ptr = _srcLocs.lock();
    std::shared_ptr<Field> src_ptr = _srcVals.lock();
    std::shared_ptr<Field> tgt_ptr = _tgtVals.lock();
    std::shared_ptr<Field> weight_ptr;
    std::shared_ptr<ScalarKernel> kernel = _kernel.lock();
    
    if (useWeights) {
        weight_ptr = _srcWeights.lock();
        if (kernel->isSingular()) {
            for (index_type i = mpiParticles.startIndex(myRank); i <= mpiParticles.endIndex(i); ++i) {
                const XyzVector tgtLoc = crd_ptr->getVec(i);
                scalar_type tgtVal = 0.0;
                for (index_type j = 0; j < i; ++j) {
                    const XyzVector srcLoc = crd_ptr->getVec(j);
                    const scalar_type strength = src_ptr->getScalar(j) * weight_ptr->getScalar(j);
                    const scalar_type kval = kernel->evaluate(tgtLoc, srcLoc);
                    tgtVal += kval * strength;
                }
                for (index_type j = i+1; j < crd_ptr->n(); ++j) {
                    const XyzVector srcLoc = crd_ptr->getVec(j);
                    const scalar_type strength = src_ptr->getScalar(j) * weight_ptr->getScalar(j);
                    const scalar_type kval = kernel->evaluate(tgtLoc, srcLoc);
                    tgtVal += kval * strength;
                }
                tgt_ptr->replace(i, tgtVal);
            }
        }
        else {
            for (index_type i = mpiParticles.startIndex(myRank); i <= mpiParticles.endIndex(myRank); ++i) {
                const XyzVector tgtLoc = crd_ptr->getVec(i);
                scalar_type tgtVal = 0.0;
                for (index_type j = 0; j < crd_ptr->n(); ++j) {
                    const XyzVector srcLoc = crd_ptr->getVec(j);
                    const scalar_type strength = src_ptr->getScalar(j) * weight_ptr->getScalar(j);
                    const scalar_type kval = kernel->evaluate(tgtLoc, srcLoc);
                    tgtVal += kval * strength;
                }
                tgt_ptr->replace(i, tgtVal);
            }
        }
    }
    else {
        if (kernel->isSingular()) {
            for (index_type i = mpiParticles.startIndex(myRank); i <= mpiParticles.endIndex(i); ++i) {
                const XyzVector tgtLoc = crd_ptr->getVec(i);
                scalar_type tgtVal = 0.0;
                for (index_type j = 0; j < i; ++j) {
                    const XyzVector srcLoc = crd_ptr->getVec(j);
                    const scalar_type strength = src_ptr->getScalar(j);
                    const scalar_type kval = kernel->evaluate(tgtLoc, srcLoc);
                    tgtVal += kval * strength;
                }
                for (index_type j = i+1; j < crd_ptr->n(); ++j) {
                    const XyzVector srcLoc = crd_ptr->getVec(j);
                    const scalar_type strength = src_ptr->getScalar(j);
                    const scalar_type kval = kernel->evaluate(tgtLoc, srcLoc);
                    tgtVal += kval * strength;
                }
                tgt_ptr->replace(i, tgtVal);
            }
        }
        else {
            for (index_type i = mpiParticles.startIndex(myRank); i <= mpiParticles.endIndex(myRank); ++i) {
                const XyzVector tgtLoc = crd_ptr->getVec(i);
                scalar_type tgtVal = 0.0;
                for (index_type j = 0; j < crd_ptr->n(); ++j) {
                    const XyzVector srcLoc = crd_ptr->getVec(j);
                    const scalar_type strength = src_ptr->getScalar(j);
                    const scalar_type kval = kernel->evaluate(tgtLoc, srcLoc);
                    tgtVal += kval * strength;
                }
                tgt_ptr->replace(i, tgtVal);
            }
        }
    }
}

void DirectSum::meshBroadcast(const MPIReplicatedData& mpiVerts, const MPIReplicatedData& mpiFaces) const {
    for (int i = 0; i < mpiVerts.getSize(); ++i) {
        scalar_type* bufStart0 = _vertTgt.lock()->getPtrToData();
        scalar_type* bufStart1 = _faceTgt.lock()->getPtrToData();
        MPI_Bcast(bufStart0 + mpiVerts.startIndex(i), mpiVerts.msgSize(i), MPI_DOUBLE, i, MPI_COMM_WORLD);
        MPI_Bcast(bufStart1 + mpiFaces.startIndex(i), mpiFaces.msgSize(i), MPI_DOUBLE, i, MPI_COMM_WORLD);
        
    }
}

void DirectSum::meshfreeBroadcast(const MPIReplicatedData& mpiParticles) const {
    for (int i = 0; i < mpiParticles.getSize(); ++i) {
        scalar_type* bufStart = _tgtVals.lock()->getPtrToData();
        MPI_Bcast(bufStart + mpiParticles.startIndex(i), mpiParticles.msgSize(i), MPI_DOUBLE, i, MPI_COMM_WORLD);
    }
}

}