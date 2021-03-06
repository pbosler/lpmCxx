#include "LpmDirectSum.h"
#include <mpi.h>
#include <exception>

namespace Lpm {

std::unique_ptr<Logger> DirectSum::log(new Logger(OutputMessage::debugPriority, "DirectSum_log"));

DirectSum::DirectSum(std::shared_ptr<MeshedParticles> pmesh, const std::shared_ptr<ScalarKernel> kernel, 
    const std::string tgtFieldName, const std::string srcFieldName) : _mesh(pmesh->meshPtr()),
     _vertTgt(pmesh->getVertexFieldPtr(tgtFieldName)), _faceTgt(pmesh->getFaceFieldPtr(tgtFieldName)), 
     _faceSrc(pmesh->getFaceFieldPtr(srcFieldName)), useWeights(false), _kernel(kernel) {};
    
DirectSum::DirectSum(std::shared_ptr<Particles> particles, const std::shared_ptr<ScalarKernel> kernel, 
    const std::string tgtFieldName, const std::string srcFieldName, const std::string srcWeightName) : 
    _srcLocs(particles->getCoordPtr(false)), _srcVals(particles->getFieldPtr(srcFieldName)),
    _srcWeights(particles->getFieldPtr(srcWeightName)), _tgtVals(particles->getFieldPtr(tgtFieldName)), 
    _kernel(kernel) { useWeights = !_srcWeights.expired();}

DirectSum::DirectSum(std::shared_ptr<MeshedParticles> pmesh, const std::shared_ptr<VectorKernel> kernel,
    const std::string tgtFieldName, const std::string srcFieldName) : _mesh(pmesh->meshPtr()),
    _vertTgt(pmesh->getVertexFieldPtr(tgtFieldName)), _faceTgt(pmesh->getFaceFieldPtr(tgtFieldName)),
    _faceSrc(pmesh->getFaceFieldPtr(srcFieldName)), useWeights(false), _vecKernel(kernel) {};

void DirectSum::meshSolve(const MPIReplicatedData& mpiVerts, const MPIReplicatedData& mpiFaces) const {
    const int myRank = mpiVerts.getRank();
    const int commSize = mpiVerts.getSize();
    
    if (!_kernel.expired()) { // Scalar solve
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
    else if (!_vecKernel.expired()) { // vector solve
        std::shared_ptr<Coords> vert_ptr = _mesh.lock()->getPhysCoords();
        std::shared_ptr<Faces> face_ptr = _mesh.lock()->getFaces();
        std::shared_ptr<Field> src_ptr = _faceSrc.lock();
        std::shared_ptr<Field> vert_tgt_ptr = _vertTgt.lock();
        std::shared_ptr<Field> face_tgt_ptr = _faceTgt.lock();
        std::shared_ptr<VectorKernel> kernel = _vecKernel.lock();
        
        const scalar_type mult = kernel->sumMultiplier();
        
        for (index_type i = mpiVerts.startIndex(myRank); i <= mpiVerts.endIndex(myRank); ++i) {
            XyzVector tgtVec;
            const XyzVector tgtLoc = vert_ptr->getVec(i);
            for (index_type j = 0; j < face_ptr->n(); ++j) {
                if (!face_ptr->hasChildren(j)) {
                    const XyzVector srcLoc = face_ptr->centroid(j);
                    const scalar_type strength = src_ptr->getScalar(j) * face_ptr->area(j);
                    XyzVector kval = kernel->evaluate(tgtLoc, srcLoc);
                    kval.scale(strength);
                    tgtVec += kval;
                }
            }
            tgtVec.scale(mult);
            vert_tgt_ptr->replace(i, tgtVec);
        }
        
        if (kernel->isSingular()) {
            for (index_type i = mpiFaces.startIndex(myRank); i <= mpiFaces.endIndex(myRank); ++i) {
                if (!face_ptr->hasChildren(i)) {
                    XyzVector tgtVec;
                    const XyzVector tgtLoc = face_ptr->centroid(i);
                    for (index_type j = 0; j < i; ++j) {
                        if (!face_ptr->hasChildren(j)) {
                            const XyzVector srcLoc = face_ptr->centroid(j);
                            const scalar_type strength = src_ptr->getScalar(j) * face_ptr->area(j);
                            XyzVector kval = kernel->evaluate(tgtLoc, srcLoc);
                            kval.scale(strength);
                            tgtVec += kval;
                        }
                    }
                    for (index_type j = i+1; j < face_ptr->n(); ++j) {
                        if (!face_ptr->hasChildren(j)) {
                            const XyzVector srcLoc = face_ptr->centroid(j);
                            const scalar_type strength = src_ptr->getScalar(j) * face_ptr->area(j);
                            XyzVector kval = kernel->evaluate(tgtLoc, srcLoc);
                            kval.scale(strength);
                            tgtVec += kval;
                        }
                    }
                    tgtVec.scale(mult);
                    face_tgt_ptr->replace(i, tgtVec);
                }
            }
        }
        else {
            for (index_type i = mpiFaces.startIndex(myRank); i <= mpiFaces.endIndex(myRank); ++i) {
                XyzVector tgtVec;
                const XyzVector tgtLoc = face_ptr->centroid(i);
                if (!face_ptr->hasChildren(i)) {
                    for (index_type j = 0; j < face_ptr->n(); ++j) {
                        if (!face_ptr->hasChildren(j)) {
                            const XyzVector srcLoc = face_ptr->centroid(j);
                            const scalar_type strength = src_ptr->getScalar(j) * face_ptr->area(j);
                            XyzVector kval = kernel->evaluate(tgtLoc, srcLoc);
                            kval.scale(strength);
                            tgtVec += kval;
                        }
                    }
                    tgtVec.scale(mult);
                    face_tgt_ptr->replace(i, tgtVec);
                }
            }
        }
    }
    else {
        OutputMessage errMsg("invalid kernel", OutputMessage::errorPriority, "DirectSum::meshSolve");
        log->logMessage(errMsg);
        throw std::logic_error("invalid kernel");
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
    if (!_kernel.expired()) { //scalar broadcast
        for (int i = 0; i < mpiVerts.getSize(); ++i) {
            scalar_type* bufStart0 = _vertTgt.lock()->getPtrToData();
            scalar_type* bufStart1 = _faceTgt.lock()->getPtrToData();
            MPI_Bcast(bufStart0 + mpiVerts.startIndex(i), mpiVerts.msgSize(i), MPI_DOUBLE, i, MPI_COMM_WORLD);
            MPI_Bcast(bufStart1 + mpiFaces.startIndex(i), mpiFaces.msgSize(i), MPI_DOUBLE, i, MPI_COMM_WORLD);
        
        }
    }
    else if (!_vecKernel.expired()) { // vector broadcast
        for (int i = 0; i < mpiVerts.getSize(); ++i) {
            scalar_type* vertBufStartComp0 = _vertTgt.lock()->getPtrToData(0);
            scalar_type* vertBufStartComp1 = _vertTgt.lock()->getPtrToData(1);
            scalar_type* vertBufStartComp2 = _vertTgt.lock()->getPtrToData(2);
            scalar_type* faceBufStartComp0 = _faceTgt.lock()->getPtrToData(0);
            scalar_type* faceBufStartComp1 = _faceTgt.lock()->getPtrToData(1);
            scalar_type* faceBufStartComp2 = _faceTgt.lock()->getPtrToData(2);
        
            MPI_Bcast(vertBufStartComp0 + mpiVerts.startIndex(i), mpiVerts.msgSize(i), MPI_DOUBLE, i, MPI_COMM_WORLD);
            MPI_Bcast(vertBufStartComp1 + mpiVerts.startIndex(i), mpiVerts.msgSize(i), MPI_DOUBLE, i, MPI_COMM_WORLD);
            MPI_Bcast(vertBufStartComp2 + mpiVerts.startIndex(i), mpiVerts.msgSize(i), MPI_DOUBLE, i, MPI_COMM_WORLD);
            MPI_Bcast(faceBufStartComp0 + mpiFaces.startIndex(i), mpiFaces.msgSize(i), MPI_DOUBLE, i, MPI_COMM_WORLD);
            MPI_Bcast(faceBufStartComp1 + mpiFaces.startIndex(i), mpiFaces.msgSize(i), MPI_DOUBLE, i, MPI_COMM_WORLD);
            MPI_Bcast(faceBufStartComp2 + mpiFaces.startIndex(i), mpiFaces.msgSize(i), MPI_DOUBLE, i, MPI_COMM_WORLD);
        }
    }
}


void DirectSum::meshfreeBroadcast(const MPIReplicatedData& mpiParticles) const {
    for (int i = 0; i < mpiParticles.getSize(); ++i) {
        scalar_type* bufStart = _tgtVals.lock()->getPtrToData();
        MPI_Bcast(bufStart + mpiParticles.startIndex(i), mpiParticles.msgSize(i), MPI_DOUBLE, i, MPI_COMM_WORLD);
    }
}

}