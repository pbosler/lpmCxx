#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmOutputMessage.h"
#include "LpmLogger.h"
#include "LpmParticles.h"
#include "LpmMeshedParticles.h"
#include "LpmPolyMesh2d.h"
#include "LpmMPIReplicatedData.h"
#include "LpmPoissonSolverDirectSum.h"
#include "LpmTimer.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <memory>
#include <cmath>
#include <mpi.h>

using namespace Lpm;

class radial2dsource : public AnalyticFunction {
    public:
        radial2dsource() {};
    
        scalar_type evaluateScalar(const scalar_type x, const scalar_type y, const scalar_type z = 0.0) const {
            const scalar_type rr = std::sqrt(x*x + y*y);
            const scalar_type rrdenom = rr / (ZERO_TOL * ZERO_TOL + rr * rr);
            return (std::abs(rr) <= 1.0 ? -0.5 * PI * (std::sin(PI * rr) * rrdenom + PI * std::cos(PI * rr)) : 0.0);
        };
        scalar_type evaluateScalar(const XyzVector& crdVec) const {
            return evaluateScalar(crdVec.x, crdVec.y);
        };
        XyzVector evaluateVector(const scalar_type x, const scalar_type y, const scalar_type z = 0.0) const {
            return XyzVector(0.0, 0.0, 0.0);
        }
        XyzVector evaluateVector(const XyzVector& crdVec) const {
            return XyzVector(0.0, 0.0, 0.0);
        }
}; 

class exact2dpotential : public AnalyticFunction {
    public:
        exact2dpotential() {};
        
        scalar_type evaluateScalar(const scalar_type x, const scalar_type y, const scalar_type z = 0.0) const {
            const scalar_type rr = std::sqrt(x*x + y*y);
            return (std::abs(rr) <= 1.0 ? 0.5 * (1.0 + std::cos(PI * rr)) : 0.0);
        };
        scalar_type evaluateScalar(const XyzVector& crdVec) const {
            return evaluateScalar(crdVec.x, crdVec.y);
        };
        XyzVector evaluateVector(const scalar_type x, const scalar_type y, const scalar_type z = 0.0) const {
            return XyzVector(0.0, 0.0, 0.0);
        }
        XyzVector evaluateVector(const XyzVector& crdVec) const {
            return XyzVector(0.0, 0.0, 0.0);
        }
};

int main (int argc, char* argv[]) {
    int mpiErrCode;
    int numProcs;
    int procRank;
    mpiErrCode = MPI_Init(&argc, &argv);
    mpiErrCode = MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    mpiErrCode = MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    
    std::unique_ptr<Logger> log(new Logger(OutputMessage::debugPriority, "Poisson solver test log", procRank));
    std::stringstream ss;
    ss << "Test info: \n \t title: " << "Poisson Solver Tests" << std::endl;
    ss << "\t objectives: " << std::endl;
    ss << "\t 1. Verify convergence of planar direct sum solver, free boundary conditions" << std::endl;
    ss << "\t 2. Verify convergence of spherical direct sum solver" << std::endl;
    OutputMessage statusMsg(ss.str(), OutputMessage::tracePriority, "main");
    log->logMessage(statusMsg);
    
    const exact2dpotential exactPotential;
    const radial2dsource source;
    
    const int maxRecursion = 6;
    const scalar_type dradius = 3.0;
    {    
        statusMsg.resetMsgString("Test 1: planar triangles");
        log->logMessage(statusMsg);
        
        std::unique_ptr<TriHexSeed> mSeed(new TriHexSeed());
        for (int k = 0; k < maxRecursion; ++k) {
            ss.str(std::string());
            ss << "TriHex recursion level " << k;
            Timer timer(ss.str());
            timer.start();
            
            MeshedParticles pmesh(*mSeed, k, false, dradius);
            
            pmesh.createVertexField("source", "n/a", 1);
            pmesh.createVertexField("potential", "n/a", 1);
            pmesh.createVertexField("exactPotential", "n/a", 1);
            pmesh.createVertexField("vertexError", "n/a", 1);
            
            pmesh.initializeVertexFieldWithFunction("source", &source);
            pmesh.initializeVertexFieldWithFunction("exactPotential", &exactPotential);
            
            pmesh.createFaceField("source", "n/a", 1);
            pmesh.createFaceField("potential", "n/a", 1);
            pmesh.createFaceField("exactPotential", "n/a", 1);
            pmesh.createFaceField("faceError", "n/a", 1);
            
            pmesh.initializeFaceFieldWithFunction("source", &source);
            pmesh.initializeFaceFieldWithFunction("exactPotential", &exactPotential);
            
            MPIReplicatedData mpiVerts(pmesh.nVertices(), procRank, numProcs);
            MPIReplicatedData mpiFaces(pmesh.nFaces(), procRank, numProcs);
            if (procRank == 0) {
                OutputMessage mpiMsg(mpiFaces.infoString(), OutputMessage::tracePriority, "main");
                log->logMessage(mpiMsg);
            }
            
            PoissonSolverDirectSum solver(pmesh.meshPtr(), pmesh.getFaceFieldPtr("source"), 
                pmesh.getVertexFieldPtr("potential"), pmesh.getFaceFieldPtr("potential"));
            solver.solve(mpiVerts, mpiFaces);
            
            if (procRank == 0 ) {
                ss.str(std::string());
                ss << "poissonTest_triHex" << k << ".vtk";
                pmesh.writeToVtkFile(ss.str(), "2d Poisson sovler test, free boundaries");
            }
            timer.end();
            OutputMessage timerMsg(timer.infoString(), OutputMessage::tracePriority, "main");
            log->logMessage(timerMsg);
        }

    }
    {
        statusMsg.resetMsgString("Test 2: planar quadrilaterals");
        log->logMessage(statusMsg);
        
        QuadRectSeed mSeed;
        for (int k = 0; k < maxRecursion; ++k) {
            ss.str(std::string());
            ss << "QuadRect_recursion_level_" << k;
            Timer timer(ss.str());
            timer.start();
            
            MeshedParticles pmesh(mSeed, k, false, dradius);
            
            pmesh.createVertexField("source", "n/a", 1);
            pmesh.createVertexField("potential", "n/a", 1);
            pmesh.createVertexField("exactPotential", "n/a", 1);
            pmesh.createVertexField("vertexError", "n/a", 1);
            
            pmesh.initializeVertexFieldWithFunction("source", &source);
            pmesh.initializeVertexFieldWithFunction("exactPotential", &exactPotential);
            
            pmesh.createFaceField("source", "n/a", 1);
            pmesh.createFaceField("potential", "n/a", 1);
            pmesh.createFaceField("exactPotential", "n/a", 1);
            pmesh.createFaceField("faceError", "n/a", 1);
            
            pmesh.initializeFaceFieldWithFunction("source", &source);
            pmesh.initializeFaceFieldWithFunction("exactPotential", &exactPotential);
            
            MPIReplicatedData mpiVerts(pmesh.nVertices(), procRank, numProcs);
            MPIReplicatedData mpiFaces(pmesh.nFaces(), procRank, numProcs);
            if (procRank == 0) {
                OutputMessage mpiMsg(mpiFaces.infoString(), OutputMessage::tracePriority, "main");
                log->logMessage(mpiMsg);
            }
            
            PoissonSolverDirectSum solver(pmesh.meshPtr(), pmesh.getFaceFieldPtr("source"), 
                pmesh.getVertexFieldPtr("potential"), pmesh.getFaceFieldPtr("potential"));;
            solver.solve(mpiVerts, mpiFaces);
            
            if (procRank == 0) {
                ss.str(std::string());
                ss << "poissonTest_quadRect" << k << ".vtk";
                pmesh.writeToVtkFile(ss.str(), "2d Poisson sovler test, free boundaries");
            }
            
            timer.end();
            OutputMessage timerMsg(timer.infoString(), OutputMessage::tracePriority, "main");
            log->logMessage(timerMsg);
        }
    }

    MPI_Finalize();
return 0;
}

