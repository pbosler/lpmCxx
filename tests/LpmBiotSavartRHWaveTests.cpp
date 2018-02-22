#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmOutputMessage.h"
#include "LpmLogger.h"
#include "LpmXyzVector.h"
#include "LpmMeshedParticles.h"
#include "LpmAnalyticFunctions.h"
#include "LpmMPIReplicatedData.h"
#include "LpmDirectSum.h"
#include "LpmTreeSum.h"
#include "LpmTimer.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <memory>
#include <cmath>
#include <mpi.h>

using namespace Lpm;

//XyzVector rh54velocity(const XyzVector& loc, const scalar_type t=0.0, const scalar_type bkgrdRotationRate = 0.0);

int main (int argc, char* argv[]) {
    int mpiErrCode;
    int numProcs;
    int procRank;
    mpiErrCode = MPI_Init(&argc, &argv);
    mpiErrCode = MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    mpiErrCode = MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    
    Timer programTimer("Biot-Savart convergence tests");
    programTimer.start();
    
    std::unique_ptr<Logger> log(new Logger(OutputMessage::debugPriority, "main", procRank));
    std::stringstream ss;
    const std::string nullstr;
    {
        ss << "Test info: \n \t title: " << "TEST_TITLE" << std::endl;
        ss << "\t objectives: " << std::endl;
        ss << "\t 1. *objective 1*" << std::endl;
        OutputMessage introMsg(ss.str(), OutputMessage::tracePriority, "main");
        log->logMessage(introMsg);
        ss.str(nullstr);
    }
///////////////////////////////////////////////////////////////////////////////
    const int maxRecursion = 6;
    const scalar_type srad = 1.0;
    const int maxSeriesOrder = 6;
    const int nParticlesPerBox = 500;
    const scalar_type mpTol = 0.33;
    std::unique_ptr<IcosTriSphereSeed> mSeed(new IcosTriSphereSeed());
    sphereHarmonic54 sphHarm;
    rossbyHaurwitz54Velocity rhwave;
    std::shared_ptr<VectorKernel> kernel(new BiotSavartSphere());
    {
        for (int k = 0; k < maxRecursion; ++k) {
            ss.str(nullstr);
            ss << "icosTri recursion level " << k;
            
            ss.str(nullstr);
            ss << "biotSavart_icosTri" << k << ".vtk";
            const std::string ofile = ss.str();
            ss.str(nullstr);
            
            Timer kTimer(ss.str());
            kTimer.start();
            
            std::shared_ptr<MeshedParticles> pmesh(new MeshedParticles(*mSeed, k, true, srad, procRank));
            
            pmesh->createVertexField("vorticity", "1/time", 1);
            pmesh->createVertexField("velocity", "length/time", 3);
            pmesh->createVertexField("treeVelocity", "length/time", 3);
            pmesh->createVertexField("exact_velocity", "length/time", 3);
            pmesh->createVertexField("velocity_error_direct", "length/time", 3);
            pmesh->createVertexField("tree_error", "length/time", 3);
            
            pmesh->createFaceField("vorticity", "1/time", 1);
            pmesh->createFaceField("velocity", "length/time", 3);
            pmesh->createFaceField("treeVelocity", "length/time", 3);
            pmesh->createFaceField("exact_velocity", "length/time", 3);
            pmesh->createFaceField("velocity_error_direct", "length/time", 3);
            pmesh->createFaceField("tree_error", "length/time", 3);
            
            pmesh->initializeVertexFieldWithFunction("vorticity", &sphHarm);
            pmesh->initializeFaceFieldWithFunction("vorticity", &sphHarm);
            pmesh->initializeVertexFieldWithFunction("exact_velocity", &rhwave);
            pmesh->initializeFaceFieldWithFunction("exact_velocity", &rhwave);
     
            MPIReplicatedData mpiVerts(pmesh->nVertices(), procRank, numProcs);
            MPIReplicatedData mpiFaces(pmesh->nFaces(), procRank, numProcs);
            
            Timer directSolveTimer("direct_sum_timer");
            directSolveTimer.start();
            DirectSum solver(pmesh, kernel, "velocity", "vorticity");
            solver.meshSolve(mpiVerts, mpiFaces);
            solver.meshBroadcast(mpiVerts, mpiFaces);
            directSolveTimer.end();
            
            std::shared_ptr<Field> errPtr = pmesh->getVertexFieldPtr("velocity_error_direct");
            errPtr->linearOp(1.0, pmesh->getVertexFieldPtr("velocity"), -1.0, pmesh->getVertexFieldPtr("exact_velocity"));
            errPtr = pmesh->getFaceFieldPtr("velocity_error_direct");
            errPtr->linearOp(1.0, pmesh->getFaceFieldPtr("velocity"), -1.0, pmesh->getFaceFieldPtr("exact_velocity"));
            
            if (procRank == 0) {
                std::cout << directSolveTimer.infoString();
                std::cout << "max rel. err. (direct) = " << pmesh->getFaceFieldPtr("velocity_error_direct")->maxMagnitude()
                    / pmesh->getFaceFieldPtr("exact_velocity")->maxMagnitude() << std::endl;
            }
            
            
                        
            Timer treecodeTimer("treecode_timer");
            treecodeTimer.start();
            TreeSum treeSolver(maxSeriesOrder, nParticlesPerBox, pmesh, kernel, "treeVelocity", "vorticity");
            treeSolver.computeMeshMoments();
            int myMacCounter = 0;
            treeSolver.meshSolve(mpiVerts, mpiFaces, mpTol, myMacCounter);
            treeSolver.meshBroadcast(mpiVerts, mpiFaces);
            treecodeTimer.end();
            int macCounter = 0;
            mpiErrCode = MPI_Reduce(&myMacCounter, &macCounter, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
            
            errPtr = pmesh->getVertexFieldPtr("tree_error");
            errPtr->linearOp(1.0, pmesh->getVertexFieldPtr("treeVelocity"), -1.0, pmesh->getVertexFieldPtr("velocity"));
            errPtr = pmesh->getFaceFieldPtr("tree_error");
            errPtr->linearOp(1.0, pmesh->getFaceFieldPtr("treeVelocity"), -1.0, pmesh->getFaceFieldPtr("velocity"));
            
            if (procRank==0) {
                std::cout << "mac: " << macCounter << "/" << pmesh->nVertices() + pmesh->nLeafFaces() << std::endl;
                std::cout <<  treecodeTimer.infoString();
                std::cout << "max rel err. treecode = " << pmesh->getFaceFieldPtr("tree_error")->maxMagnitude() /
                    pmesh->getFaceFieldPtr("exact_velocity")->maxMagnitude() << std::endl;
            }
            
            ss << "tree_n0" << nParticlesPerBox<< "_icosTri" << k << ".vtk";
            std::string treefile(ss.str());
            ss.str(nullstr);
            if (procRank == 0) treeSolver.writeToVtk(treefile);
            
            
            
            if (procRank ==0)
                pmesh->writeToVtkFile(ofile);
            
            kTimer.end();
            if (procRank==0)
                std::cout << kTimer.infoString();
        }
    }
///////////////////////////////////////////////////////////////////////////////
    programTimer.end();
    if (procRank==0) {
        OutputMessage finalMsg("PROGAM COMPLETE: " + programTimer.infoString(), OutputMessage::remarkPriority, "main");
        log->logMessage(finalMsg);
    }
    
    MPI_Finalize();
return 0;
}

// XyzVector rh54velocity(const XyzVector& loc, const scalar_type t, const scalar_type bkgrdRotationRate) {
//     const scalar_type lat = latitude(loc);
//     const scalar_type lon = longitude(loc);
//     
//     const scalar_type u = 0.5 * std::cos(4.0 * lon) * std::pow(std::cos(lat), 3) * 
//         (5.0 * std::cos(2.0 * lat) - 3.0);
//     const scalar_type v = 4.0 * std::pow(std::cos(lat), 3) * std::sin(lat) * std::sin(4.0 * lon);
//     
//     const scalar_type uu = -bkgrdRotationRate * loc.y - u * std::sin(lon) - v * std::sin(lat) * std::cos(lon);
//     const scalar_type vv =  bkgrdRotationRate * loc.x + u * std::cos(lon) - v * std::sin(lat) * std::sin(lon);
//     const scalar_type ww = v * std::cos(lat);
//     
//     return XyzVector(uu, vv, ww);
// }

