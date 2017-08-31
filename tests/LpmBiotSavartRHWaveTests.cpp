#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmOutputMessage.h"
#include "LpmLogger.h"
#include "LpmXyzVector.h"
#include "LpmMeshedParticles.h"
#include "LpmAnalyticFunctions.h"
#include "LpmMPIReplicatedData.h"
#include "LpmDirectSum.h"
#include "LpmTimer.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <memory>
#include <cmath>
#include <mpi.h>

using namespace Lpm;

XyzVector rh54velocity(const XyzVector& loc, const scalar_type t, const scalar_type bkgrdRotationRate = 0.0);

int main (int argc, char* argv[]) {
    int mpiErrCode;
    int numProcs;
    int procRank;
    mpiErrCode = MPI_Init(&argc, &argv);
    mpiErrCode = MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    mpiErrCode = MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    
    Timer programTimer("Biot-Savart convergence tests");
    programTimer.start();
    
    std::unique_ptr<Logger> log(new Logger(OutputMessage::debugPriority));
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
    const int maxRecursion = 2;
    const scalar_type srad = 1.0;
    std::unique_ptr<IcosTriSphereSeed> mSeed(new IcosTriSphereSeed());
    sphereHarmonic54 sphHarm;
    std::shared_ptr<VectorKernel> kernel(new BiotSavartSphere());
    {
        for (int k = 0; k < maxRecursion; ++k) {
            ss.str(nullstr);
            ss << "icosTri recursion level " << k;
            Timer kTimer(ss.str());
            kTimer.start();
            
            std::shared_ptr<MeshedParticles> pmesh(new MeshedParticles(*mSeed, k, true, srad, procRank));
            
            pmesh->createVertexField("vorticity", "1/time", 1);
            pmesh->createVertexField("velocity", "length/time", 3);
            pmesh->createVertexField("exact_velocity", "length/time", 3);
            pmesh->createVertexField("abs(velocity_error)", "length/time", 1);
            
            pmesh->createFaceField("vorticity", "1/time", 1);
            pmesh->createFaceField("velocity", "length/time", 3);
            pmesh->createFaceField("exact_velocity", "length/time", 3);
            pmesh->createFaceField("abs(velocity_error)", "length/time", 1);
            
            pmesh->initializeVertexFieldWithFunction("vorticity", &sphHarm);
            
            pmesh->initializeFaceFieldWithFunction("vorticity", &sphHarm);
    
            MPIReplicatedData mpiVerts(pmesh->nVertices(), procRank, numProcs);
            MPIReplicatedData mpiFaces(pmesh->nFaces(), procRank, numProcs);
            
            DirectSum solver(pmesh, kernel, "velocity", "vorticity");
            solver.meshSolve(mpiVerts, mpiFaces);
            solver.meshBroadcast(mpiVerts, mpiFaces);
                        
            kTimer.end();
        }
    }
///////////////////////////////////////////////////////////////////////////////
    programTimer.end();
    
    OutputMessage finalMsg("PROGAM COMPLETE: " + programTimer.infoString(), OutputMessage::remarkPriority, "main");
    log->logMessage(finalMsg);
    
    MPI_Finalize();
return 0;
}

XyzVector rh54velocity(const XyzVector& loc, const scalar_type t, const scalar_type bkgrdRotationRate) {
    const scalar_type lat = latitude(loc);
    const scalar_type lon = longitude(loc);
    
    const scalar_type u = 0.5 * std::cos(4.0 * lon) * std::pow(std::cos(lat), 3) * 
        (5.0 * std::cos(2.0 * lat) - 3.0);
    const scalar_type v = 4.0 * std::pow(std::cos(lat), 3) * std::sin(lat) * std::sin(4.0 * lon);
    
    const scalar_type uu = -bkgrdRotationRate * loc.y - u * std::sin(lon) - v * std::sin(lat) * std::cos(lon);
    const scalar_type vv =  bkgrdRotationRate * loc.x + u * std::cos(lon) - v * std::sin(lat) * std::sin(lon);
    const scalar_type ww = v * std::cos(lat);
    
    return XyzVector(uu, vv, ww);
}

