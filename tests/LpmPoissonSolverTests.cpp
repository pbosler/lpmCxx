#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmOutputMessage.h"
#include "LpmLogger.h"
#include "LpmAnalyticFunctions.h"
#include "LpmParticles.h"
#include "LpmMeshedParticles.h"
#include "LpmPolyMesh2d.h"
#include "LpmMPIReplicatedData.h"
#include "LpmDirectSum.h"
#include "LpmTimer.h"
#include "LpmScalarKernel.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <memory>
#include <cmath>
#include <mpi.h>

using namespace Lpm;

scalar_type appxConvergenceRate(const scalar_type dx1, const scalar_type err1, const scalar_type dx2, const scalar_type err2);

int main (int argc, char* argv[]) {
    int mpiErrCode;
    int numProcs;
    int procRank;
    mpiErrCode = MPI_Init(&argc, &argv);
    mpiErrCode = MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    mpiErrCode = MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    
    Timer programTimer("Poisson solver convergence tests");
    programTimer.start();
    
    std::unique_ptr<Logger> log(new Logger(OutputMessage::debugPriority, "Poisson solver test log", procRank));
    std::stringstream ss;
    ss << "Test info: \n \t title: " << "Poisson Solver Tests" << std::endl;
    ss << "\t objectives: " << std::endl;
    ss << "\t 1. Verify convergence of planar direct sum solver, free boundary conditions" << std::endl;
    ss << "\t 2. Verify convergence of spherical direct sum solver" << std::endl;
    OutputMessage statusMsg(ss.str(), OutputMessage::tracePriority, "main");
    log->logMessage(statusMsg);
    
    const int maxRecursion = 3;
    { // planar tests
    const exact2dpotential exactPotential;
    const radial2dsource source;
    std::shared_ptr<ScalarKernel> kernel(new PlanarGreensFnFreeBoundaries);

    const scalar_type dradius = 3.0;
    std::vector<scalar_type> triHexTimes(maxRecursion, 0.0);
    std::vector<scalar_type> quadRectTimes(maxRecursion, 0.0);
    std::vector<scalar_type> triHexLinfVerts(maxRecursion, 0.0);
    std::vector<scalar_type> triHexLinfFaces(maxRecursion, 0.0);
    std::vector<scalar_type> quadRectLinfVerts(maxRecursion, 0.0);
    std::vector<scalar_type> quadRectLinfFaces(maxRecursion, 0.0);
    std::vector<scalar_type> triHexAvgMeshSize(maxRecursion, 0.0);
    std::vector<scalar_type> quadRectAvgMeshSize(maxRecursion, 0.0);
    std::vector<scalar_type> triHexConvRateFaces(maxRecursion-1, 0.0);
    std::vector<scalar_type> triHexConvRateVerts(maxRecursion-1, 0.0);
    std::vector<scalar_type> quadRectConvRateFaces(maxRecursion-1, 0.0);
    std::vector<scalar_type> quadRectConvRateVerts(maxRecursion-1,0.0);
    {    
        statusMsg.resetMsgString("Test 1: planar triangles");
        log->logMessage(statusMsg);
        
        std::unique_ptr<TriHexSeed> mSeed(new TriHexSeed());
        for (int k = 0; k < maxRecursion; ++k) {
            ss.str(std::string());
            ss << "TriHex recursion level " << k;
            Timer timer(ss.str());
            timer.start();
            
            std::shared_ptr<MeshedParticles> pmesh(new MeshedParticles(*mSeed, k, false, dradius, procRank));
            
            pmesh->createVertexField("source", "n/a", 1);
            pmesh->createVertexField("potential", "n/a", 1);
            pmesh->createVertexField("exactPotential", "n/a", 1);
            pmesh->createVertexField("vertexError", "n/a", 1);
            
            pmesh->initializeVertexFieldWithFunction("source", &source);
            pmesh->initializeVertexFieldWithFunction("exactPotential", &exactPotential);
            
            pmesh->createFaceField("source", "n/a", 1);
            pmesh->createFaceField("potential", "n/a", 1);
            pmesh->createFaceField("exactPotential", "n/a", 1);
            pmesh->createFaceField("faceError", "n/a", 1);
            
            pmesh->initializeFaceFieldWithFunction("source", &source);
            pmesh->initializeFaceFieldWithFunction("exactPotential", &exactPotential);
            
            MPIReplicatedData mpiVerts(pmesh->nVertices(), procRank, numProcs);
            MPIReplicatedData mpiFaces(pmesh->nFaces(), procRank, numProcs);
            
            DirectSum solver(pmesh, kernel, "potential", "source");
            solver.meshSolve(mpiVerts, mpiFaces);
            solver.meshBroadcast(mpiVerts, mpiFaces);
            
            pmesh->getVertexFieldPtr("vertexError")->update(1.0, pmesh->getVertexFieldPtr("potential"), 
                                                          -1.0, pmesh->getVertexFieldPtr("exactPotential"));
            pmesh->getVertexFieldPtr("vertexError")->abs();
            
            pmesh->getFaceFieldPtr("faceError")->update(1.0, pmesh->getFaceFieldPtr("potential"),
                                                      -1.0, pmesh->getFaceFieldPtr("exactPotential"));
            pmesh->getFaceFieldPtr("faceError")->abs();
            
            if (procRank == 0 ) {
                ss.str(std::string());
                ss << "poissonTest_triHex" << k << ".vtk";
                pmesh->writeToVtkFile(ss.str(), "2d Poisson sovler test, free boundaries");
            }
            timer.end();
            OutputMessage timerMsg(timer.infoString(), OutputMessage::tracePriority, "main");
            log->logMessage(timerMsg);
            triHexTimes[k] = timer.elapsed();
            triHexLinfVerts[k] = pmesh->getVertexFieldPtr("vertexError")->maxScalarVal();
            triHexLinfFaces[k] = pmesh->getFaceFieldPtr("faceError")->maxScalarVal();
            triHexAvgMeshSize[k] = pmesh->meshPtr()->avgMeshSize();  
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
            
            std::shared_ptr<MeshedParticles> pmesh(new MeshedParticles(mSeed, k, false, dradius, procRank));
            
            pmesh->createVertexField("source", "n/a", 1);
            pmesh->createVertexField("potential", "n/a", 1);
            pmesh->createVertexField("exactPotential", "n/a", 1);
            pmesh->createVertexField("vertexError", "n/a", 1);
            
            pmesh->initializeVertexFieldWithFunction("source", &source);
            pmesh->initializeVertexFieldWithFunction("exactPotential", &exactPotential);
            
            pmesh->createFaceField("source", "n/a", 1);
            pmesh->createFaceField("potential", "n/a", 1);
            pmesh->createFaceField("exactPotential", "n/a", 1);
            pmesh->createFaceField("faceError", "n/a", 1);
            
            pmesh->initializeFaceFieldWithFunction("source", &source);
            pmesh->initializeFaceFieldWithFunction("exactPotential", &exactPotential);
            
            MPIReplicatedData mpiVerts(pmesh->nVertices(), procRank, numProcs);
            MPIReplicatedData mpiFaces(pmesh->nFaces(), procRank, numProcs);
            
            DirectSum solver(pmesh, kernel, "potential", "source");
            solver.meshSolve(mpiVerts, mpiFaces);
            solver.meshBroadcast(mpiVerts, mpiFaces);
            
            pmesh->getVertexFieldPtr("vertexError")->update(1.0, pmesh->getVertexFieldPtr("potential"), 
                                                          -1.0, pmesh->getVertexFieldPtr("exactPotential"));
            pmesh->getVertexFieldPtr("vertexError")->abs();
            
            pmesh->getFaceFieldPtr("faceError")->update(1.0, pmesh->getFaceFieldPtr("potential"),
                                                      -1.0, pmesh->getFaceFieldPtr("exactPotential"));
            pmesh->getFaceFieldPtr("faceError")->abs();
            
            if (procRank == 0) {
                ss.str(std::string());
                ss << "poissonTest_quadRect" << k << ".vtk";
                pmesh->writeToVtkFile(ss.str(), "2d Poisson sovler test, free boundaries");
            }
            
            timer.end();
            OutputMessage timerMsg(timer.infoString(), OutputMessage::tracePriority, "main");
            log->logMessage(timerMsg);
            quadRectTimes[k] = timer.elapsed();
            quadRectLinfVerts[k] = pmesh->getVertexFieldPtr("vertexError")->maxScalarVal();
            quadRectLinfFaces[k] = pmesh->getFaceFieldPtr("faceError")->maxScalarVal();
            quadRectAvgMeshSize[k] = pmesh->meshPtr()->avgMeshSize();
        }
    }
    
    mpiErrCode = MPI_Barrier(MPI_COMM_WORLD);
    
    for (int k = 1; k < maxRecursion; ++k) {
        triHexConvRateFaces[k-1] = appxConvergenceRate(triHexAvgMeshSize[k-1], triHexLinfFaces[k-1],
            triHexAvgMeshSize[k], triHexLinfFaces[k]);
        triHexConvRateVerts[k-1] = appxConvergenceRate(triHexAvgMeshSize[k-1], triHexLinfVerts[k-1],
            triHexAvgMeshSize[k], triHexLinfVerts[k]);
        quadRectConvRateFaces[k-1] = appxConvergenceRate(quadRectAvgMeshSize[k-1], quadRectLinfFaces[k-1],
            quadRectAvgMeshSize[k], quadRectLinfFaces[k]);
        quadRectConvRateVerts[k-1] = appxConvergenceRate(quadRectAvgMeshSize[k-1], quadRectLinfVerts[k-1],
            quadRectAvgMeshSize[k], quadRectLinfVerts[k]);
    }   
    
    ss.str(std::string());
    ss << "SUMMARY (nProc = " << numProcs << ")" << std::endl;
    ss << "TriHex results: " << std::endl;
    ss << "\tmesh size: ";
    for (int k = 0; k < maxRecursion; ++k)
        ss << triHexAvgMeshSize[k] << " ";
    ss << std::endl;
    ss << "\tvertex linf error: ";
    for (int k = 0; k < maxRecursion; ++k)
        ss << triHexLinfVerts[k] << " ";
    ss << std::endl;
    ss << "\tface linf error: ";
    for (int k = 0; k < maxRecursion; ++k)
        ss << triHexLinfFaces[k] << " ";
    ss << std::endl;
    ss << "\ttime (seconds) :";
    for (int k = 0; k < maxRecursion; ++k)
        ss << triHexTimes[k] << " ";
    ss << std::endl;
    ss << "\tappx conv. rates : " << std::endl;
    ss << "\t\tFaces: -- ";
    for (int i = 0; i < maxRecursion - 1; ++i) 
        ss << triHexConvRateFaces[i] << " ";
    ss << std::endl;
    ss << "\t\tVerts: -- ";
    for (int k = 0; k < maxRecursion - 1; ++k) 
        ss << triHexConvRateVerts[k] << " ";
    ss << std::endl;
    ss << "QuadRect results:" << std::endl;
    ss << "\tmesh size: ";
    for (int k = 0; k < maxRecursion; ++k)
        ss << quadRectAvgMeshSize[k] << " ";
    ss << std::endl;
    ss << "\tvertex linf error: ";
    for (int k = 0; k < maxRecursion; ++k)
        ss << quadRectLinfVerts[k] << " ";
    ss << std::endl;
    ss << "\tface linf error: ";
    for (int k = 0; k < maxRecursion; ++k)
        ss << quadRectLinfFaces[k] << " ";
    ss << std::endl;
    ss << "\ttime (seconds) :";
    for (int k = 0; k < maxRecursion; ++k)
        ss << quadRectTimes[k] << " ";
    ss << std::endl;
    ss << "\tappx. conv. rates : " << std::endl;
    ss << "\t\tFaces: -- ";
    for (int k = 0; k < maxRecursion - 1; ++k) 
        ss << quadRectConvRateFaces[k] << " ";
    ss << std::endl;
    ss << "\t\tVerts: -- ";
    for (int k =0; k < maxRecursion - 1; ++k)
        ss << quadRectConvRateVerts[k] << " ";
    ss << std::endl;
    
    OutputMessage summaryMsg(ss.str(), OutputMessage::remarkPriority, "main");
    log->logMessage(summaryMsg);
    }// end planar tests
    
    {// spherical tests
        const scalar_type sradius = 1.0;
        sphereHarmonic54 sphHarm;
        std::shared_ptr<ScalarKernel> kernel(new SphereGreensFn(sradius));
        
        std::vector<scalar_type> icosTriTimes(maxRecursion, 0.0);
        std::vector<scalar_type> icosTriLinfFaces(maxRecursion, 0.0);
        std::vector<scalar_type> icosTriLinfVerts(maxRecursion, 0.0);
        std::vector<scalar_type> icosTriAvgMeshSize(maxRecursion, 0.0);
        std::vector<scalar_type> icosTriConvRateFaces(maxRecursion-1, 0.0);
        std::vector<scalar_type> icosTriConvRateVerts(maxRecursion-1, 0.0);
        
        std::vector<scalar_type> cubedSphereTimes(maxRecursion, 0.0);
        std::vector<scalar_type> cubedSphereLinfFaces(maxRecursion, 0.0);
        std::vector<scalar_type> cubedSphereLinfVerts(maxRecursion, 0.0);
        std::vector<scalar_type> cubedSphereAvgMeshSize(maxRecursion, 0.0);
        std::vector<scalar_type> cubedSphereConvRateFaces(maxRecursion-1, 0.0);
        std::vector<scalar_type> cubedSphereConvRateVerts(maxRecursion-1, 0.0);
        { // icos tri section
            statusMsg.resetMsgString("Test 3: spherical triangles");
            log->logMessage(statusMsg);
            
            std::unique_ptr<IcosTriSphereSeed> mSeed(new IcosTriSphereSeed());
            for (int k = 0; k < maxRecursion; ++k) {
                ss.str(std::string());
                ss << "icosTri recursion level " << k;
                Timer timer(ss.str());
                timer.start();
            
                std::shared_ptr<MeshedParticles> pmesh(new MeshedParticles(*mSeed, k, false, sradius, procRank));
            
                pmesh->createVertexField("source", "n/a", 1);
                pmesh->createVertexField("potential", "n/a", 1);
                pmesh->createVertexField("exactPotential", "n/a", 1);
                pmesh->createVertexField("vertexError", "n/a", 1);
            
                pmesh->initializeVertexFieldWithFunction("source", &sphHarm);
                pmesh->initializeVertexFieldWithFunction("exactPotential", &sphHarm);
                pmesh->getVertexFieldPtr("exactPotential")->scale(1.0 / 30.0);
            
                pmesh->createFaceField("source", "n/a", 1);
                pmesh->createFaceField("potential", "n/a", 1);
                pmesh->createFaceField("exactPotential", "n/a", 1);
                pmesh->createFaceField("faceError", "n/a", 1);
            
                pmesh->initializeFaceFieldWithFunction("source", &sphHarm);
                pmesh->initializeFaceFieldWithFunction("exactPotential", &sphHarm);
                pmesh->getFaceFieldPtr("exactPotential")->scale(1.0 / 30.0);
            
                MPIReplicatedData mpiVerts(pmesh->nVertices(), procRank, numProcs);
                MPIReplicatedData mpiFaces(pmesh->nFaces(), procRank, numProcs);
    //             if (procRank == 0) {
    //                 OutputMessage mpiMsg(mpiFaces.infoString(), OutputMessage::tracePriority, "main");
    //                 log->logMessage(mpiMsg);
    //             }
            
                DirectSum solver(pmesh, kernel, "potential", "source");
                solver.meshSolve(mpiVerts, mpiFaces);
                solver.meshBroadcast(mpiVerts, mpiFaces);
            
                pmesh->getVertexFieldPtr("vertexError")->update(1.0, pmesh->getVertexFieldPtr("potential"), 
                                                              -1.0, pmesh->getVertexFieldPtr("exactPotential"));
                pmesh->getVertexFieldPtr("vertexError")->abs();
            
                pmesh->getFaceFieldPtr("faceError")->update(1.0, pmesh->getFaceFieldPtr("potential"),
                                                          -1.0, pmesh->getFaceFieldPtr("exactPotential"));
                pmesh->getFaceFieldPtr("faceError")->abs();
            
                if (procRank == 0 ) {
                    ss.str(std::string());
                    ss << "poissonTest_icosTriSphere" << k << ".vtk";
                    pmesh->writeToVtkFile(ss.str(), "Spherical Poisson sovler test");
                }
                timer.end();
                OutputMessage timerMsg(timer.infoString(), OutputMessage::tracePriority, "main");
                log->logMessage(timerMsg);
                icosTriTimes[k] = timer.elapsed();
                icosTriLinfVerts[k] = pmesh->getVertexFieldPtr("vertexError")->maxScalarVal();
                icosTriLinfFaces[k] = pmesh->getFaceFieldPtr("faceError")->maxScalarVal();
                icosTriAvgMeshSize[k] = pmesh->meshPtr()->avgMeshSize();  
            
            }
        } // end icos tri section
        { // cubed sphere section
            statusMsg.resetMsgString("Test 4: spherical quadrilaterals");
            log->logMessage(statusMsg);
            
            std::unique_ptr<CubedSphereSeed> mSeed(new CubedSphereSeed());
            for (int k = 0; k < maxRecursion; ++k) {
                ss.str(std::string());
                ss << "cubedSphere recursion level " << k;
                Timer timer(ss.str());
                timer.start();
            
                std::shared_ptr<MeshedParticles> pmesh(new MeshedParticles(*mSeed, k, false, sradius, procRank));
            
                pmesh->createVertexField("source", "n/a", 1);
                pmesh->createVertexField("potential", "n/a", 1);
                pmesh->createVertexField("exactPotential", "n/a", 1);
                pmesh->createVertexField("vertexError", "n/a", 1);
            
                pmesh->initializeVertexFieldWithFunction("source", &sphHarm);
                pmesh->initializeVertexFieldWithFunction("exactPotential", &sphHarm);
                pmesh->getVertexFieldPtr("exactPotential")->scale(1.0 / 30.0);
            
                pmesh->createFaceField("source", "n/a", 1);
                pmesh->createFaceField("potential", "n/a", 1);
                pmesh->createFaceField("exactPotential", "n/a", 1);
                pmesh->createFaceField("faceError", "n/a", 1);
            
                pmesh->initializeFaceFieldWithFunction("source", &sphHarm);
                pmesh->initializeFaceFieldWithFunction("exactPotential", &sphHarm);
                pmesh->getFaceFieldPtr("exactPotential")->scale(1.0 / 30.0);
            
                MPIReplicatedData mpiVerts(pmesh->nVertices(), procRank, numProcs);
                MPIReplicatedData mpiFaces(pmesh->nFaces(), procRank, numProcs);
    //             if (procRank == 0) {
    //                 OutputMessage mpiMsg(mpiFaces.infoString(), OutputMessage::tracePriority, "main");
    //                 log->logMessage(mpiMsg);
    //             }
            
                DirectSum solver(pmesh, kernel, "potential", "source");
                solver.meshSolve(mpiVerts, mpiFaces);
                solver.meshBroadcast(mpiVerts, mpiFaces);
            
                pmesh->getVertexFieldPtr("vertexError")->update(1.0, pmesh->getVertexFieldPtr("potential"), 
                                                              -1.0, pmesh->getVertexFieldPtr("exactPotential"));
                pmesh->getVertexFieldPtr("vertexError")->abs();
            
                pmesh->getFaceFieldPtr("faceError")->update(1.0, pmesh->getFaceFieldPtr("potential"),
                                                          -1.0, pmesh->getFaceFieldPtr("exactPotential"));
                pmesh->getFaceFieldPtr("faceError")->abs();
            
                if (procRank == 0 ) {
                    ss.str(std::string());
                    ss << "poissonTest_cubedSphere" << k << ".vtk";
                    pmesh->writeToVtkFile(ss.str(), "Spherical sovler test, free boundaries");
                }
                timer.end();
                OutputMessage timerMsg(timer.infoString(), OutputMessage::tracePriority, "main");
                log->logMessage(timerMsg);
                cubedSphereTimes[k] = timer.elapsed();
                cubedSphereLinfVerts[k] = pmesh->getVertexFieldPtr("vertexError")->maxScalarVal();
                cubedSphereLinfFaces[k] = pmesh->getFaceFieldPtr("faceError")->maxScalarVal();
                cubedSphereAvgMeshSize[k] = pmesh->meshPtr()->avgMeshSize();  
            
            }
        } // end cubed sphere section
        for (int k = 1; k < maxRecursion; ++k) {
            icosTriConvRateFaces[k-1] = appxConvergenceRate(icosTriAvgMeshSize[k-1], icosTriLinfFaces[k-1],
                icosTriAvgMeshSize[k], icosTriLinfFaces[k]);
            icosTriConvRateVerts[k-1] = appxConvergenceRate(icosTriAvgMeshSize[k-1], icosTriLinfVerts[k-1],
                icosTriAvgMeshSize[k], icosTriLinfVerts[k]);
            cubedSphereConvRateFaces[k-1] = appxConvergenceRate(cubedSphereAvgMeshSize[k-1], cubedSphereLinfFaces[k-1],
                cubedSphereAvgMeshSize[k], cubedSphereLinfFaces[k]);
            cubedSphereConvRateVerts[k-1] = appxConvergenceRate(cubedSphereAvgMeshSize[k-1], cubedSphereLinfVerts[k-1],
                cubedSphereAvgMeshSize[k], cubedSphereLinfVerts[k]);
        }
        ss.str(std::string());
        ss << "SUMMARY (nProc = " << numProcs << ")" << std::endl;
        ss << "IcosTriSphere results: " << std::endl;
        ss << "\tmesh size: ";
        for (int k = 0; k < maxRecursion; ++k)
            ss << icosTriAvgMeshSize[k] << " ";
        ss << std::endl;
        ss << "\tvertex linf error: ";
        for (int k = 0; k < maxRecursion; ++k)
            ss << icosTriLinfVerts[k] << " ";
        ss << std::endl;
        ss << "\tface linf error: ";
        for (int k = 0; k < maxRecursion; ++k)
            ss << icosTriLinfFaces[k] << " ";
        ss << std::endl;
        ss << "\ttime (seconds) :";
        for (int k = 0; k < maxRecursion; ++k)
            ss << icosTriTimes[k] << " ";
        ss << std::endl;
        ss << "\tappx conv. rates : " << std::endl;
        ss << "\t\tFaces: -- ";
        for (int i = 0; i < maxRecursion - 1; ++i) 
            ss << icosTriConvRateFaces[i] << " ";
        ss << std::endl;
        ss << "\t\tVerts: -- ";
        for (int k = 0; k < maxRecursion - 1; ++k) 
            ss << icosTriConvRateVerts[k] << " ";
        ss << std::endl;
        ss << "CubedSphere results:" << std::endl;
        ss << "\tmesh size: ";
        for (int k = 0; k < maxRecursion; ++k)
            ss << cubedSphereAvgMeshSize[k] << " ";
        ss << std::endl;
        ss << "\tvertex linf error: ";
        for (int k = 0; k < maxRecursion; ++k)
            ss << cubedSphereLinfVerts[k] << " ";
        ss << std::endl;
        ss << "\tface linf error: ";
        for (int k = 0; k < maxRecursion; ++k)
            ss << cubedSphereLinfFaces[k] << " ";
        ss << std::endl;
        ss << "\ttime (seconds) :";
        for (int k = 0; k < maxRecursion; ++k)
            ss << cubedSphereTimes[k] << " ";
        ss << std::endl;
        ss << "\tappx. conv. rates : " << std::endl;
        ss << "\t\tFaces: -- ";
        for (int k = 0; k < maxRecursion - 1; ++k) 
            ss << cubedSphereConvRateFaces[k] << " ";
        ss << std::endl;
        ss << "\t\tVerts: -- ";
        for (int k =0; k < maxRecursion - 1; ++k)
            ss << cubedSphereConvRateVerts[k] << " ";
        ss << std::endl;
    
        OutputMessage summaryMsg(ss.str(), OutputMessage::remarkPriority, "main");
        log->logMessage(summaryMsg);
    }
    programTimer.end();
    
    OutputMessage finalMsg("PROGAM COMPLETE: " + programTimer.infoString(), OutputMessage::remarkPriority, "main");
    log->logMessage(finalMsg);
    
    MPI_Finalize();
return 0;
}

scalar_type appxConvergenceRate(const scalar_type dx1, const scalar_type err1, const scalar_type dx2, const scalar_type err2) {
    const scalar_type numer = std::log(err2) - std::log(err1);
    const scalar_type denom = std::log(dx2) - std::log(dx1);
    return numer / denom;
}