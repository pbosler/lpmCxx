#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmOutputMessage.h"
#include "LpmLogger.h"
#include "LpmSphericalCoords.h"
#include "LpmTimer.h"
#include "LpmScalarKernel.h"
#include "LpmTaylorSeries3d.h"
#include "LpmAnalyticFunctions.h"
#include "LpmMeshedParticles.h"
#include "LpmMeshSeed.h"
#include "LpmDirectSum.h"
#include "LpmMPIReplicatedData.h"
#include "LpmTreeSum.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <memory>
#include <mpi.h>

using namespace Lpm;

int main (int argc, char* argv[]) {
    int mpiErrCode;
    int numProcs;
    int procRank;
    mpiErrCode = MPI_Init(&argc, &argv);
    mpiErrCode = MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    mpiErrCode = MPI_Comm_rank(MPI_COMM_WORLD, &procRank);

    std::unique_ptr<Logger> log(new Logger(OutputMessage::debugPriority, "TreeSumTest_log", procRank));
    std::stringstream ss;
    const std::string nullstr;
    {
        ss << "Test info: \n \t title: " << "Tree Sum tests" << std::endl;
        ss << "\t objectives: " << std::endl;
        ss << "\t 1. Verify basic functionality of TreeSum class and SumNode struct." << std::endl;
        ss << "\t 2. Verify accuracy and order of accuracy of treecode accelerated summation." << std::endl;
        ss << "\t 3. Verify efficiency gain of treecode accelerated summation vs. direct summation." << std::endl;
        OutputMessage introMsg(ss.str(), OutputMessage::tracePriority, "main");
        log->logMessage(introMsg);
        ss.str(nullstr);
    }
    
    const int maxRecursion = 3;
    const scalar_type sphRadius = 1.0;
    const sphereHarmonic54 sphHarm;
    std::shared_ptr<ScalarKernel> kernel(new SphereGreensFn(sphRadius));
    
    IcosTriSphereSeed seed;
    
    std::shared_ptr<MeshedParticles> sphere(new MeshedParticles(seed, maxRecursion, false, sphRadius, procRank));
    
    sphere->createVertexField("source", "n/a", 1);
    sphere->createVertexField("exactPotential", "n/a", 1);
    sphere->createVertexField("direct_potential", "n/a", 1);
    sphere->createVertexField("tree_potential", "n/a", 1);
    sphere->createVertexField("direct_error", "n/a", 1);
    sphere->createVertexField("tree_error", "n/a", 1);
    sphere->createVertexField("tree-direct", "n/a", 1);
    
    sphere->createFaceField("source", "n/a", 1);
    sphere->createFaceField("exactPotential", "n/a", 1);
    sphere->createFaceField("direct_potential", "n/a", 1);
    sphere->createFaceField("tree_potential", "n/a", 1);
    sphere->createFaceField("direct_error", "n/a", 1);
    sphere->createFaceField("tree_error", "n/a", 1);
    sphere->createFaceField("tree-direct", "n/a", 1);
    
    sphere->initializeVertexFieldWithFunction("source", &sphHarm);
    sphere->initializeVertexFieldWithFunction("exactPotential", &sphHarm);
    sphere->getVertexFieldPtr("exactPotential")->scale(1.0 / 30.0);
    
    sphere->initializeFaceFieldWithFunction("source", &sphHarm);
    sphere->initializeFaceFieldWithFunction("exactPotential", &sphHarm);
    sphere->getFaceFieldPtr("exactPotential")->scale(1.0 / 30.0);
    
    MPIReplicatedData mpiVerts(sphere->nVertices(), procRank, numProcs);
    MPIReplicatedData mpiFaces(sphere->nFaces(), procRank, numProcs);
    
    Timer directTimer("Direct sum timer");
    directTimer.start();
    DirectSum directSolver(sphere, kernel, "direct_potential", "source");
    directSolver.meshSolve(mpiVerts, mpiFaces);
    directSolver.meshBroadcast(mpiVerts, mpiFaces);
    directTimer.end();
    
    sphere->getVertexFieldPtr("direct_error")->update(1.0, sphere->getVertexFieldPtr("direct_potential"), 
        -1.0, sphere->getVertexFieldPtr("exactPotential"));
    sphere->getVertexFieldPtr("direct_error")->abs();
    //sphere->getVertexFieldPtr("direct_error")->scale(1.0 / sphere->getVertexFieldPtr("exactPotential")->maxScalarVal());
    sphere->getFaceFieldPtr("direct_error")->update(1.0, sphere->getFaceFieldPtr("direct_potential"), 
        -1.0, sphere->getFaceFieldPtr("exactPotential"));
    sphere->getFaceFieldPtr("direct_error")->abs();
    //sphere->getFaceFieldPtr("direct_error")->scale(1.0 / sphere->getFaceFieldPtr("exactPotential")->maxScalarVal());
    
    ss << "Direct sum complete." << std::endl;
    ss << directTimer.infoString();
    ss << "\tl_inf error (vertices) : " << sphere->getVertexFieldPtr("direct_error")->maxScalarVal() << std::endl;
    ss << "\tl_inf error (faces) : " << sphere->getFaceFieldPtr("direct_error")->maxScalarVal() << std::endl;
    OutputMessage directSolveMsg(ss.str(), OutputMessage::remarkPriority, "main");
    log->logMessage(directSolveMsg);
    ss.str(nullstr);
    
    Timer treeSetupTimer("Tree setup timer");
    const scalar_type maxAspect = 1.5;
    const int maxP = 4;
    const scalar_type farFieldParam = 0.5;
    const scalar_type seriesParam = 0.0;
    const index_type maxCoordsPerNode = 50;
    treeSetupTimer.start();
    TreeSum tree(sphere, kernel, "tree_potential", "source", maxAspect, maxP, seriesParam, procRank, farFieldParam);
//     OutputMessage debugMsg("returned from tree constructor: "  + tree.infoString(), OutputMessage::tracePriority, "main");
//     log->logMessage(debugMsg);
    tree.buildTreeFromMesh(maxCoordsPerNode);
    treeSetupTimer.end();
    
    Timer treeSolveTimer("Tree solve timer");
    treeSolveTimer.start();
    tree.meshSolve(mpiVerts, mpiFaces);
    tree.meshBroadcast(mpiVerts, mpiFaces);
    treeSolveTimer.end();
    
    sphere->getVertexFieldPtr("tree_error")->update(1.0, sphere->getVertexFieldPtr("tree_potential"), 
        -1.0, sphere->getVertexFieldPtr("exactPotential"));
    sphere->getVertexFieldPtr("tree_error")->abs();
    sphere->getFaceFieldPtr("tree_error")->update(1.0, sphere->getFaceFieldPtr("tree_potential"),
        -1.0, sphere->getFaceFieldPtr("exactPotential"));
    sphere->getFaceFieldPtr("tree_error")->abs();
    sphere->getVertexFieldPtr("tree-direct")->update(1.0, sphere->getVertexFieldPtr("tree_potential"), 
        -1.0, sphere->getVertexFieldPtr("direct_potential"));
    sphere->getFaceFieldPtr("tree-direct")->update(1.0, sphere->getFaceFieldPtr("tree_potential"), 
        -1.0, sphere->getFaceFieldPtr("direct_potential"));
    
    ss << "Tree sum complete." << std::endl;
    ss << treeSetupTimer.infoString();
    ss << treeSolveTimer.infoString();
    ss << "\ttotal time: " << treeSetupTimer.elapsed() + treeSolveTimer.elapsed() << std::endl;
    ss << "\tl_inf error (vertices) : " << sphere->getVertexFieldPtr("tree_error")->maxScalarVal() << std::endl;
    ss << "\tl_inf error (faces) : " << sphere->getFaceFieldPtr("tree_error")->maxScalarVal() << std::endl;
    OutputMessage treeSolveMsg(ss.str(), OutputMessage::remarkPriority, "main");
    log->logMessage(treeSolveMsg);
    ss.str(nullstr);
    
    if (procRank == 0) {
        ss << "treecodeTest_icosTriSphere" << maxRecursion << ".vtk";
        sphere->writeToVtkFile(ss.str(), "Poisson solver, direct sum vs. treecode");
    }
    MPI_Finalize();
return 0;
}

