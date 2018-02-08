#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmOutputMessage.h"
#include "LpmLogger.h"
#include "LpmSakajo2009.h"
#include "LpmTimer.h"
#include "LpmVtkFileIO.h"
#include "LpmMultiIndex.h"
#include "LpmOctree.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <memory>
#include <fstream>
#include <mpi.h>
#include <cmath>

using namespace Lpm;

int main (int argc, char* argv[]) {
    int mpiErrCode;
    int numProcs;
    int procRank;
    mpiErrCode = MPI_Init(&argc, &argv);
    mpiErrCode = MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    mpiErrCode = MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    
    Timer programTimer("Sakajo 2009 Unit Tests");
    programTimer.start();

    std::unique_ptr<Logger> log(new Logger(OutputMessage::debugPriority, "SakajoUnitTests_log", procRank));
    std::stringstream ss;
    const std::string nullstr;
    {
        ss << "Test info: \n \t title: " << "Unit tests for Sakajo 2009 data structures" << std::endl;
//         ss << "\t objectives: " << std::endl;
//         ss << "\t 1. *objective 1*" << std::endl;
        OutputMessage introMsg(ss.str(), OutputMessage::tracePriority, "main");
        log->logMessage(introMsg);
        ss.str(nullstr);
    }

    const int nLat = 12;
    const int nLon = 12;
    const int seriesOrder = 2;
    const int treeDepthParam = 1;
    const scalar_type meshSize = 1.0 / (std::pow(2, treeDepthParam));
    const scalar_type nu = 1.0 / 6;
    const scalar_type delta = 0.1;
    
    const scalar_type dz = 1.0 / (0.5 * nLat);
    
    std::shared_ptr<SphericalCoords> sc(new SphericalCoords(4096));
    for (int i=0; i<nLat; ++i) {
        const scalar_type zz = 1.0 - (i+1) * dz;
        for (int j=0; j < nLon; ++j) {
            const scalar_type x = std::sqrt(1.0-square(zz)) * std::cos(2.0 * PI * j / nLon);
            const scalar_type y = std::sqrt(1.0-square(zz)) * std::sin(2.0 * PI * j / nLon);
            sc->insert(x,y,zz);
//             std::cout << "getVec(" << i+j << ") : (" << x << ", " << y << ", " << zz << ")  =  " << sc->getVec(i+j) << std::endl;
        }
    }
    std::cout << "sc->n(): " << sc->n() << std::endl;
    std::shared_ptr<Field> circ(new Field(sc->n(), 1, "Gamma", "1/time"));
    for (index_type i=0; i<sc->n(); ++i) {
        circ->insert(std::pow(-1.0, i));
    }
    
    
    SakajoTree tree(seriesOrder, treeDepthParam, 1.0, delta, procRank);
    
    std::cout << tree.Tree::infoString();
    
    std::shared_ptr<Field> velDirect(new Field(sc->n(), 3, "direct_sum_velocity", "dist/time"));
    velDirect->initializeToConstant(sc.get());
    Timer directSumTimer("Direct sum timer");
    directSumTimer.start();
    for (index_type i=0; i<sc->n(); ++i) {
        const XyzVector tgtVec = sc->getVec(i);
        XyzVector vel(0.0, 0.0, 0.0);
        for (index_type j=0; j<sc->n(); ++j) {
            const XyzVector srcVec = sc->getVec(j);
            const scalar_type Gamma = circ->getScalar(j);   
            const XyzVector kernel = tree.biotSavart(tgtVec, srcVec, delta);
            vel += kernel.scalarMultiply(Gamma);
        }
//         for (index_type j=i+1; j<sc->n(); ++j) {
//             const XyzVector srcVec = sc->getVec(j);
//             const scalar_type Gamma = circ->getScalar(j);
//             const XyzVector kernel = tree.biotSavart(tgtVec, srcVec, delta);
//             vel += kernel.scalarMultiply(Gamma);
//         }
        velDirect->replace(i, vel);
    }
    directSumTimer.end();
    std::cout << directSumTimer.infoString();
    
    std::shared_ptr<Field> velTree(new Field(sc->n(), 3, "tree_sum_velocity", "dist/time"));
    velTree->initializeToConstant(sc.get());
    
    
    tree.writeToVtk("sakajoTree.vtk");
    Timer treecodeTimer("Tree sum timer");
    treecodeTimer.start();
    tree.computeMoments(sc, circ);
    
    tree.computeVelocity(velTree, sc, circ, meshSize, nu);
    tree.printAll();
    
    treecodeTimer.end();
    std::cout << treecodeTimer.infoString();
    
    std::shared_ptr<Field> treecodeError(new Field(sc->n(), 3, "treecode_error", "dist/time"));
    treecodeError->initializeToConstant(sc.get());
    treecodeError->update(1.0, velDirect, -1.0, velTree);
    std::cout << "max(error) = " << treecodeError->maxMagnitude() << std::endl;
   
    std::vector<std::shared_ptr<Field>> fields = {circ, velDirect, velTree, treecodeError};
    const std::string outfname = "sakajoUnitTest.csv";
    std::ofstream fs(outfname);
    VtkWriter writer;
    writer.writePointsAndFieldsToCSV(fs, sc, fields);  
    fs.close();  
    
    programTimer.end();
    OutputMessage finalMsg("PROGAM COMPLETE: " + programTimer.infoString(), OutputMessage::remarkPriority, "main");
    log->logMessage(finalMsg);
    MPI_Finalize();
return 0;
}

