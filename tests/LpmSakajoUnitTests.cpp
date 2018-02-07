#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmOutputMessage.h"
#include "LpmLogger.h"
#include "LpmSakajo2009.h"
#include "LpmTimer.h"
#include "LpmVtkFileIO.h"
#include "LpmMultiIndex.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <memory>
#include <fstream>
#include <mpi.h>

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

    const int nn = 10000;
    const scalar_type lat1 = PI/4.0;
    const scalar_type lat2 = PI/6.0;
    const scalar_type dlambda = 2.0*PI / nn;
    const scalar_type perturb_amp = 0.02;
    
    std::shared_ptr<SphericalCoords> sc(new SphericalCoords(2*nn));
    for (int i=0; i < nn; ++i) {
        scalar_type x1, x2, y1, y2, z1, z2;
        llToXyz(x1, y1, z1, i*dlambda, lat1);
        llToXyz(x2, y2, z2, i*dlambda, lat2);
        sc->insert(x1, y1, z1 + perturb_amp * std::sin(i * dlambda));
        sc->insert(x2, y2, z2);
    }
    const scalar_type meshSize = sphereDistance(sc->getVec(0), sc->getVec(1));
    const scalar_type nu = 1.0;
    
    std::shared_ptr<Field> circ(new Field(sc->n(), 1, "Gamma", "1/time"));
    const scalar_type Gamma = 1.0 / nn;
    for (int i=0; i < nn; ++i){
        circ->insert(Gamma);
        circ->insert(Gamma);
    }
    
    const int seriesOrder = 6;
    const int treeDepth = 4;
    const scalar_type delta = 0.1;
    SakajoTree tree(seriesOrder, treeDepth, 1.0, delta, procRank);
    
    
    std::shared_ptr<Field> velDirect(new Field(sc->n(), 3, "direct_sum_velocity", "dist/time"));
    velDirect->initializeToConstant(sc.get());
    Timer directSumTimer("Direct sum timer");
    directSumTimer.start();
    for (index_type i=0; i<sc->n(); ++i) {
        const XyzVector tgtVec = sc->getVec(i);
        XyzVector vel(0.0, 0.0, 0.0);
        for (index_type j=0; j<i; ++j) {
            const XyzVector srcVec = sc->getVec(j);
            const scalar_type Gamma = circ->getScalar(j);   
            const XyzVector kernel = tree.biotSavart(tgtVec, srcVec, delta);
            vel += kernel.scalarMultiply(Gamma);
        }
        for (index_type j=i+1; j<sc->n(); ++j) {
            const XyzVector srcVec = sc->getVec(j);
            const scalar_type Gamma = circ->getScalar(j);
            const XyzVector kernel = tree.biotSavart(tgtVec, srcVec, delta);
            vel += kernel.scalarMultiply(Gamma);
        }
        velDirect->replace(i, vel);
    }
    directSumTimer.end();
    std::cout << directSumTimer.infoString();
    
    std::vector<std::shared_ptr<Field>> fields;
    fields.push_back(circ);
    
    SakajoNode* root = tree.getRoot();
    const std::vector<MultiIndex> keys = root->getKeys();
//     for (int i=0; i<keys.size(); ++i) 
//         std::cout << keys[i] << "   " << std::endl;
        
//     std::cout << "Coefficients" << std::endl;
//     std::cout << root->coeffString();
    
//     tree.computeMoments(sc, circ);
//     std::cout << root->momentString();
    
    std::shared_ptr<Field> velTree(new Field(sc->n(), 3, "tree_sum_velocity", "dist/time"));
    velTree->initializeToConstant(sc.get());
    Timer treecodeTimer("Tree sum timer");
    treecodeTimer.start();
    tree.computeVelocity(velTree, sc, circ, meshSize, nu);
    treecodeTimer.end();
    std::cout << treecodeTimer.infoString();
    
    std::shared_ptr<Field> treecodeError(new Field(sc->n(), 3, "treecode_error", "dist/time"));
    treecodeError->initializeToConstant(sc.get());
    treecodeError->update(1.0, velDirect, -1.0, velTree);
    std::cout << "max(error) = " << treecodeError->maxMagnitude() << std::endl;
    
    
    const std::string outfname = "sakajoUnitTest.csv";
    std::ofstream fs(outfname);
    VtkWriter writer;
    writer.writePointsAndScalarFieldsToCSV(fs, sc, fields);
//     writer.writeVTKHeader(fs, "sakajo-2009");
//     writer.writeCoordsToVTKPoints(fs, sc);
//     writer.writeVTKPointDataHeader(fs, sc->n());
//     writer.writeFieldToVTKData(fs, circ);
    fs.close();
    
    
    
    tree.writeToVtk("sakajoTree.vtk");
    
    
    programTimer.end();
    OutputMessage finalMsg("PROGAM COMPLETE: " + programTimer.infoString(), OutputMessage::remarkPriority, "main");
    log->logMessage(finalMsg);
    MPI_Finalize();
return 0;
}

