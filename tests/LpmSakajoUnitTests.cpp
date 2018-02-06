#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmOutputMessage.h"
#include "LpmLogger.h"
#include "LpmSakajo2009.h"
#include "LpmTimer.h"
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

    const int nn = 1000;
    const scalar_type lat1 = PI/4.0;
    const scalar_type lat2 = PI/6.0;
    const scalar_type dlambda = 2.0*PI / nn;
    
    std::shared_ptr<SphericalCoords> sc(new SphericalCoords(2*nn));
    for (int i=0; i < nn; ++i) {
        scalar_type x1, x2, y1, y2, z1, z2;
        llToXyz(x1, y1, z1, i*dlambda, lat1);
        llToXyz(x2, y2, z2, i*dlambda, lat2);
        sc->insert(x1, y1, z1);
        sc->insert(x2, y2, z2);
    }
    std::ofstream fs("sakajoUnitTest.csv");
    sc->writeCoordsCSV(fs);
    fs.close();
    
    const int seriesOrder = 3;
    const int treeDepth = 4;
    const scalar_type delta = 0.0;
    
    SakajoTree tree(seriesOrder, treeDepth, 1.0, delta, procRank);
    tree.writeToVtk("sakajoTree.vtk");
    
    
    programTimer.end();
    OutputMessage finalMsg("PROGAM COMPLETE: " + programTimer.infoString(), OutputMessage::remarkPriority, "main");
    log->logMessage(finalMsg);
    MPI_Finalize();
return 0;
}

