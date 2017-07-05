#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmOutputMessage.h"
#include "LpmLogger.h"
#include "LpmTimer.h"
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

    std::unique_ptr<Logger> log(new Logger(OutputMessage::debugPriority, "LpmTimerUnitTests_main_log", procRank));
    {
        std::stringstream ss;
        ss << "Test info: \n \t title: " << "LpmTimer Unit Tests" << std::endl;
        ss << "\t objectives: " << std::endl;
        ss << "\t 1. verify functionality of timer class" << std::endl;
        OutputMessage introMsg(ss.str(), OutputMessage::tracePriority, "main");
        log->logMessage(introMsg);
    }

    Timer timer("unitTestTimer");
    timer.start();
    const index_type maxIter = 1e6;
    index_type j = 0;
    for (index_type i = 0; i < maxIter; ++i) {
        j += (j+1) * i;
    }
    timer.end();
    
    if (procRank == 0) {
        std::cout << "elapsed time = " << timer.elapsed() << std::endl;
    }
    
    OutputMessage statusMsg(timer.infoString(), OutputMessage::remarkPriority, "main");
    log->logMessage(statusMsg);
    
    MPI_Finalize();
return 0;
}

