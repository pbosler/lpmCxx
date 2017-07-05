#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmOutputMessage.h"
#include "LpmLogger.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <memory>
#include <mpi.h>
#include "LpmMPIReplicatedData.h"

using namespace Lpm;

int main (int argc, char* argv[]) {
    int mpiErrCode;
    int numProcs;
    int procRank;
    mpiErrCode = MPI_Init(&argc, &argv);
    mpiErrCode = MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    mpiErrCode = MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    
    std::unique_ptr<Logger> log(new Logger(OutputMessage::debugPriority, "LpmMPIReplicatedDataUnitTests_main_log", procRank));
    {
        std::stringstream ss;
        ss << "Test info: \n \t title: " << "MPIReplicatedData unit tests" << std::endl;
        ss << "\t objectives: " << std::endl;
        ss << "\t 1. Verify functionality of MPIReplicatedData class" << std::endl;
        OutputMessage introMsg(ss.str(), OutputMessage::tracePriority, "main");
        log->logMessage(introMsg);
    }
    
    const int nn = 1055;
    
    MPIReplicatedData mpi(nn, procRank, numProcs);
    
    OutputMessage statusMsg(mpi.infoString(), OutputMessage::remarkPriority, "main");
    log->logMessage(statusMsg);

    MPI_Finalize();
return 0;
}

