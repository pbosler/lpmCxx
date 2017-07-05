#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmLogger.h"
#include "LpmOutputMessage.h"
#include <memory>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <mpi.h>

using namespace Lpm;

int main (int argc, char* argv[] ) {
    int mpiErrCode;
    int numProcs;
    int procRank;
    mpiErrCode = MPI_Init(&argc, &argv);
    mpiErrCode = MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    mpiErrCode = MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    std::unique_ptr<Logger> log = std::unique_ptr<Logger>(new Logger(OutputMessage::debugPriority,"Logger_unitTest_log", procRank));
    
    {
    std::stringstream ss;
    std::string s;
    ss << "-- LPM Version " << LPM_VERSION_MAJOR << "." << LPM_VERSION_MINOR << " -- \n ";
    ss << "   Unit Test " << argv[0] << ": covers Logger.h";
    OutputMessage introMsg( ss.str(), OutputMessage::remarkPriority, "main");

    log->logMessage(introMsg);
    }
    
    {
    OutputMessage statusMsg("... end test.", OutputMessage::remarkPriority, "main");
    log->logMessage(statusMsg);
    }
    
    {
        OutputMessage errMsg("test error ... just a test.", OutputMessage::errorPriority, "main");
        log->logMessage(errMsg);
    }
    
    MPI_Finalize();
return 0;
}
