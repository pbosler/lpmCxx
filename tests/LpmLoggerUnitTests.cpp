#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmLogger.h"
#include "LpmOutputMessage.h"
#include <memory>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>

using namespace Lpm;

int main (int argc, char* argv[] ) {
    std::unique_ptr<Logger> log = std::unique_ptr<Logger>(new Logger(OutputMessage::debugPriority,"Logger_unitTest_log"));
    
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
return 0;
}
