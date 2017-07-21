#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmOutputMessage.h"
#include "LpmLogger.h"
#include "LpmMultiIndex.h"
#include "LpmXyzVector.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <memory>

using namespace Lpm;

int main (int argc, char* argv[]) {

    std::unique_ptr<Logger> log(new Logger(OutputMessage::debugPriority));
    std::stringstream ss;
    const std::string nullstr;
    {
        
        ss << "Test info: \n \t title: " << "MultiIndex Unit Tests" << std::endl;
        ss << "\t objectives: " << std::endl;
        ss << "\t 1. Verify class MultiIndex implements the standard multi-index operations" << std::endl;
        OutputMessage introMsg(ss.str(), OutputMessage::tracePriority, "main");
        log->logMessage(introMsg);
    }

    const MultiIndex k(4, 3, 1);
    const XyzVector vec(2.0, 3.0, 0.5);
    
    ss.str(nullstr);
    ss << "MultiIndex k = " << k << ":" << std::endl;
    ss << "\tmagnitude (8) = " << k.magnitude() << std::endl;
    ss << "\tvectorPower (216): " << vec << "^k = " << k.vectorPower(vec) << std::endl;
    ss << "\tfactorial (144) = " << k.factorial() << std::endl;
    
    OutputMessage statusMsg(ss.str(), OutputMessage::remarkPriority, "main");
    log->logMessage(statusMsg);
    
    
return 0;
}

