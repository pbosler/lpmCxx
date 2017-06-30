#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmOutputMessage.h"
#include "LpmLogger.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <memory>

using namespace Lpm;

int main (int argc, char* argv[]) {

    std::unique_ptr<Logger> log(new Logger(OutputMessage::debugPriority));
    {
        std::stringstream ss;
        ss << "Test info: \n \t title: " << "TEST_TITLE" << std::endl;
        ss << "\t objectives: " << std::endl;
        ss << "\t 1. *objective 1*" << std::endl;
        OutputMessage introMsg(ss.str(), OutputMessage::tracePriority, "main");
        log->logMessage(introMsg);
    }

return 0;
}

