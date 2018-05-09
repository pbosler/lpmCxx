#include <iostream>
#include <sstream>
#include <iomanip>
#include <memory>
#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmOutputMessage.h"
#include "LpmLogger.h"
#include "LpmAosTypes.hpp"

using namespace Lpm;

typedef Vec<2> vec2type;
typedef Vec<3> vec3type;

int main (int argc, char* argv[]) {

    std::unique_ptr<Logger> log(new Logger(OutputMessage::debugPriority));
    std::stringstream ss;
    const std::string nullstr;
    {
        ss << "Test info: \n \t title: " << "LPM Array of Structures types unit test" << std::endl;
        ss << "\t objectives: " << std::endl;
        ss << "\t 1. verify template classes' basic functions" << std::endl;
        OutputMessage introMsg(ss.str(), OutputMessage::tracePriority, "main");
        log->logMessage(introMsg);
        ss.str(nullstr);
    }

    scalar_type xx[3] = {0.0, 1.0, 2.0};

    Vec<2> v2(xx);
    
    Vec<3> v3(xx);
    
    std::cout << "vec2 = " << v2 << std::endl;
    std::cout << "vec3 = " << v3 << std::endl;

return 0;
}

