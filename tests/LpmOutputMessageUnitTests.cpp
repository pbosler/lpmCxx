#include "LpmTypeDefs.h"
#include "LpmConfig.h"
#include "LpmOutputMessage.h"
#include <iostream>
#include <vector>
#include <string>
#include <sstream>

using namespace Lpm;

int main(int argc, char* argv[]) {

    std::stringstream ss;
    std::string s;
    ss << "-- LPM Version " << LPM_VERSION_MAJOR << "." << LPM_VERSION_MINOR << " -- \n ";
    ss << "   Unit Test " << argv[0] << ": covers OutputMessage.h and Cmake Config.";
    OutputMessage introMsg( ss.str(), OutputMessage::remarkPriority, "main");
    
    std::cout << introMsg;
    
    std::vector<std::string> content = {"string 1", "string 2"};
    LongMessage longMsg("LongMessage", OutputMessage::remarkPriority, "main", content);
    
    std::cout << longMsg;   

return 0;
}