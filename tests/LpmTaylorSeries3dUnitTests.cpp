#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmOutputMessage.h"
#include "LpmLogger.h"
#include "LpmTaylorSeries3d.h"
#include "LpmXyzVector.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <memory>

using namespace Lpm;

int main (int argc, char* argv[]) {

    std::unique_ptr<Logger> log(new Logger(OutputMessage::debugPriority));
    {
        std::stringstream ss;
        ss << "Test info: \n \t title: " << "Taylor series 3d unit tests" << std::endl;
        ss << "\t objectives: " << std::endl;
        ss << "\t 1. *objective 1*" << std::endl;
        OutputMessage introMsg(ss.str(), OutputMessage::tracePriority, "main");
        log->logMessage(introMsg);
    }

    const int maxOrder = 4;
    
    
    XyzVector tgtVec(1.0, 1.0, 1.0);
    XyzVector srcVec(-1.0, 1.0, 1.0);
    tgtVec.normalize();
    srcVec.normalize();

    std::cout << "tgt vector = " << tgtVec << std::endl;
    std::cout << "src vector = " << srcVec << std::endl;
    
    SphereGreensCoeffs ts(maxOrder);
    ts.computeCoeffs(tgtVec, srcVec);
    
    std::cout << ts.infoString();
return 0;
}

