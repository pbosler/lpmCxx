#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmOutputMessage.h"
#include "LpmLogger.h"
#include "LpmSphericalCoords.h"
#include "LpmScalarKernel.h"
#include "LpmTreeSum.h"
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
        ss << "Test info: \n \t title: " << "Tree sum unit tests" << std::endl;
        ss << "\t objectives: " << std::endl;
        ss << "\t 1. Verify Treesum data structure" << std::endl;
        ss << "\t 2. Verify Taylor coefficient computation" << std::endl;
        ss << "\t 3. Verify moment computations " << std::endl;
        ss << "\t 4. Verify tree summation of a scalar" << std::endl;
        OutputMessage introMsg(ss.str(), OutputMessage::tracePriority, "main");
        log->logMessage(introMsg);
        ss.str(nullstr);
    }
    
    const int nMax = 5000;
    std::shared_ptr<SphericalCoords> sc(new SphericalCoords(nMax));

    sc->initRandom();
    
    const scalar_type maxAspectRatio = 1.5;
    const int maxNCoordsPerCluster = 15;
    const int maxSeriesOrder = 4;
    const scalar_type sphRadius = 1.0;
    const SphereGreensFn kernel(sphRadius);
    
    std::shared_ptr<TreeSumNode> tree(new TreeSumNode(sc, maxAspectRatio, kernel, maxAspectRatio));
    
    generateTree(tree, sc, maxNCoordsPerCluster);

return 0;
}

