#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmOutputMessage.h"
#include "LpmLogger.h"
#include "LpmTaylorSeries3d.h"
#include "LpmXyzVector.h"
#include "LpmField.h"
#include "LpmSphericalCoords.h"
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
        ss << "\t 1. Check basic functions of TaylorSeries3d class" << std::endl;
        OutputMessage introMsg(ss.str(), OutputMessage::tracePriority, "main");
        log->logMessage(introMsg);
    }

    const int maxOrder = 4;
    
    XyzVector tgtVec(1.0, 1.0, 1.0);
    XyzVector srcVec(-1.0, 1.0, 1.0);
    tgtVec.normalize();
    srcVec.normalize();
    const XyzVector cntd = sphereMidpoint(tgtVec, srcVec);
    const std::vector<index_type> crdInds = {0, 1};
    
    std::shared_ptr<Coords> sc(new SphericalCoords(2));
    sc->insert(tgtVec);
    sc->insert(srcVec);
    std::cout << "Coords initialized: " << std::endl;
    std::cout << sc->listAllCoords();
    
    std::shared_ptr<Field> one(new Field(2, 1, "const_1", "n_a"));
    one->initializeToConstant(sc.get(), 2 * PI);
    std::cout << "Field initialized: " << std::endl;
    std::cout << one->listFieldValues();
    
    SphereGreensSeries ts(maxOrder);

    ts.computeCoeffs(tgtVec, srcVec);

    ts.computeMoments(sc, crdInds, cntd,  one);

    std::cout << ts.infoString();
    
    std::cout << "Series sum = " << ts.sum() << std::endl;
    
    
return 0;
}

