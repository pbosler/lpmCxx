#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmOutputMessage.h"
#include "LpmLogger.h"
#include "LpmParticles.h"
#include "LpmPolyMesh2d.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <memory>

using namespace Lpm;

class radial2dsource : public AnalyticFunction {
    public:
        radial2dsource() {};
    
        scalar_type evaluateScalar(const scalar_type x, const scalar_type y, const scalar_type z = 0.0) const {
            const scalar_type rr = std::sqrt(x*x + y*y);
            return (std::abs(rr) <= 1.0 ? -0.5 * PI *(std::sin(PI * rr) / rr + std::cos(PI * rr)) : 0.0);
        };
        scalar_type evaluateScalar(const XyzVector& crdVec) const {
            return evaluateScalar(crdVec.x, crdVec.y);
        };
        XyzVector evaluateVector(const scalar_type x, const scalar_type y, const scalar_type z = 0.0) const {
            return XyzVector(0.0, 0.0, 0.0);
        }
        XyzVector evaluateVector(const XyzVector& crdVec) const {
            return XyzVector(0.0, 0.0, 0.0);
        }
}; 

class exact2dpotential : public AnalyticFunction {
    public:
        exact2dpotential() {};
        
        scalar_type evaluateScalar(const scalar_type x, const scalar_type y, const scalar_type z = 0.0) const {
            const scalar_type rr = std::sqrt(x*x + y*y);
            return (std::abs(rr) <= 1.0 ? 0.5 * (1.0 + std::cos(PI * rr)) : 0.0);
        };
        scalar_type evaluateScalar(const XyzVector& crdVec) const {
            return evaluateScalar(crdVec.x, crdVec.y);
        };
        XyzVector evaluateVector(const scalar_type x, const scalar_type y, const scalar_type z = 0.0) const {
            return XyzVector(0.0, 0.0, 0.0);
        }
        XyzVector evaluateVector(const XyzVector& crdVec) const {
            return XyzVector(0.0, 0.0, 0.0);
        }
};

int main (int argc, char* argv[]) {

    std::unique_ptr<Logger> log(new Logger(OutputMessage::debugPriority));
    std::stringstream ss;
    ss << "Test info: \n \t title: " << "Poisson Solver Tests" << std::endl;
    ss << "\t objectives: " << std::endl;
    ss << "\t 1. Verify convergence of planar direct sum solver, free boundary conditions" << std::endl;
    ss << "\t 2. Verify convergence of spherical direct sum solver" << std::endl;
    OutputMessage statusMsg(ss.str(), OutputMessage::tracePriority, "main");
    log->logMessage(statusMsg);
    
    const exact2dpotential exactPotential;
    const radial2dsource source;
    
    const int maxRecursion = 6;
    const scalar_type dradius = 3.0;
    {    
        statusMsg.resetMsgString("Test 1: planar triangles");
        log->logMessage(statusMsg);
        
        std::unique_ptr<TriHexSeed> mSeed(new TriHexSeed());
        for (int k = 0; k < maxRecursion; ++k) {
            std::shared_ptr<PolyMesh2d> mesh(new PolyMesh2d(*mSeed, k, false, dradius));
            
            Particles particles(mesh);
            particles.createField("source", "n/a", 1);
            particles.createField("potential", "n/a", 1);
            particles.createField("exactPotential", "n/a", 1);
            particles.createField("error", "n/a", 1);
            
            particles.initializeFieldWithFunction("source", &source);
            particles.initializeFieldWithFunction("exactPotential", &exactPotential);
        }

    }
    {
        statusMsg.resetMsgString("Test 2: planar quadrilaterals");
        log->logMessage(statusMsg);
        
        QuadRectSeed mSeed;
        for (int k = 0; k < maxRecursion; ++k) {
            std::cout << "Mesh recursion " << k << std::endl;
            std::shared_ptr<PolyMesh2d> mesh(new PolyMesh2d(mSeed, k, false, dradius));
            
            Particles particles(mesh);
            particles.createField("source", "n/a", 1);
            particles.createField("potential", "n/a", 1);
            particles.createField("exactPotential", "n/a", 1);
            particles.createField("error", "n/a", 1);
            
            particles.initializeFieldWithFunction("source", &source);
            particles.initializeFieldWithFunction("exactPotential", &exactPotential);
        }
    }

return 0;
}

