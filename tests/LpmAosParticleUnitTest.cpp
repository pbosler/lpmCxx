#include <iostream>
#include <sstream>
#include <iomanip>
#include <memory>
#include <vector>
#include <exception>
#include <cmath>
#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmOutputMessage.h"
#include "LpmLogger.h"
#include "LpmAosTypes.hpp"
#include "LpmAosParticle.hpp"
#include "LpmAosParticleFactory.hpp"
#include "LpmAosSWEParticle.hpp"
#include "LpmAosParticleSet.hpp"
#ifdef HAVE_VTK

#endif

using namespace Lpm::Aos;

using Lpm::index_type;
using Lpm::scalar_type;
using Lpm::Logger;
using Lpm::OutputMessage;
using Lpm::ZERO_TOL;
using Lpm::PI;

int main (int argc, char* argv[]) {

    std::unique_ptr<Logger> log(new Logger(OutputMessage::debugPriority));
    std::stringstream ss;
    const std::string nullstr;
    {
        ss << "Test info: \n \t title: " << "LPM Array of Structures particle types unit test" << std::endl;
        ss << "\t objectives: " << std::endl;
        ss << "\t 1. verify particle aos classes' basic functions" << std::endl;
        ss << "\t 2. verify particleSet aos classes' basic functions" << std::endl;
        OutputMessage introMsg(ss.str(), OutputMessage::tracePriority, "main");
        log->logMessage(introMsg);
        ss.str(nullstr);
    }

    BasicParticleFactory<2> planeFactory;

    std::shared_ptr<ParticleFactory<3>> sphereFactory(new SWEParticleFactory<3>());

    std::unique_ptr<Particle<2>> base_particle = planeFactory.createParticle();
    std::unique_ptr<Particle<3>> swe_particle = sphereFactory->createParticle();

    Vec<2> oneone(1.0, 1.0);

    Vec<3> sphpoint(1.0, 1.0, 1.0);
    sphpoint.normalizeInPlace();

    std::cout << "base_particle->infoString()" << std::endl;
    base_particle->init(oneone);
    std::cout << base_particle->infoString() << std::endl;

    std::cout << "swe_particle->infoString()" << std::endl;
    swe_particle->init(sphpoint);
    std::cout << swe_particle->infoString() << std::endl;

    ParticleSet<3> cubedSphere(sphereFactory, 100);
    cubedSphere.initFromParticleSetFile("cubed_sphere/cubedSphere_1.vtk");
    std::cout << cubedSphere.infoString() << std::endl;

    if (std::abs(cubedSphere.totalWeight() - 4.0*PI) > ZERO_TOL) {
        throw std::runtime_error("Surface area unit test FAILED");
    }

#ifdef HAVE_VTK
	vtkSmartPointer<vtkPoints> vpts = cubedSphere.toVtkPoints();
	scalar_type bounds[6];
	vpts->ComputeBounds();
	vpts->GetBounds(bounds);
	std::cout << "vtk computed bounds: (xmin, xmax), (ymin, ymax), (zmin, zmax) = " << std::endl;
	std::cout << "\t(" << bounds[0] << ", " << bounds[1] << ") (" << bounds[2] << ", " << bounds[3] 
		<< ") (" << bounds[4] << ", " << bounds[5] << ")" << std::endl;
#endif

//     std::vector<std::string> allinfo = cubedSphere.particlesInfoStrings();
//     for (int i=0; i<allinfo.size(); ++i) {
//         std::cout << "particle " << i << std::endl;
//         std::cout << allinfo[i];
//     }

return 0;
}
