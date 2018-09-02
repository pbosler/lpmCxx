#ifndef LPM_AOS_MESHSEED_FACTORY_HPP
#define LPM_AOS_MESHSEED_FACTORY_HPP

#include "LpmAosMeshSeed.hpp"
#include <memory>
#include <string>
#include <exception>
#include <iostream>
#include <sstream>

namespace Lpm {
namespace Aos {

class MeshSeedFactory {
    public :
        std::unique_ptr<MeshSeed> createSeed(const std::string id) const {
            if (id == "triHexPlane") {
                return std::unique_ptr<MeshSeed>(new TriHexSeed());
            }
            else if (id == "quadRectPlane" ) {
                return std::unique_ptr<MeshSeed>(new QuadRectSeed());
            }
            else if (id == "icosTriSphere") {
                return std::unique_ptr<MeshSeed>(new IcosTriSphereSeed());
            }
            else if (id == "cubedSphere") {
                return std::unique_ptr<MeshSeed>(new CubedSphereSeed());
            }
            else {
                std::ostringstream ss;
                ss << "error: unrecognized seed id " << id;
                throw std::invalid_argument(ss.str());
            }
        }
};
  
}
}
#endif
