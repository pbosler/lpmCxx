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

template <int ndim> class MeshSeedFactory {
    public :
        std::unique_ptr<MeshSeed> createSeed(const std::string id) const {
            switch (id) {
                case "triHexPlane" {
                    return std::unique_ptr<MeshSeed>(new TriHexSeed());
                    break;
                }
                case "quadRectPlane" {
                    return std::unique_ptr<MeshSeed>(new QuadRectSeed());
                    break;
                }
                case "icosTriSphere" {
                    return std::unique_ptr<MeshSeed>(new IcosTriSphereSeed());
                    break;
                }
                case "cubedSphere" {
                    return std::unique_ptr<MeshSeed>(new CubedSphereSeed());
                    break;
                }
                default : {
                    std::ostringstream ss;
                    ss << "error: unrecognized seed id " << id;
                    throw std::invalid_argument(ss.str());
                    break;
                }
            }
        }
};
  
}
}
#endif
