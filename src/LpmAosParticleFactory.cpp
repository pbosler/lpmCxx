#include "LpmAosParticleFactory.hpp"

namespace Lpm {
namespace Aos {
    
template class ParticleFactory<2>;
template class ParticleFactory<3>;
template class SWEParticleFactory<3>;
template class SWEParticleFactory<2>;
template class BasicParticleFactory<3>;
template class BasicParticleFactory<2>;

}
}