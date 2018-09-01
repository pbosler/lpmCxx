#ifndef LPM_AOS_QG_PARTICLE_HPP
#define LPM_AOS_QG_PARTICLE_HPP

#include "LpmAosTypes.hpp"
#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmAosParticle.hpp"

namespace Lpm {
namespace Aos {

class QGParticle : public Particle<2> {
    public:
        QGParticle() : Particle<2>() {
            this->registerVectorField("velocity");
            this->registerScalarField("potvort");
            this->registerScalarField("relvort");
            this->registerScalarField("streamfn");
        }
        
        QGParticle(const Vec<2>& pos, const scalar_type ar=0.0) : Particle<2>(pos, ar);
        
        QGParticle(const Vec<2>& xx, const Vec<2>& aa, const scalar_type ar=0.0) : Particle<2>(xx, aa, ar);
        

};

inline scalar_type coriolis(const Vec<2>& pos, const scalar_type f0 = 0.0, const scalar_type beta = 0.0) {
    return f0 + beta*pos.x[1];
}

}
}
#endif
