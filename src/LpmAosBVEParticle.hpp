#ifndef LPM_AOS_BVE_PARTICLE_HPP
#define LPM_AOS_BVE_PARTICLE_HPP

#include "LpmAosTypes.hpp"
#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmAosParticle.hpp"
#include <cmath>

namespace Lpm {
namespace Aos {

template <int ndim=3> class BVEParticle : public Particle<ndim> {
    public:
        BVEParticle(const bool isvert=false) : 
        	Particle<ndim>(isvert) {
            if (ndim > 1) {
                this->registerVectorField("velocity");
            }
            else {
                this->registerScalarField("velocity");
            }
            this->registerScalarField("potvort");
            this->registerScalarField("relvort");
        }

        BVEParticle(const Vec<ndim>& pos, const scalar_type wgt = 0.0, const bool isvert=false) : 
        	Particle<ndim>(pos, wgt, isvert) {
            if (ndim > 1) {
                this->registerVectorField("velocity");
            }
            else {
                this->registerScalarField("velocity");
            }
            this->registerScalarField("potvort");
            this->registerScalarField("relvort");
        }
        
        BVEParticle(const Vec<ndim>& xx, const Vec<ndim>& aa, const scalar_type wgt = 0.0, const bool isvert=false) :
            Particle<ndim>(xx, aa, wgt, isvert) {
            if (ndim > 1) {
                this->registerVectorField("velocity");
            }
            else {
                this->registerScalarField("velocity");
            }
            this->registerScalarField("potvort");
            this->registerScalarField("relvort");
        }
};

inline scalar_type coriolis(const Vec<2>& pos, const scalar_type f0=0.0, const scalar_type beta=0.0) {return f0 + beta*pos.x[1];}

inline scalar_type coriolis(const Vec<3>& pos, const scalar_type Omega=0.0) {return 2.0*Omega*std::sin(pos.latitude());}

}
}
#endif