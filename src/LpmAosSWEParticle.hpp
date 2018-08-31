#ifndef LPM_AOS_SWE_PARTICLE_HPP
#define LPM_AOS_SWE_PARTICLE_HPP

#include "LpmAosTypes.hpp"
#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmAosParticle.hpp"

namespace Lpm {

template <int ndim=3> class SWEParticle : public Particle<ndim> {
    public:
        SWEParticle() : Particle<ndim>() {
            if (ndim > 1) {
                this->registerVectorField("velocity");
            }
            else {
                this->registerScalarField("velocity");
            }
            this->registerScalarField("potvort");
            this->registerScalarField("relvort");
            this->registerScalarField("divergence");
            this->registerScalarField("depth");
            this->registerScalarField("surf_height");
            this->registerScalarField("bottom_height");
        }

        SWEParticle(const Vec<ndim>& pos, const scalar_type len = 0.0, const scalar_type ar=0.0, const scalar_type vv=0.0) : Particle<ndim>(pos, len, ar, vv) {
            if (ndim > 1) {
                this->registerVectorField("velocity");
            }
            else {
                this->registerScalarField("velocity");
            }
            this->registerScalarField("potvort");
            this->registerScalarField("relvort");
            this->registerScalarField("divergence");
            this->registerScalarField("depth");
            this->registerScalarField("surf_height");
            this->registerScalarField("bottom_height");
        }

        SWEParticle(const Vec<ndim>& xx, const Vec<ndim>& aa, const scalar_type len = 0.0, const scalar_type ar = 0.0, const scalar_type vv=0.0) :
            Particle<ndim>(xx, aa, len, ar, vv) {
            if (ndim > 1) {
                this->registerVectorField("velocity");
            }
            else {
                this->registerScalarField("velocity");
            }
            this->registerScalarField("potvort");
            this->registerScalarField("relvort");
            this->registerScalarField("divergence");
            this->registerScalarField("depth");
            this->registerScalarField("surf_height");
            this->registerScalarField("bottom_height");
        }

        scalar_type calcPV(const scalar_type f) const {return (this->relvort + f)/this->depth;}

        void setDepthFromPV(const scalar_type f) {
            this->depth = (this->_sFields["relvort"] + f)/this->_sFields["potvort"];
            this->surf_height = this->_sFields["bottom_height"] + this->_sFields["depth"];
        }

        void setRelvortFromPV(const scalar_type f) {this->_sFields["relvort"] =
            this->_sFields["potvort"] * this->_sFields["depth"] - f;
        }

};

inline scalar_type coriolis(const Vec<2>& pos, const scalar_type f0=0.0, const scalar_type beta=0.0) {return f0 + beta*pos.x[1];}

inline scalar_type coriolis(const Vec<3>& pos, const scalar_type Omega=0.0) {return 2.0*Omega*pos.latitude();}




}
#endif
