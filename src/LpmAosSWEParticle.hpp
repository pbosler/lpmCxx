#ifndef LPM_AOS_SWE_PARTICLE_HPP
#define LPM_AOS_SWE_PARTICLE_HPP

#include "LpmAosTypes.hpp"
#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmAosParticle.hpp"
#include "LpmAosBVEParticle.hpp"
#include <cmath>

namespace Lpm {
namespace Aos {

template <int ndim=3> class SWEParticle : public BVEParticle<ndim> {
    public:
        SWEParticle(const bool isvert=false) : 
        	BVEParticle<ndim>(isvert) {
            this->registerScalarField("divergence");
            this->registerScalarField("depth");
            this->registerScalarField("surf_height");
            this->registerScalarField("bottom_height");
        }

        SWEParticle(const Vec<ndim>& pos, const scalar_type wgt = 0.0, const bool isvert=false) : 
        	BVEParticle<ndim>(pos, wgt, isvert) {
            this->registerScalarField("divergence");
            this->registerScalarField("depth");
            this->registerScalarField("surf_height");
            this->registerScalarField("bottom_height");
        }
        
        SWEParticle(const Vec<ndim>& xx, const Vec<ndim>& aa, const scalar_type wgt = 0.0, const bool isvert=false) :
            BVEParticle<ndim>(xx, aa, wgt, isvert) {
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

}
}
#endif
