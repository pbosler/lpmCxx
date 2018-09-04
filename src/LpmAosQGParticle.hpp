#ifndef LPM_AOS_QG_PARTICLE_HPP
#define LPM_AOS_QG_PARTICLE_HPP

#include "LpmAosTypes.hpp"
#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmAosParticle.hpp"

namespace Lpm {
namespace Aos {

template <int ndim> class QGParticle : public Particle<ndim> {
    public:
        QGParticle() : Particle<ndim>() {
            this->registerVectorField("velocity");
            this->registerScalarField("potvort");
            this->registerScalarField("relvort");
            this->registerScalarField("streamfn");
            this->_wgt_name="area";
        }
        
        QGParticle(const Vec<ndim>& pos, const scalar_type ar=0.0) : Particle<ndim>(pos, ar, "area") {}
        
        QGParticle(const Vec<ndim>& xx, const Vec<ndim>& aa, const scalar_type ar=0.0) : 
            Particle<ndim>(xx, aa, ar, "area"){}
        
        void init(const Vec<ndim>& xx, const Vec<ndim>& aa, const scalar_type wgt=0.0, const std::string& wname="area") override {
            this->_physCrd=xx;
            this->_lagCrd=aa;
            this->_weight=wgt;
            this->_wgt_name=wname;
        }
        
        
        

};

template <int ndim> inline scalar_type coriolis(const Vec<ndim>& pos, const scalar_type f0 = 0.0, const scalar_type beta = 0.0) {
    return f0 + beta*pos.x[1];
}

}
}
#endif
