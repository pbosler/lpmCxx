#ifndef LPM_AOS_SWE_PARTICLE_HPP
#define LPM_AOS_SWE_PARTICLE_HPP

#include "LpmAosTypes.hpp"
#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmAosParticle.hpp"

namespace Lpm {

template <int ndim=3> struct SWEParticle : public Particle<ndim> {
    scalar_type velocity[ndim];
    scalar_type potvort;
    scalar_type relvort;
    scalar_type depth;
    scalar_type surf_height;
    scalar_type bottom_height;
    scalar_type divergence;

    SWEParticle() : Particle(), area(0.0), volume(0.0), potvort(0.0), relvort(0.0), depth(0.0), surf_height(0.0),
        bottom_height(0.0), divergence(0.0) {
        for (int i=0; i<ndim; ++i)
            velocity[i] = 0.0;
    };

    SWEParticle(const Vec<ndim> initCrd, const scalar_type area=0.0, const scalar_type volume=0.0,
        const scalar_type pv=0.0, const scalar_type rv=0.0, const scalar_type hh=0.0,
        const scalar_type sh=0.0, const scalar_type bh=0.0, const scalar_type div=0.0) :
        Particle(initCrd, area, volume), potvort(pv), relvort(rv), depth(hh), surf_height(sh), bottom_height(bh),
        divergence(div) {
            for (int i=0; i<ndim; ++i)
            velocity[i] = 0.0;
        }

    SWEParticle(const Vec<ndim> xx, const Vec<ndim> aa, const scalar_type area=0.0, const scalar_type volume=0.0,
        const scalar_type pv=0.0, const scalar_type rv=0.0, const scalar_type hh=0.0,
        const scalar_type sh=0.0, const scalar_type bh=0.0, const scalar_type div=0.0) :
        Particle(xx, aa, area, volume), potvort(pv), relvort(rv), depth(hh), surf_height(sh), bottom_height(bh),
        divergence(div) {
            for (int i=0; i<ndim; ++i)
            velocity[i] = 0.0;
        }

    scalar_type calcPV(const scalar_type f) const {return (this->relvort + f)/this->depth;}

    void setDepthFromPV(const scalar_type f) {this->depth = (this->relvort + f)/this->potvort}

    void setRelvortFromPV(const scalar_type f) {this->relvort = this->potvort * this->depth - f;}

    void setVel(const scalar_type* vv) {
        for (int i=0; i<ndim; ++i)
            this->velocity[i] = vv[i];
    }

    void setVel(const std::vector<scalar_type>& vv) {
        for (int i=0; i<ndim; ++i)
            this->velocity[i] = vv[i];
    }

    void setVel(const scalar_type u, const scalar_type v) {
        this->velocity[0] = u;
        this->velocity[1] = v;
    }

    void setVel(const scalar_type u, const scalar_type v, const scalar_type w) {
        this->velocity[0] = u;
        this->velocity[1] = v;
        this->velocity[2] = w;
    }
}

scalar_type coriolis(const Vec<2>& pos, const scalar_type f0=0.0, const scalar_type beta=0.0) {return f0 + beta*pos.x[1];}

scalar_type coriolis(const Vec<3>& pos, const scalar_type Omega=0.0) {return 2.0*Omega*pos.latitude();}




}
#endif
