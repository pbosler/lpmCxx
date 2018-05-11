#ifndef _LPM_AOS_PARTICLE_FACTORY_HPP
#define _LPM_AOS_PARTICLE_FACTORY_HPP

#include <memory>
#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmAosParticle.hpp"
#include "LpmAosSWEParticle.hpp"

namespace Lpm {

template <int ndim> class ParticleFactory {
    public:
        virtual std::unique_ptr<Particle<ndim>> createParticle() const = 0;
    protected:
};

template <int ndim> class BasicParticleFactory : public ParticleFactory<ndim> {
    public:
        std::unique_ptr<Particle<ndim>> createParticle() const {
            return std::unique_ptr<Particle<ndim>>(new Particle<ndim>());
        }
    protected:
};

template <int ndim> class SWEParticleFactory : public ParticleFactory<ndim> {
    public:
        std::unique_ptr<Particle<ndim>> createParticle() const {
            return std::unique_ptr<Particle<ndim>>(new SWEParticle<ndim>());
        }
    protected:
};

}

#endif