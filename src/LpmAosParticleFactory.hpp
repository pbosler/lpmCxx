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
        virtual std::unique_ptr<Particle<ndim>> createParticle(const Vec<ndim>& pos, const scalar_type aa=0.0,
            const scalar_type vv=0.0) const = 0;
        virtual std::unique_ptr<Particle<ndim>> createParticle(const Vec<ndim>& xx, const Vec<ndim>& aa,
            const scalar_type ar=0.0, const scalar_type vol=0.0) const = 0;
    protected:
};

template <int ndim> class BasicParticleFactory : public ParticleFactory<ndim> {
    public:
        std::unique_ptr<Particle<ndim>> createParticle() const {
            return std::unique_ptr<Particle<ndim>>(new Particle<ndim>());
        }

        std::unique_ptr<Particle<ndim>> createParticle(const Vec<ndim>& pos, const scalar_type aa=0.0,
            const scalar_type vv=0.0) const {
            return std::unique_ptr<Particle<ndim>>(new Particle<ndim>(pos, aa, vv));
        }

        std::unique_ptr<Particle<ndim>> createParticle(const Vec<ndim>& xx, const Vec<ndim>& aa,
            const scalar_type ar=0.0, const scalar_type vol=0.0) const {
            return std::unique_ptr<Particle<ndim>>(new Particle<ndim>(xx, aa, ar, vol));
        }
    protected:
};

template <int ndim> class SWEParticleFactory : public ParticleFactory<ndim> {
    public:
        std::unique_ptr<Particle<ndim>> createParticle() const {
            return std::unique_ptr<Particle<ndim>>(new SWEParticle<ndim>());
        }

        std::unique_ptr<Particle<ndim>> createParticle(const Vec<ndim>& pos, const scalar_type aa=0.0,
            const scalar_type vv=0.0) const {
            return std::unique_ptr<Particle<ndim>>(new SWEParticle<ndim>(pos, aa, 0.0));
        }

        std::unique_ptr<Particle<ndim>> createParticle(const Vec<ndim>& xx, const Vec<ndim>& aa,
            const scalar_type ar=0.0, const scalar_type vol=0.0) const {
            return std::unique_ptr<Particle<ndim>>(new SWEParticle<ndim>(xx, aa, ar, vol));
        }
    protected:
};

}

#endif
