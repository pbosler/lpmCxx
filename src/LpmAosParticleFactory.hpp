#ifndef _LPM_AOS_PARTICLE_FACTORY_HPP
#define _LPM_AOS_PARTICLE_FACTORY_HPP

#include <memory>
#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmAosParticle.hpp"
#include "LpmAosSWEParticle.hpp"
#include "LpmAosQGParticle.hpp"

namespace Lpm {
namespace Aos {

template <int ndim> class ParticleFactory {
    public:
        virtual std::unique_ptr<Particle<ndim>> createParticle() const = 0;
        virtual std::unique_ptr<Particle<ndim>> createParticle(const Vec<ndim>& pos, const scalar_type wgt = 0.0,
            const bool isvert=false) const = 0;
        virtual std::unique_ptr<Particle<ndim>> createParticle(const Vec<ndim>& xx, const Vec<ndim>& aa, 
        	const scalar_type wgt=0.0, const bool isvert=false) const=0;
};

template <int ndim> class BasicParticleFactory : public ParticleFactory<ndim> {
    public:
         std::unique_ptr<Particle<ndim>> createParticle() const {
            return std::unique_ptr<Particle<ndim>>(new Particle<ndim>());
        }

         std::unique_ptr<Particle<ndim>> createParticle(const Vec<ndim>& pos, const scalar_type wgt=0.0, 
            const bool isvert=false) const {
            return std::unique_ptr<Particle<ndim>>(new Particle<ndim>(pos, wgt, isvert));
        }
        
        std::unique_ptr<Particle<ndim>> createParticle(const Vec<ndim>& xx, const Vec<ndim>& aa, 
            const scalar_type wgt=0.0, const bool isvert=false) const {
        	return std::unique_ptr<Particle<ndim>>(new Particle<ndim>(xx, aa, wgt, isvert));
        }
};

template <int ndim> class BVEParticleFactory : public ParticleFactory<ndim> {
	public:
		std::unique_ptr<Particle<ndim>> createParticle() const {
			return std::unique_ptr<Particle<ndim>>(new BVEParticle<ndim>());
		}
		
		std::unique_ptr<Particle<ndim>> createParticle(const Vec<ndim>& pos, const scalar_type wgt = 0.0, 
			const bool isvert=false) const {
			return std::unique_ptr<Particle<ndim>>(new BVEParticle<ndim>(pos, wgt, isvert));
		}
		
		std::unique_ptr<Particle<ndim>> createParticle(const Vec<ndim>& xx, const Vec<ndim>& aa, const scalar_type wgt = 0.0, 
			const bool isvert=false) const {
			return std::unique_ptr<Particle<ndim>>(new BVEParticle<ndim>(xx, aa, wgt, isvert));
		}
};

template <int ndim> class SWEParticleFactory : public ParticleFactory<ndim> {
    public:
         std::unique_ptr<Particle<ndim>> createParticle() const {
            return std::unique_ptr<Particle<ndim>>(new SWEParticle<ndim>());
        }
		
		 std::unique_ptr<Particle<ndim>> createParticle(const Vec<ndim>& pos, const scalar_type wgt=0.0, 
		 	const bool isvert=false) const {
			return std::unique_ptr<Particle<ndim>>(new SWEParticle<ndim>(pos, wgt,isvert));
		}
		
		std::unique_ptr<Particle<ndim>> createParticle(const Vec<ndim>& xx, const Vec<ndim>& aa, 
		    const scalar_type wgt=0.0, const bool isvert=false) const {
			return std::unique_ptr<Particle<ndim>>(new SWEParticle<ndim>(xx, aa, wgt, isvert));
		}
		
};



}//namespace Aos



}//namespace Lpm
#endif
