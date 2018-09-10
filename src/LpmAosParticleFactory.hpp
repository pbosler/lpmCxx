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
        	const std::string wname="area", const bool isvert=false) const = 0;
        virtual std::unique_ptr<Particle<ndim>> createParticle(const Vec<ndim>& xx, const Vec<ndim>& aa, 
        	const scalar_type wgt=0.0, const bool isvert=false) const;
        virtual std::unique_ptr<Particle<ndim>> createParticle(const Vec<ndim>& xx, const Vec<ndim>& aa,
            const scalar_type wgt = 0.0, const std::string wname="area", const bool isvert=false) const = 0;
        virtual std::unique_ptr<Particle<ndim>> createParticle(const Vec<ndim>& xx, const scalar_type wgt=0.0, const bool isvert=false) const=0;
};

template <int ndim> class BasicParticleFactory : public ParticleFactory<ndim> {
    public:
         std::unique_ptr<Particle<ndim>> createParticle() const {
            return std::unique_ptr<Particle<ndim>>(new Particle<ndim>("area"));
        }

         std::unique_ptr<Particle<ndim>> createParticle(const Vec<ndim>& pos, const scalar_type wgt=0.0, const std::string wname="area", const bool isvert=false) const {
            return std::unique_ptr<Particle<ndim>>(new Particle<ndim>(pos, wgt, wname, isvert));
        }
        
        std::unique_ptr<Particle<ndim>> createParticle(const Vec<ndim>& xx, const Vec<ndim>& aa, const scalar_type wgt=0.0, const bool isvert=false) const {
        	return std::unique_ptr<Particle<ndim>>(new Particle<ndim>(xx, aa, wgt, "area", isvert));
        }

         std::unique_ptr<Particle<ndim>> createParticle(const Vec<ndim>& xx, const Vec<ndim>& aa, const scalar_type wgt=0.0, const std::string wname="area"
         	, const bool isvert=false) const {
            return std::unique_ptr<Particle<ndim>>(new Particle<ndim>(xx, aa, wgt, wname, isvert));
        }
        
        std::unique_ptr<Particle<ndim>> createParticle(const Vec<ndim>& xx, const scalar_type wgt=0.0, const bool isvert=false) const {
        	return std::unique_ptr<Particle<ndim>>(new Particle<ndim>(xx, wgt, "area", isvert));
        }
};

template <int ndim> class SWEParticleFactory : public ParticleFactory<ndim> {
    public:
         std::unique_ptr<Particle<ndim>> createParticle() const {
            return std::unique_ptr<Particle<ndim>>(new SWEParticle<ndim>((ndim==1 ? "length" : "area")));
        }
		
		 std::unique_ptr<Particle<ndim>> createParticle(const Vec<ndim>& pos, const scalar_type wgt=0.0, const std::string wname=(ndim==1 ? "length" : "area"), 
		 	const bool isvert=false) const {
			return std::unique_ptr<Particle<ndim>>(new SWEParticle<ndim>(pos, wgt, (ndim==1 ? "length" : "area"), isvert));
		}
		
		std::unique_ptr<Particle<ndim>> createParticle(const Vec<ndim>& xx, const Vec<ndim>& aa, const scalar_type wgt=0.0, 
			const bool isvert=false) const {
			return std::unique_ptr<Particle<ndim>>(new SWEParticle<ndim>(xx, aa, wgt, (ndim==1 ? "length" : "area"), isvert));
		}
		
		 std::unique_ptr<Particle<ndim>> createParticle(const Vec<ndim>& pos, const Vec<ndim>& aa, const scalar_type wgt=0.0, 
		 	const std::string wname=(ndim==1 ? "length" : "area"), const bool isvert=false) const {
			return std::unique_ptr<Particle<ndim>>(new SWEParticle<ndim>(pos, aa, wgt, wname, isvert));
		}
		
		std::unique_ptr<Particle<ndim>> createParticle(const Vec<ndim>& xx, const scalar_type wgt=0.0, const bool isvert=false) const {
			return std::unique_ptr<Particle<ndim>>(new SWEParticle<ndim>(xx, wgt, (ndim==1 ? "length" : "area"), isvert));
		}
		
//         std::unique_ptr<Particle<ndim>> createParticle(const Vec<ndim>& pos, const scalar_type wgt =0.0) const {
//             return std::unique_ptr<Particle<ndim>>(new SWEParticle<ndim>(pos, wgt, ));
//         }
// 
//         std::unique_ptr<Particle<ndim>> createParticle(const Vec<ndim>& xx, const Vec<ndim>& aa, const scalar_type wgt=0.0) const {
//             return std::unique_ptr<Particle<ndim>>(new SWEParticle<ndim>(xx, aa, wgt, (ndim==1 ? "length" : "area")));
//         }

};

template <int ndim> class QGParticleFactory : public ParticleFactory<ndim> {
	public:
		 std::unique_ptr<Particle<ndim>> createParticle() const {
			return std::unique_ptr<Particle<ndim>>(new QGParticle<ndim>((ndim==1 ? "length" : "area")));
		}
		
		 std::unique_ptr<Particle<ndim>> createParticle(const Vec<ndim>& pos, const scalar_type wgt=0.0, const std::string wname="area",
		  const bool isvert=false) const {
			return std::unique_ptr<Particle<ndim>>(new QGParticle<ndim>(pos, wgt, wname, isvert));
		}
		
		 std::unique_ptr<Particle<ndim>> createParticle(const Vec<ndim>& xx, const Vec<ndim>& aa, const scalar_type wgt=0.0, 
		 	const std::string wname="area", const bool isvert=false) const {
			return std::unique_ptr<Particle<ndim>>(new QGParticle<ndim>(xx, aa, wgt, wname, false));
		}
		
		std::unique_ptr<Particle<ndim>> createParticle(const Vec<ndim>& xx, const Vec<ndim>& aa, const scalar_type wgt=0.0, 
			const bool isvert =false) const {
			return std::unique_ptr<Particle<ndim>>(new QGParticle<ndim>(xx, aa, wgt, "area", isvert));
		}
		
// 		std::unique_ptr<Particle<ndim>> createParticle(const Vec<ndim>& pos, const scalar_type ar=0.0) const {
// 			return std::unique_ptr<Particle<ndim>>(new QGParticle<ndim>(pos,ar, "area"));
// 		}
// 		
// 		std::unique_ptr<Particle<ndim>> createParticle(const Vec<ndim>& xx, const Vec<ndim>& aa, const scalar_type ar=0.0) const {
// 			return std::unique_ptr<Particle<ndim>>(new QGParticle<ndim>(xx,aa,ar, "area"));
// 		}
};

}
}
#endif
