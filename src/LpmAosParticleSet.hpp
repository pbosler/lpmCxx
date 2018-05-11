#ifndef _LPM_AOS_PARTICLE_SET_HPP
#define _LPM_AOS_PARTICLE_SET_HPP

#include <vector>
#include <memory>
#include <map>
#include <exception>
#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmAosParticleFactory.hpp"

namespace Lpm {

template <int ndim> class ParticleSet {
    public:
        ParticleSet(const std::shared_ptr<ParticleFactory<ndim>> factory) : _factory(factory) {}

        index_type nMax() const {return this->_nMax;}
        index_type n() const {return this->_particles.size();}
        index_type nActive() const {return this->_nActive;}
        index_type nPassive() const {return this->n() - this->_nActive;}
        
        void initFromParticleSetFile(const std::string& fname);
        
    protected:
        index_type _nMax;
        index_type _nActive;
        std::shared_ptr<ParticleFactory<ndim>> _factory;
        std::vector<std::unique_ptr<Particle<ndim>>> _particles;
};

// class CubedSphereParticles : public ParticleSet<3> {
//     public:
//         
//     protected:
// };
// 
// class IcosTriSphereParticles : public ParticleSet<3> {
//     public:
//         
//     protected:
// };
// 
// class CubedSphereCubicParticles : public ParticleSet<3> {
//     public:
//         
//     protected:
// };
// 
// class PlanarTriParticles : public ParticleSet<2> {
//     public:
//         
//     protected:
// };
// 
// class PlanarQuadParticles : public ParticleSet<2> {
//     public:
//         
//     protected:
// };
// 
// class PlanarCubicQuadParticles :: public ParticleSet<2> {
//     public:
//         
//     protected:
// };


}
#endif