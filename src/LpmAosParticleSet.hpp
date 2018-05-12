#ifndef _LPM_AOS_PARTICLE_SET_HPP
#define _LPM_AOS_PARTICLE_SET_HPP

#include <vector>
#include <memory>
#include <string>
#include <exception>
#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmAosParticleFactory.hpp"
#include "LpmAnalyticFunctions.h"

namespace Lpm {

template <int ndim> class ParticleSet {
    public:
        ParticleSet(const std::shared_ptr<ParticleFactory<ndim>> factory) : _factory(factory) {}

        virtual ~ParticleSet() {}

        index_type nMax() const {return this->_nMax;}
        index_type n() const {return this->_particles.size();}
        index_type nActive() const {return this->_nActive;}
        index_type nPassive() const {return this->n() - this->_nActive;}

        scalar_type totalVolume() const;
        scalar_type totalArea() const;

        std::vector<std::string> fieldNames() const;

        virtual scalar_type scalarIntegral(const std::string& field_name) const;

        void initFromParticleSetFile(const std::string& fname);

        virtual std::string infoString() const;

        std::vector<std::string> particlesInfoStrings() const;

        void initScalarFieldFromFn(const std::string& field_name, const AnalyticFunction* fn);
        void initVectorFieldFromFn(const std::string& field_name, const AnalyticFunction* fn);

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
