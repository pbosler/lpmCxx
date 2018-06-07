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
        ParticleSet(const std::shared_ptr<ParticleFactory<ndim>> factory, const index_type nMax) : 
            _factory(factory), _nMax(nMax), _nActive(0) {}

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

        void insert(const Vec<ndim>& xx, const Vec<ndim>& aa, const scalar_type area=0.0, const scalar_type vol=0.0);
        void insert(const Vec<ndim>& xx, const scalar_type area = 0.0, const scalar_type vol = 0.0);

        virtual std::string infoString() const;

        std::vector<std::string> particlesInfoStrings() const;

        inline Particle<ndim>* getPtr(const index_type ind) const {return _particles[ind].get();}

        void initScalarFieldFromFn(const std::string& field_name, const AnalyticFunction* fn);
        void initVectorFieldFromFn(const std::string& field_name, const AnalyticFunction* fn);

        inline Vec<ndim> physCrd(const index_type ind) const {return _particles[ind]->physCrd();}
        inline Vec<ndim> lagCrd(const index_type ind) const {return _particles[ind]->lagCrd();}

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
