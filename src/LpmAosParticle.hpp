#ifndef _LPM_AOS_PARTICLE_HPP
#define _LPM_AOS_PARTICLE_HPP

#include "LpmAosTypes.hpp"
#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#ifdef USE_KOKKOS
#include "Kokkos_Core.hpp"
#endif
#include <vector>
#include <memory>

namespace Lpm {

template <int ndim=3> struct Particle {
    Vec<ndim> physCrd;
    Vec<ndim> lagCrd;
    scalar_type area;
    scalar_type volume;

    Particle() : physCrd(), lagCrd(), area(0.0), volume(0.0) {};

    Particle(const Vec<ndim>& initCrd, const scalar_type aa=0.0, const scalar_type vv=0.0) : physCrd(initCrd),
        lagCrd(initCrd), area(aa), volume(vv) {};

    Particle(const Vec<ndim>& xx, const Vec<ndim>& aa, const scalar_type ar=0.0, const scalar_type vv=0.0) :
        physCrd(xx), lagCrd(aa), area(ar), volume(vv) {};
};

enum PhysicsType {BASIC, SHALLOW_WATER};

}

#endif
