#ifndef LPM_AOS_TIME_INT_HPP
#define LPM_AOS_TIME_INT_HPP

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmUtilities.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <vector>
#include "LpmAosPolyMesh2d.hpp"

#ifdef HAVE_KOKKOS
#include "Kokkos_Core.hpp"
#include "Kokkos_View.hpp"


namespace Lpm {
namespace Aos {

typedef Kokkos::TeamPolicy<> team_policy;
typedef typename team_policy::member_type member_type;
typedef Kokkos::View<scalar_type*[3]> view_3d;
typedef Kokkos::View<scalar_type*[2]> view_2d;
typedef Kokkos::View<scalar_type*> view_1d;
typedef Kokkos::View<const scalar_type*[3]> const_view_3d;
typedef Kokkos::View<const scalar_type*[2]> const_view_2d;
typedef Kokkos::View<const scalar_type*> const_view_1d;
typedef typename view_3d::HostMirror host_view3d;
typedef typename view_2d::HostMirror host_view2d;
typedef typename view_1d::HostMirror host_view1d;

struct SphereBiotSavartKernel {
    const_view_3d src_crds;
    const_view_1d src_vort;
    const_view_1d src_area;
    typedef scalar_type value_type[];
    value_type tgt_crd;
    typedef view_3d::size_type size_type;
    static constexpr size_type value_count = 3;
    static constexpr scalar_type OO4PI = 1.0 / (4.0 * PI);
    static constexpr scalar_type FLOAT_ZERO = 1.0e-14;
    
    KOKKOS_INLINE_FUNCTION
    SphereBiotSavartKernel(value_type tgt_crd, const_view_3d sc, const_view_1d zeta, const_view_1d a) :
        src_crds(sc), src_vort(zeta), src_area(a) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator() (const index_type j, value_type sum) const {
        const scalar_type dotprod = tgt_crd[0]*src_crds(j,0) + tgt_crd[1]*src_crds(j,1) + tgt_crd[2]*src_crds(j,2);
        const scalar_type denom = 1.0 - dotprod;
        if (std::abs(denom) > FLOAT_ZERO) {
            const scalar_type str = src_vort(j) * src_area(j) * OO4PI / denom;
            sum[0] += (tgt_crd[1] * src_crds(j,2) - tgt_crd[2] * src_crds(j,1)) * str;
            sum[1] += (tgt_crd[2] * src_crds(j,0) - tgt_crd[0] * src_crds(j,2)) * str;
            sum[2] += (tgt_crd[0] * src_crds(j,1) - tgt_crd[1] * src_crds(j,0)) * str;
        }
    }
    
    KOKKOS_INLINE_FUNCTION
    void join(volatile value_type dst, const volatile value_type src) const {
        for (int j=0; j<3; ++j)
            dst[j] += src[j];
    }
    
    KOKKOS_INLINE_FUNCTION
    void init(value_type sum) {
        for (int j=0; j<3; ++j) 
            sum[j] = 0.0;
    }
};

template <typename InnerReduce> struct Outer3dVectorSum {
    view_3d result;
    const_view_3d tgt_crds;
    const_view_3d src_crds;
    const_view_1d f_convol;
    const_view_1d src_wgt;
    
    Outer3dVectorSum(view_3d uvw, const_view_3d t, const_view_3d s, const_view_1d ff, const_view_1d wgt) : 
        result(uvw), tgt_crds(t), src_crds(s), f_convol(ff), src_wgt(wgt) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator() (const member_type& team_member) const {
        const index_type tgt_ind = team_member.league_rank();
        const scalar_type tgt_crd[3] = {tgt_crds(tgt_ind,0), tgt_crds(tgt_ind, 1), tgt_crds(tgt_ind,2)};
        scalar_type uvw[3];
        Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team_member, src_crds.extent(0)), 
            InnerReduce(tgt_crd, src_crds, f_convol, src_wgt), uvw);
        for (int j=0; j<3; ++j) 
            result(tgt_ind, j) = uvw[j];
    }
};


#endif // HAVE_KOKKOS
}
}
#endif
