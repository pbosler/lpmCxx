#include "LpmConfig.h"
#include "LpmAosTypes.hpp"
#include "LpmTypeDefs.h"
#include "LpmMPIReplicatedData.h"
#include <iostream>
#include "Kokkos_Core.hpp"
#include "Kokkos_View.hpp"
#include <mpi.h>
#include <cstdio>
#include <typeinfo>

using Lpm::index_type;
using Lpm::scalar_type;

using namespace Lpm::Aos;

typedef Kokkos::View<index_type*> view_type;
typedef Kokkos::View<const index_type*> const_view_type;
typedef typename view_type::HostMirror host_view_type;
typedef Kokkos::pair<index_type, index_type> pair_type;

typedef Kokkos::TeamPolicy<> team_policy;
typedef typename team_policy::member_type member_type;

struct reduce_functor {
    const_view_type src;
    typedef index_type value_type;
    
    KOKKOS_INLINE_FUNCTION
    reduce_functor(const_view_type s) : src(s) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator() (const index_type j, index_type& sum) const {
        sum += src(j);
    }
    
    KOKKOS_INLINE_FUNCTION
    void join(volatile index_type& dest, const volatile index_type& src) const {
        dest += src;
    }
    
    KOKKOS_INLINE_FUNCTION
    void init(index_type& val) const {
        val = 0;
    }
};

template <typename Inner> struct Outer {
    view_type tgt;
    const_view_type src;
    
    Outer(view_type t, const_view_type s) : tgt(t), src(s) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator() (const member_type& mbr) const {
        const index_type i = mbr.league_rank();
        index_type sum = 0;
        Kokkos::parallel_reduce(Kokkos::TeamThreadRange(mbr, src.extent(0)), Inner(src), sum);
        tgt(i) = sum;
    }
};

int main(int argc, char* argv[]) {
    // initialize mpi
    int mpiErrCode;
    int numProcs = 1;
    int procRank = 0;
    mpiErrCode = MPI_Init(&argc, &argv);
    mpiErrCode = MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    mpiErrCode = MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    // initialize kokkos
    Kokkos::initialize(argc, argv);
    {
    if (procRank == 0) {
        printf("MPI/Kokkos nested parallelism example on ExeSpace = %s\n", 
            typeid(Kokkos::DefaultExecutionSpace).name());
    }
    
    const index_type nn = 10045;
    // create views on device
    view_type sv("src", nn);
    view_type tv("tgt", nn);
    
    // partition work across mpi ranks
    Lpm::MPIReplicatedData mpi(nn, procRank, numProcs);
    //if (procRank == 0) std::cout << mpi.infoString();
    
    // initialize data on host
    host_view_type hsv = Kokkos::create_mirror_view(sv);
    for (index_type i=0; i<nn; ++i) 
        hsv(i) = i;
    // copy data to device
    Kokkos::deep_copy(sv, hsv);
    
    // compute partitioned work on device
    const index_type my_work_start = mpi.startIndex(procRank);
    const index_type my_work_end = mpi.endIndex(procRank) + 1;
    const index_type work_size = mpi.msgSize(procRank);
    view_type work_view = Kokkos::subview(tv, pair_type(my_work_start, my_work_end));
    
    // Launch 1 thread team per rank-local target
    const team_policy policy(work_size, Kokkos::AUTO);
    Kokkos::parallel_for(policy, Outer<reduce_functor>(work_view, sv));
    Kokkos::fence();
    
    
    // copy results from device to host
    host_view_type htv = Kokkos::create_mirror_view(tv);
    host_view_type work_host = Kokkos::subview(htv, pair_type(my_work_start, my_work_end));
    Kokkos::deep_copy(work_host, work_view);
    
    // broadcast results to all mpi ranks
    for (int i=0; i<numProcs; ++i) {
        index_type* bufstart = htv.data() + mpi.startIndex(i);
        MPI_Bcast(bufstart, mpi.msgSize(i), MPI_INT, i, MPI_COMM_WORLD);
    }
    
    // check for errors 
    bool test_pass = true;
    const index_type true_sum = nn*(nn-1)/2;
    for (int k = 0; k<numProcs; ++k) {
        if (procRank == k) {
            for (index_type i=0; i<nn; ++i) {
                if (htv(i) != true_sum) {
                    std::cout << "error at mpi rank " << procRank << ": htv(" << i << ") = " << htv(i) 
                              << " != " << true_sum << std::endl;
                    test_pass = false;
                }
            }
        }
        mpiErrCode = MPI_Barrier(MPI_COMM_WORLD);
    }
    if (procRank==0) {
        std::cout << (test_pass ? "test passed." : "FAILURE detected.") << std::endl;
    }
    }
    Kokkos::finalize();
    MPI_Finalize();
return 0;
}