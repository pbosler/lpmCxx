#include <iostream>
#include <string>
#include <sstream>
#include <memory>
#include <exception>
#include <cmath>
#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmAosTypes.hpp"
#include "LpmAosPolyMesh2d.hpp"
#include "LpmAosParticleFactory.hpp"
#include "LpmAosBVEParticle.hpp"
#include "LpmAosPolyMesh2d.hpp"
#include "LpmAosFaceFactory.hpp"
#include "LpmAosParticleSet.hpp"
#include "LpmAnalyticFunctions.h"

#ifdef HAVE_VTK
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPolyDataWriter.h"
#endif

#ifdef HAVE_KOKKOS
#include "Kokkos_Core.hpp"
#include "Kokkos_View.hpp"
#endif


using namespace Lpm::Aos;
using Lpm::index_type;
using Lpm::scalar_type;
using Lpm::PI;
using Lpm::AnalyticFunction;

#ifdef HAVE_KOKKOS
    typedef Kokkos::View<scalar_type*[3]> vec_view_type;
    typedef typename vec_view_type::HostMirror vec_host_view_type;
    typedef Kokkos::View<scalar_type*> scalar_view_type;
    typedef typename scalar_view_type::HostMirror scalar_host_view_type;
#endif

struct Input {
    std::string prog_name;
    int face_type;
    int init_nest;
    int max_nest;
    int remesh_interval;
    bool use_amr;
    int amr_limit;
    scalar_type dt;
    scalar_type tfinal;
    std::string output_fname;
    int write_interval;
    
    Input(int argc, char* argv[]);
    std::string infoString() const;
    std::shared_ptr<MeshSeed<3>> getSeed() const;
    std::shared_ptr<FaceFactory<3>> getFaceFactory() const;
    std::string make_vtk_fname(const std::string& froot, const int fc) const;
};

struct BVEVelocityDirectSumActive {
    vec_view_type crds;
    scalar_view_type relvort;
    scalar_view_type area;
    vec_view_type velocity;
    
    BVEVelocityDirectSumActive(vec_view_type av, scalar_view_type zeta, 
        scalar_view_type ar, vec_view_type u) :
        crds(av), passive_crds(pv), relvort(zeta), area(ar), velocity(u) {}
        
    KOKKOS_INLINE_FUNCTION
    void operator() (const index_type& i) {
        for (int k=0; k<3; ++k) {
            velocity(i,k) = 0.0;
        }
        for (index_type j=0; j<i; ++j) {
            scalar_type dotprod 0.0;
            for (int k=0; k<3; ++k) {
                dotprod += crds(i,k) * crds(j,k);
            }
            const scalar_type ustr =  -relvort(j) * area(j) / (4.0*PI * (1.0 - dotprod));
            velocity(i, 0) += (crds(i,1)*crds(j,2) - crds(i,2)*crds(j,1))*ustr;
            velocity(i, 1) += (crds(i,2)*crds(j,0) - crds(i,0)*crds(j,2))*ustr;
            velocity(i, 2) += (crds(i,0)*crds(j,1) - crds(i,1)*crds(j,0))*ustr;
        }
        for (index_type j=i+1; j<crds.dimension_0(); ++j) {
            scalar_type dotprod 0.0;
            for (int k=0; k<3; ++k) {
                dotprod += crds(i,k) * crds(j,k);
            }
            const scalar_type ustr =  -relvort(j) * area(j) / (4.0*PI * (1.0 - dotprod));
            velocity(i, 0) += (crds(i,1)*crds(j,2) - crds(i,2)*crds(j,1))*ustr;
            velocity(i, 1) += (crds(i,2)*crds(j,0) - crds(i,0)*crds(j,2))*ustr;
            velocity(i, 2) += (crds(i,0)*crds(j,1) - crds(i,1)*crds(j,0))*ustr;
        }
    }
};

struct BVEVelocityDirectSumPassive {
    vec_view_type acrds; // active
    vec_view_type pcrds; // passive
    scalar_view_type relvort; // active
    scalar_view_type area; // active
    vec_view_type velocity; // passive particles
    
    BVEVelocityDirectSumPassive(vec_view_type acrds_, vec_view_type pcrds_, scalar_view_type zeta, 
        scalar_view_type ar, vec_view_type u) :
        acrds(acrds_), pcrds(pcrds_), relvort(zeta), area(ar), velocity(u) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()  (const index_type& i) {
        for (int k=0; k<3; ++k) {
            velocity(i,k) = 0.0;
        }
        for (index_type j=0; j<acrds.dimension_0(); ++j) {
            scalar_type dotprod = 0.0;
            for (int k=0; k<3; ++k) {
                dotprod += pcrds(i,k)*acrds(j,k);
            }
            const scalar_type ustr = -relvort(j)*area(j) / (4.0*PI*(1.0-dotprod));
            velocity(i,0) += (pcrds(i,1)*acrds(j,2) - pcrds(i,2)*acrds(j,1))*ustr;
            velocity(i,1) += (pcrds(i,2)*acrds(j,0) - pcrds(i,0)*acrds(j,2))*ustr;
            velocity(i,2) += (pcrds(i,0)*acrds(j,1) - pcrds(i,1)*acrds(j,0))*ustr;
        }
    }
};

struct BVEStreamFnDirectSumActive {
    vec_view_type crds;
    scalar_view_type relvort;
    scalar_view_type area;
    scalar_view_type stream;
    
    BVEStreamFnDirectSumActive(vec_view_type c, scalar_view_type zeta, scalar_view_type ar, 
        scalar_view_type psi) : crds(c), relvort(zeta), area(ar), stream(psi) {}
    
    void operator() (const index_type& i) {
        stream(i) = 0.0;
        for (index_type j=0; j<i; ++j) {
            scalar_type dotprod = 0.0;
            for (int k=0; k<3; ++k) {
                dotprod += crds(i,k) * crds(j,k);
            }
            stream(i) += std::log(1.0-dotprod)*relvort(j)*area(j)/(4.0*PI);
        }
        for (index_type j=i+1; j<crds.dimension_0(); ++j) {
            scalar_type dotprod = 0.0;
            for (int k=0; k<3; ++k) {
                dotprod += crds(i,k) * crds(j,k);
            }
            stream(i) -= std::log(1.0-dotprod)*relvort(j)*area(j)/(4.0*PI);
        }
    }
};

struct BVEStreamFnDirectSumPassive {
    vec_view_type acrds;
    vec_view_type pcrds;
    scalar_view_type relvort;
    scalar_view_type area;
    scalar_view_type stream;
    
    BVEStreamFnDirectSumPassive(vec_view_type ac, vec_view_type pc, scalar_view_type zeta, scalar_view_type ar,
        scalar_view_type psi) : acrds(ac), pcrds(pc), relvort(zeta), area(ar), stream(psi) {}
    
    void operator() (const index_type& i) {
        stream(i) = 0.0;
        for (index_type j=0; j<acrds.dimension_0(); ++j) {
            scalar_type dotprod = 0.0;
            for (int k=0; k<3; ++k) {
                dotprod += pcrds(i,k)*acrds(j,k);
            }
            stream(i) -= std::log(1.0-dotprod)*relvort(j)*area(j)/(4.0*PI);
        }
    }
};


int main(int argc, char* argv[]) {
#ifdef HAVE_KOKKOS
    Kokkos::initialize(argc, argv);
    {
    //
    //  problem setup, memory allocation
    //
    const scalar_type SphRadius = 1.0;
    const scalar_type Omega = 2.0 * PI / SphRadius;  // sphere's angular rotation rate about z-axis
    Input input(argc, argv);
    std::cout << input.infoString();
    
    std::shared_ptr<MeshSeed<3>> seed_ptr = input.getSeed();
    std::shared_ptr<EdgeFactory<3>> efac(new LinearEdgeFactory<3>());
    std::shared_ptr<FaceFactory<3>> ffac = input.getFaceFactory();
    std::shared_ptr<ParticleFactory<3>> pfac(new BVEParticleFactory<3>());
    
    PolyMesh2d<3> sphere(seed_ptr, pfac, efac, ffac, input.init_nest, input.max_nest, input.amr_limit, SphRadius);
    sphere.initStaggeredVerticesAndFacesFromSeed();
    std::cout << sphere.infoString();
    
    //
    //  initial conditions
    //
    const Lpm::solidBodyVorticity vorfn(Omega);
    const Lpm::solidBodyVelocity velfn(Omega);
    ParticleSet<3>* pptr = sphere.particle_set_raw_ptr();
    
    pptr->registerScalarField("stream_error");
    pptr->registerScalarField("relvort_error");
    pptr->registerScalarField("position_error");
    pptr->registerScalarField("velocity_error");
    
    for (index_type i=0; i<pptr->n(); ++i) {
        const Vec<3> pc = pptr->physCrd(i);
        pptr->initScalarFieldValue("relvort", i, vorfn.evaluateScalar(pc));
        pptr->initVectorFieldValue("velocity", i, velfn.evaluateVector(pc));
        pptr->initScalarFieldValue("stream_fn", i, vorfn.evaluateScalar(pc));
    }
    
    
    
#ifdef HAVE_VTK
    std::ostringstream vtkss;
    int frame_counter = 0;
    vtkss << input.output_fname << "_" <<  (input.face_type == 3 ? "icos" : "cubs") << input.init_nest << "_to_"
         << input.max_nest << "_dt" << input.dt << "_";
    const std::string vtk_froot = vtkss.str();
    {
        const std::string ofname = input.make_vtk_fname(vtk_froot, frame_counter++);
        vtkSmartPointer<vtkPolyData> pd = sphere.toVtkPolyData();
        vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
        writer->SetInputData(pd);
        writer->SetFileName(ofname.c_str());
        writer->Write();
        std::cout << "wrote output file: " << ofname << std::endl;
    
    }
#endif    
    
    //
    //  time integration
    //
    const index_type ntimesteps = std::floor(input.tfinal/input.dt);
    const scalar_type dt = input.tfinal/scalar_type(ntimesteps);
    
    vec_view_type active_view;
    vec_host_view_type active_host;
    vec_view_type passive_view;
    vec_host_view_type passive_host;
    vec_view_type vel_view;
    vec_host_view_type vel_host_view;
    scalar_view_type vort_view;
    scalar_host_view_type vort_host_view;
    scalar_view_type area_view;
    scalar_host_view_type area_host;
    
    pptr->init_pack_active_passive_coords(active_view, active_host, passive_view, 
        passive_host, area_view, area_host);
    pptr->init_pack_scalar_field(vort_view, vort_host_view, "relvort");
    pptr->init_pack_vector_field(vel_view, vel_host_view, "velocity");
    
    }
    Kokkos::finalize();
    return 0;
#else
    std::cout << "BVESolidBody requires Kokkos." << std::endl;
    return 1;
#endif    
}

Input::Input(int argc, char* argv[]) {
    prog_name = argv[0];
    face_type = 3;
    max_nest = 2;
    init_nest = 2;
    amr_limit = 0;
    dt = 0.1;
    tfinal = 1.0;
    remesh_interval = -1;
    write_interval = 1;
    output_fname = "./bve/solid_body";    
    for (int i=1; i<argc; ++i) {
        const std::string& token = argv[i];
        if (token == "-i" || token == "--init_nest") {
            init_nest = std::stoi(argv[++i]);
        }
        else if (token == "-o") {
            output_fname = argv[++i];
        }
        else if (token == "-f") {
            face_type = std::stoi(argv[++i]);
        }
        else if (token == "-m" || token == "--max_nest") {
            max_nest = std::stoi(argv[++i]);
        }
        else if (token == "--amr_limit") {
            amr_limit = std::stoi(argv[++i]);
        }
        else if (token == "-dt") {
            dt = std::stod(argv[++i]);
        }
        else if (token == "-tf") {
            tfinal = std::stod(argv[++i]);
        }
        else if (token == "-wi") {
            write_interval = std::stoi(argv[++i]);
        }
        else if (token == "-rm") {
            remesh_interval = std::stoi(argv[++i]);
        }
    }
    use_amr = (amr_limit > 0 && max_nest > init_nest);
    if (face_type != 3 && face_type != 4) {
        throw std::runtime_error("BVE Solid Body ERROR: invalid face type.");
    }
    if (max_nest < init_nest) {
        max_nest = init_nest;
        std::cout << "BVE Solid Body WARNING: max_nest must be >= init_nest." << std::endl;
        std::cout << "... setting max_nest = " << init_nest << std::endl;
    }
    if (amr_limit > 0 && max_nest == init_nest) {
        throw std::logic_error("BVE Solid Body ERROR: to use amr, max_nest must be larger than init_nest.");
    }
    if (write_interval < remesh_interval) {
        std::cout << "I/O Performance note: I/O can only be done from the host, not device.  " 
            << "You have requested I/O between remeshing steps, which will decrease performance.";
    }
}

std::shared_ptr<MeshSeed<3>> Input::getSeed() const {
    if (face_type == 3) {
        return std::shared_ptr<MeshSeed<3>>(new IcosTriSphereSeed());
    }
    else if (face_type == 4) {
        return std::shared_ptr<MeshSeed<3>>(new CubedSphereSeed());
    }
}

std::shared_ptr<FaceFactory<3>> Input::getFaceFactory() const {
    if (face_type == 3) {
        return std::shared_ptr<FaceFactory<3>>(new TriFaceFactory<3>());
    }
    else if (face_type == 4) {
        return std::shared_ptr<FaceFactory<3>>(new QuadFaceFactory<3>());
    }
}

std::string Input::make_vtk_fname(const std::string& fr, const int fc) const {
    std::ostringstream ss;
    ss << fr << std::setw(4) << std::setfill('0') << fc << ".vtk";
    return ss.str();
}

std::string Input::infoString() const {
    std::ostringstream ss;
    ss << prog_name << " input summary:" << std::endl;
    ss << "face_type: " << (face_type == 3 ? "triangular" : "quadrilateral") << std::endl;
    ss << "max_nest = " << max_nest << std::endl;
    ss << "init_nest = " << init_nest << std::endl;
    ss << "amr_limit = " << amr_limit << std::endl;
    ss << "use_amr?: " << (use_amr ? "yes" : "no") << std::endl;
    ss << "dt = " << dt << std::endl;
    ss << "tfinal = " << tfinal << std::endl;
    ss << "remesh_interval = " << remesh_interval << std::endl;
    ss << "write_interval = " << write_interval << std::endl;
    ss << "output_fname = " << output_fname << std::endl;
    return ss.str();
}