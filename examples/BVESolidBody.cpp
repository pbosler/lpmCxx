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
    
    sphere.registerScalarField("stream_error");
    sphere.registerScalarField("relvort_error");
    sphere.registerScalarField("position_error");
    sphere.registerScalarField("velocity_error");
    
    for (index_type i=0; i<sphere.n(); ++i) {
        const Vec<3> pc = sphere.physCrd(i);
        sphere.setScalarFieldValue("relvort", i, vorfn.evaluateScalar(pc));
        sphere.setVectorFieldValue("velocity", i, velfn.evaluateVector(pc));
        sphere.setScalarFieldValue("stream_fn", i, vorfn.evaluateScalar(pc));
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
    vec_view_type crds;
    vec_host_view_type hcrds;
    scalar_view_type wgts;
    scalar_host_view_type hwgts;
    sphere.init_pack_all_coords(crds, hcrds, wgts, hwgts);
    
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
        if (token == "--init_nest" || token == "-in") {
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