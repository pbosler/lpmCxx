#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmOutputMessage.h"
#include "LpmLogger.h"
#include "LpmSakajo2009.h"
#include "LpmTimer.h"
#include "LpmVtkFileIO.h"
#include "LpmMultiIndex.h"
#include "LpmOctree.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <memory>
#include <fstream>
#include <mpi.h>
#include <cmath>

using namespace Lpm;

inline scalar_type legendre54(scalar_type z) {
    return z*square(square(z)-1.0);
}

struct Input {
    Input(int argc, char* argv[]);

    scalar_type delta;
    int tree_depth_l;
    int series_order;
    std::string program_name;
    index_type nLon;
    index_type nLat;
    scalar_type nu_exp;
    
    std::string infoString() const;
};

int main (int argc, char* argv[]) {
    int mpiErrCode;
    int numProcs;
    int procRank;
    mpiErrCode = MPI_Init(&argc, &argv);
    mpiErrCode = MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    mpiErrCode = MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    
    Input input(argc, argv);
    
    Timer programTimer("Sakajo 2009 Unit Tests");
    programTimer.start();

    std::unique_ptr<Logger> log(new Logger(OutputMessage::debugPriority, "SakajoUnitTests_log", procRank));
    std::stringstream ss;
    const std::string nullstr;
    {
        ss << "Test info: \n \t title: " << "Unit tests for Sakajo 2009 data structures" << std::endl;
//         ss << "\t objectives: " << std::endl;
//         ss << "\t 1. *objective 1*" << std::endl;
        OutputMessage introMsg(ss.str(), OutputMessage::tracePriority, "main");
        log->logMessage(introMsg);
        ss.str(nullstr);
    }
    {
        const int ts_deg = 3;
        Box3d testBox(0.37, 0.38, 0.525, 0.55, 0.74, 0.76);
        SakajoNode node(testBox, NULL, ts_deg+1);
        
        const scalar_type h = 0.25;
        const scalar_type nu = 1.0;
        const XyzVector tgt(-std::sqrt(3.0)/3, -std::sqrt(3.0)/3, -std::sqrt(3.0)/3);
        std::vector<XyzVector> ccs = testBox.corners();
        std::cout << "|corners| = ";
        for (int i=0; i<8; ++i) {
            std::cout << ccs[i].magnitude() << "  ";
        }
        std::cout << std::endl;
        const scalar_type dist0 = distance(testBox.centroid, tgt);
        std::cout << "box intersects sphere: " << (testBox.intersectsSphere() ? "yes" : "no") << std::endl;
        std::cout << "centroid = " << testBox.centroid << ", norm = " << testBox.centroid.magnitude() << std::endl;
        std::cout << "dist(testBox.centroid, tgt) = " << dist0 << std::endl;
        std::cout << "h**nu * dist = " << std::pow(h,nu) * dist0 << std::endl;
        std::cout << "box radius = " << node.box.maxRadius << std::endl;
        std::cout << "MAC = " << (node.multipoleAcceptance(tgt, h, nu) ? "true" : "false") << std::endl;
        std::cout << "D(x,y_T) = " << 1.0/(1 - tgt.dotProduct(testBox.centroid)) << std::endl;
        
        const scalar_type delta = 0.0;
        const scalar_type sph_radius = 1.0;
        node.calc_coeffs(tgt, sph_radius, delta);
        
        XyzVector srcVec = XyzVector(0.376867, 0.538382, 0.753735);
        srcVec.normalize();
        std::cout << srcVec - testBox.centroid << "," << srcVec.magnitude() << std::endl;

        std::cout << "Box contains source " << srcVec << " ? " << (testBox.containsPoint(srcVec) ? "yes" : "no") << std::endl;
        
        const scalar_type gamma = 1.0;
        node.calc_moments(srcVec, gamma);
        
        std::cout << node.infoString();
        
        XyzVector nodesum(0.0, 0.0, 0.0);
        for (auto& elem : node.coeffs) {
            const std::vector<scalar_type> mm = node.moments.at(elem.first);
            nodesum.x += elem.second * (tgt.y*mm[2] - tgt.z*mm[1]);
            nodesum.y += elem.second * (tgt.z*mm[0] - tgt.x*mm[2]);
            nodesum.z += elem.second * (tgt.x*mm[1] - tgt.y*mm[0]);
        }
        
        const XyzVector bs = biotSavart(tgt, srcVec, sph_radius, delta);
        const XyzVector err = nodesum - bs;

        std::cout << (delta>0.0 ? "Smoothed " : "Singular ") << "point source at " << srcVec << ", induced velocity at " << tgt << std::endl;
        std::cout << "Appx. velocity = " << nodesum.scalarMultiply(-1.0/(4.0*PI*sph_radius))  
                  << " kernel velocity = " << bs.scalarMultiply(-1.0/(4.0*PI*sph_radius)) 
                  << " |rel. err.| = " << err.magnitude()/bs.magnitude() << std::endl;
    }
    const int nLat = input.nLat;
    const int nLon = input.nLon;
    const int seriesOrder = input.series_order;
    const int treeDepthParam = input.tree_depth_l;
    
    const scalar_type nu = input.nu_exp;
    const scalar_type delta = input.delta;
    
    const scalar_type meshSize = 1.0 / std::pow(2.0, input.tree_depth_l);
    
    std::cout << input.infoString();
    
    const scalar_type dz = 1.0 / (0.5 * nLat);
    
    std::shared_ptr<SphericalCoords> sc(new SphericalCoords(nLat*nLon));
    for (int i=1; i<nLat; ++i) {
        const scalar_type zz = 1.0 - (i) * dz;
        for (int j=0; j < nLon; ++j) {
            const scalar_type x = std::sqrt(1.0-square(zz)) * std::cos(2.0 * PI * j / nLon);
            const scalar_type y = std::sqrt(1.0-square(zz)) * std::sin(2.0 * PI * j / nLon);
            sc->insert(x,y,zz);
//             std::cout << "getVec(" << i+j << ") : (" << x << ", " << y << ", " << zz << ")  =  " << sc->getVec(i+j) << std::endl;
        }
    }
    std::cout << "sc->n(): " << sc->n() << std::endl;
    std::shared_ptr<Field> circ(new Field(sc->n(), 1, "Gamma", "1/time"));
    for (index_type i=0; i<sc->n(); ++i) {
        //circ->insert( (i%2 == 0 ? 1.0 : -1.0));
        const XyzVector cvec = sc->getVec(i);
        circ->insert(30.0 * std::cos(4.0 * longitude(cvec)) * legendre54(cvec.z) * 4.0 * PI / (sc->n()));
    }
    
    std::shared_ptr<Field> exactVelocity(new Field(sc->n(), 3, "u", "length/time"));
    for (index_type i=0; i<sc->n(); ++i) {
        const XyzVector cvec = sc->getVec(i);
        const scalar_type lat = latitude(cvec);
        const scalar_type lon = longitude(cvec);
        const scalar_type u = 0.5 * std::cos(4.0*lon) * cube(std::cos(lat))*(5.0 * std::cos(2.0*lat) - 3.0);
        const scalar_type v = 4.0 * cube(std::cos(lat))*std::sin(lat)*std::sin(4.0*lon);
        exactVelocity->insert(-u*std::sin(lon) - v*std::sin(lat)*std::cos(lon), u*std::cos(lon) - v*std::sin(lat)*std::sin(lon),
            v*std::cos(lat));
    }
    
    Timer buildTimer("Tree build timer");
    buildTimer.start();
    SakajoTree tree(seriesOrder, treeDepthParam, 1.0, delta, procRank);
    buildTimer.end();
    std::cout << buildTimer.infoString();
    
    std::cout << "treecode params: h = " << meshSize << ", nu = " << nu << ", delta = " << delta << std::endl;
    std::cout << tree.Tree::infoString();
    
    std::shared_ptr<Field> velDirect(new Field(sc->n(), 3, "direct_sum_velocity", "dist/time"));
    velDirect->initializeToConstant(sc.get());
    Timer directSumTimer("Direct sum timer");
    directSumTimer.start();
    for (index_type i=0; i<sc->n(); ++i) {
        const XyzVector tgtVec = sc->getVec(i);
        XyzVector vel(0.0, 0.0, 0.0);
        for (index_type j=0; j<sc->n(); ++j) {
            //if (i != j ) {
            const XyzVector srcVec = sc->getVec(j);
            const scalar_type Gamma = circ->getScalar(j);   
            const XyzVector kernel = biotSavart(tgtVec, srcVec, 1.0, delta);
            vel += kernel.scalarMultiply(-Gamma / (4.0 * PI ));
            //}
        }
//         for (index_type j=i+1; j<sc->n(); ++j) {
//             const XyzVector srcVec = sc->getVec(j);
//             const scalar_type Gamma = circ->getScalar(j);
//             const XyzVector kernel = tree.biotSavart(tgtVec, srcVec, delta);
//             vel += kernel.scalarMultiply(Gamma);
//         }
        velDirect->replace(i, vel);
    }
    directSumTimer.end();
    std::cout << directSumTimer.infoString();
    
    std::shared_ptr<Field> velTree(new Field(sc->n(), 3, "tree_sum_velocity", "dist/time"));
    velTree->initializeToConstant(sc.get());
    
    
    tree.writeToVtk("sakajoTree.vtk");
    Timer treecodeTimer("Tree sum timer");
    treecodeTimer.start();
    tree.computeMoments(sc, circ);
    
    tree.computeVelocity(velTree, sc, circ, meshSize, nu);
    //tree.printAll();
    
    treecodeTimer.end();
    std::cout << treecodeTimer.infoString();
    
    std::shared_ptr<Field> directSumError(new Field(sc->n(), 3, "direct_sum_error", "dist/time"));
    directSumError->initializeToConstant(sc.get());
    directSumError->update(1.0, exactVelocity, -1.0, velDirect);
    std::cout << "direct sum max err = " << directSumError->maxMagnitude() << std::endl;
    
    std::shared_ptr<Field> treecodeError(new Field(sc->n(), 3, "treecode_error", "dist/time"));
    treecodeError->initializeToConstant(sc.get());
    treecodeError->update(1.0, velDirect, -1.0, velTree);
    std::cout << "max(error) = " << treecodeError->maxMagnitude() << std::endl;
   
    std::vector<std::shared_ptr<Field>> fields = {circ, velDirect, velTree, treecodeError, exactVelocity, directSumError};
    const std::string outfname = "sakajoUnitTest.csv";
    std::ofstream fs(outfname);
    VtkWriter writer;
    writer.writePointsAndFieldsToCSV(fs, sc, fields);  
    fs.close();  
    
    programTimer.end();
    OutputMessage finalMsg("PROGAM COMPLETE: " + programTimer.infoString(), OutputMessage::remarkPriority, "main");
    log->logMessage(finalMsg);
    MPI_Finalize();
return 0;
}

Input::Input(int argc, char* argv[]) {
    program_name = argv[0];
    delta = 0.0;
    series_order = 3;
    nLat = 64;
    nLon = 64;
    nu_exp = 1.0;
    tree_depth_l = 4;
    for (int i=1; i<argc; ++i) {
        const std::string& token = argv[i];
        if (token == "--delta") {
            delta = std::stod(argv[++i]);
        }
        else if (token == "--order" || token == "--lambda" || token == "--lam" ) {
            series_order = std::stoi(argv[++i]);
        }
        else if (token == "--nLat") {
            nLat = std::stoi(argv[++i]);
        }
        else if (token == "--nLon") {
            nLon = std::stoi(argv[++i]);
        }
        else if (token == "--tree" || token == "--tree_depth_param" || token == "--l") {
            tree_depth_l = std::stoi(argv[++i]);
        }
        else if (token == "--nu") {
            nu_exp = std::stod(argv[++i]);
        }
    }
}

std::string Input::infoString() const {
    std::stringstream ss;
    ss << program_name << " --delta " << delta << " --order " << series_order << " --nLat " << nLat << " --nLon " << nLon
       << " --tree " << tree_depth_l << " --nu " << nu_exp << std::endl;
    ss << "\t\tS2009: log_2(N) = tree_depth_l = " << std::log2(nLat*nLon) << " nu = " << 1.0/std::log2(nLat*nLon) 
       << " particle_dist = " << std::sqrt(4.0*PI/nLat/nLon) << std::endl;
    
    return ss.str();
}