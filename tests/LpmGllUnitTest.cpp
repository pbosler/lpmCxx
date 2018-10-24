#include "LpmGll.hpp"
#include "LpmAosTypes.hpp"
#include "LpmUtilities.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <exception>

using namespace Lpm;

int main(int argc, char* argv[]) {
    {// CubicGLL<2> test
    CubicGLL<2> gll;

    std::cout << "One-d cubic reference element:" << std::endl;
    std::cout << "quadrature points:" << std::endl;
    std::cout << "\t";
    for (int j=0; j<4; ++j) {
        std::cout << gll.qp4(j) << " ";
    }
    std::cout << std::endl;
    std::cout << "quadrature weights:" << std::endl;
    std::cout << "\t";
    scalar_type sum1d = 0.0;
    for (int j=0; j<4; ++j) {
        std::cout << gll.qw4(j) << " ";
        sum1d += gll.qw4(j);
    }
    std::cout << std::endl;
    std::cout << "sum1d = " << sum1d << std::endl;
    
    if (std::abs(sum1d - 2.0) > 1.0e-15) {
        throw std::runtime_error("sum1d is incorrect.");
    }
    
    Aos::Vec<2> refqp;
    scalar_type sum2d = 0.0;
    std::cout << "Two-d quad16 reference element edge points, weights:" << std::endl;
    for (int j=0; j<12; ++j) {
        refqp = gll.quad16edgeqp(j);
        sum2d += gll.quad16edgeqw(j);
        std::cout << "\t" << refqp << " , " << gll.quad16edgeqw(j) << std::endl;
    }
    
    std::cout << "Two-d quad16 reference element interior points, weights:" << std::endl;
    for (int j=0; j<4; ++j) {
        refqp = gll.quad16centerqp(j);
        sum2d += gll.quad16centerqw(j);
        std::cout << "\t" << refqp << " , " << gll.quad16centerqw(j) << std::endl;
    }
    std::cout << "sum2d = " << sum2d << std::endl;
    
    if (std::abs(sum2d - 4.0) > 1.0e-15) {
        throw std::runtime_error("sum2d is incorrect.");
    }
    
    
    scalar_type r1, r2;
    const scalar_type a = 0.5;
    const scalar_type b = 3.0;
    const scalar_type c = 2.0;
    quadraticRoots(r1, r2, a, b, c);
    std::cout << "x^2/2 + 3x + 2 = 0 at x = -3 +/- sqrt{5}: x1 = -5.23607, x2 = -0.763932" << std::endl;
    std::cout << "\tcomputed roots: x1 = " << r1 << ", x2 = " << r2 << std::endl;
    std::cout << "\tpicked root = " << gll.pickRoot(r1, r2) << std::endl;
    if (std::abs(r1 + 0.76393202250021030359) > 1.0e-15) {
        throw std::runtime_error("quadratic root error.");
    }
    
    std::vector<Aos::Vec<2>> corners(4);
    corners[0] = Aos::Vec<2>(2,3);
    corners[1] = Aos::Vec<2>(1.5, 1);
    corners[2] = Aos::Vec<2>(2.25, 1.125);
    corners[3] = Aos::Vec<2>(2.5, 2.75);
    const scalar_type quad_area = 1.15625;
    std::cout << "Corners = " << std::endl;
    for (int i=0; i<4; ++i) {
        std::cout << "Vec: " << corners[i] << std::endl;
        std::cout << "comps: (" << corners[i][0] << ", " << corners[i][1] << ")" << std::endl;
    }
    
    Aos::Vec<2> mappt;
    Aos::Vec<2> invpt;
    scalar_type quad_ar = 0.0;
    for (int i=0; i<12; ++i) {
        const Aos::Vec<2> refcrd = gll.quad16edgeqp(i);
        mappt = gll.bilinearMap(corners, refcrd[0], refcrd[1]);
        invpt = gll.invertBilinearMap(corners, mappt);
        const scalar_type jj = gll.bilinearMapJacobian(refcrd, corners);
        quad_ar += jj*gll.quad16edgeqw(i);
        
        std::cout << "refcrd   " << refcrd << " maps to " << mappt << std::endl;
        std::cout << "invmap = " << invpt << ", |err| = " << (refcrd - invpt).mag() << std::endl;
        std::cout << "jac = " << jj << std::endl;
        if ( (refcrd - invpt).mag() > 1.0e-15) {
            std::ostringstream ss;
            ss << "edge pt " << i << ": inverse map error magnitude = " << (refcrd - invpt).mag();
            throw std::runtime_error(ss.str());
        }
    }
    
    for (int i=0; i<4; ++i) {
        const Aos::Vec<2> refcrd = gll.quad16centerqp(i);
        mappt = gll.bilinearMap(corners, refcrd[0], refcrd[1]);
        invpt = gll.invertBilinearMap(corners, mappt);
        const scalar_type jj = gll.bilinearMapJacobian(refcrd, corners);
        quad_ar += jj*gll.quad16centerqw(i);
        
        std::cout << "refcrd   " << refcrd << " maps to " << mappt << std::endl;
        std::cout << "invmap = " << invpt << ", |err| = " << (refcrd - invpt).mag() <<  std::endl;
        std::cout << "jac = " << jj << std::endl;
        if ( (refcrd - invpt).mag() > 1.0e-15) {
            std::ostringstream ss;
            ss << "interior pt " << i << ": inverse map error magnitude = " << (refcrd - invpt).mag();
            throw std::runtime_error(ss.str());
        }
    }
    std::cout << "area = " << quad_ar << ", |err| = " << std::abs(quad_ar-quad_area) << std::endl;
    if (std::abs(quad_area - quad_ar) > 1.0e-15) {
        throw std::runtime_error("Area integral error.");
    }
    }// <2>
    {// CubicGLL<3> test
    }
return 0;
}

