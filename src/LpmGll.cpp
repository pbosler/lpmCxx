#include "LpmGll.hpp"
#include "LpmUtilities.h"
#include <limits>

namespace Lpm {

template <int ndim> scalar_type CubicGLL<ndim>::qp4(const int id) const {
    switch (id) {
        case(0): {return qp0; break;} 
        case(1): {return qp1; break;} 
        case(2): {return qp2; break;} 
        case(3): {return qp3; break;}
    }    
}

template <int ndim> scalar_type CubicGLL<ndim>::qw4(const int id) const {
    switch (id) {
        case(0): {return qw0; break;} 
        case(1): {return qw1; break;} 
        case(2): {return qw2; break;} 
        case(3): {return qw3; break;}
    }    
}

/** Inverts bilinearMap to find reference coordinates (r1,r2) in [-1,1]^2
 Folows C. Hua, 1990, Finite Elem. Anal. Design 7:159--166.
*/
template <int ndim> Aos::Vec<ndim> CubicGLL<ndim>::quad16edgeqp(const int id) const {
    switch (id) {
        case(0) : {
            return Aos::Vec<ndim>(qp4(0), qp4(3));
            break;
        }
        case(1) : {
            return Aos::Vec<ndim>(qp4(0), qp4(2));
            break;
        }
        case(2) : {
            return Aos::Vec<ndim>(qp4(0), qp4(1));
            break;
        }
        case(3) : {
            return Aos::Vec<ndim>(qp4(0), qp4(0));
            break;
        }
        case(4) : {
            return Aos::Vec<ndim>(qp4(1), qp4(0));
            break;
        }
        case(5) : {
            return Aos::Vec<ndim>(qp4(2), qp4(0));
            break;
        }
        case(6) : {
            return Aos::Vec<ndim>(qp4(3), qp4(0));
            break;
        }
        case(7) : {
            return Aos::Vec<ndim>(qp4(3), qp4(1));
            break;
        }
        case(8) : {
            return Aos::Vec<ndim>(qp4(3), qp4(2));
            break;
        }
        case(9) : {
            return Aos::Vec<ndim>(qp4(3), qp4(3));
            break;
        }
        case(10) : {
            return Aos::Vec<ndim>(qp4(2), qp4(3));
            break;
        }
        case(11) : {
            return Aos::Vec<ndim>(qp4(1), qp4(3));
            break;
        }
    }
}

template <int ndim> Aos::Vec<ndim> CubicGLL<ndim>::quad16centerqp(const int id) const {
    switch (id) {
        case(0) : {
            return Aos::Vec<ndim>(qp4(1), qp4(2));
            break;
        }
        case(1) : {
            return Aos::Vec<ndim>(qp4(1), qp4(1));
            break;
        }
        case(2) : {
            return Aos::Vec<ndim>(qp4(2), qp4(1));
            break;
        }
        case(3) : {
            return Aos::Vec<ndim>(qp4(2), qp4(2));
            break;
        }
    }
}

template <int ndim> scalar_type CubicGLL<ndim>::quad16edgeqw(const int id) const {
    switch (id) {
        case(0): {
            return qw4(0)*qw4(3);
            break;
        }
        case(1): {
            return qw4(0)*qw4(2);
            break;
        }
        case(2): {
            return qw4(0)*qw4(1);
            break;
        }
        case(3): {
            return qw4(0)*qw4(0);
            break;
        }
        case(4): {
            return qw4(1)*qw4(0);
            break;
        }
        case(5): {
            return qw4(2)*qw4(0);
            break;
        }
        case(6): {
            return qw4(3)*qw4(0);
            break;
        }
        case(7): {
            return qw4(3)*qw4(1);
            break;
        }
        case(8): {
            return qw4(3)*qw4(2);
            break;
        }
        case(9): {
            return qw4(3)*qw4(3);
            break;
        }
        case(10): {
            return qw4(2)*qw4(3);
            break;
        }
        case(11): {
            return qw4(1)*qw4(3);
            break;
        }
    }
}

template <int ndim> scalar_type CubicGLL<ndim>::quad16centerqw(const int id) const {
    switch (id) {
        case (0) : {
            return qw4(1)*qw4(2);
            break;
        }
        case (1) : {
            return qw4(1)*qw4(1);
            break;
        }
        case (2) : {
            return qw4(2)*qw4(1);
            break;
        }
        case (3) : {
            return qw4(2)*qw4(2);
            break;
        }
    }
}


template <int ndim> Aos::Vec<ndim> CubicGLL<ndim>::bilinearMap(const std::vector<Aos::Vec<ndim>>& corners, 
    const scalar_type s1, const scalar_type s2) const {
    return (corners[0].scale((1-s1)*(1+s2)) + corners[1].scale((1-s1)*(1-s2)) 
        + corners[2].scale((1+s1)*(1-s2)) + corners[3].scale((1+s1)*(1+s2))).scale(0.25);
}

template <int ndim> Aos::Vec<ndim> CubicGLL<ndim>::sphereBilinearMap(const std::vector<Aos::Vec<ndim>>& corners, 
    const scalar_type s1, const scalar_type s2, const scalar_type radius) const {
    return bilinearMap(corners, s1, s2).normalize().scale(radius);
}

template <int ndim> void CubicGLL<ndim>::pickBilinearIJ(int& i, int& j, const std::vector<Aos::Vec<ndim>>& corners) const {
    i = 0; 
    j = 1;
    if (ndim >2) {
        scalar_type minx = std::numeric_limits<scalar_type>::max();
        scalar_type maxx = std::numeric_limits<scalar_type>::lowest();
        scalar_type miny = minx;
        scalar_type maxy = maxx;
        scalar_type minz = minx;
        scalar_type maxz = maxx;
        for (int i=0; i<4; ++i) {
            maxx = std::max(corners[i][0], maxx);
            minx = std::min(corners[i][0], minx);
            maxy = std::max(corners[i][1], maxy);
            miny = std::min(corners[i][1], miny);
            maxz = std::max(corners[i][2], maxz);
            minz = std::min(corners[i][2], minz);
        }
        const scalar_type dx = maxx - minx;
        const scalar_type dy = maxy - miny;
        const scalar_type dz = maxz - minz;
    
        if (dy < dx && dy < dz) 
            j=2;
        else if (dx < dy && dx < dz) {
            i = 1;
            j = 2;
        }
    }
}

template <int ndim> scalar_type CubicGLL<ndim>::pickRoot(const scalar_type r1, const scalar_type r2) const {
    if (-1.0 <= r1 && r1 <= 1.0) 
        return r1;
    else
        return r2;
}

template <int ndim> Aos::Vec<ndim> CubicGLL<ndim>::invertBilinearMap(const std::vector<Aos::Vec<ndim>>& corners, const Aos::Vec<ndim>& queryPt) const {
    Aos::Vec<ndim> aa, bb, cc, dd;
    int i, j;
    scalar_type xi, eta, a_b, a_c, a_d, b_c, b_d, c_d;
    scalar_type r1, r2;
    
    aa = corners[0].scale(-1) + corners[1] - corners[2] + corners[3];
    bb = corners[0].scale(-1) - corners[1] + corners[2] + corners[3];
    cc =  corners[0] - corners[1] - corners[2] + corners[3];
    dd = queryPt.scale(4.0);
    for (int ii=0; ii<4; ++ii) {
        dd -= corners[ii];
    }
    
    pickBilinearIJ(i, j, corners);
    xi = std::numeric_limits<scalar_type>::max();
    eta = std::numeric_limits<scalar_type>::max();
    
    a_b = twoByTwoDeterminant(aa[i], aa[j], bb[i], bb[j]);
    a_c = twoByTwoDeterminant(aa[i], aa[j], cc[i], cc[j]);
    a_d = twoByTwoDeterminant(aa[i], aa[j], dd[i], dd[j]);
    b_c = twoByTwoDeterminant(bb[i], bb[j], cc[i], cc[j]);
    b_d = twoByTwoDeterminant(bb[i], bb[j], dd[i], dd[j]);
    c_d = twoByTwoDeterminant(cc[i], cc[j], dd[i], dd[j]);
    
    if (std::abs(aa[i]) < ZERO_TOL) {
        // case I
        if (std::abs(aa[j]) < ZERO_TOL) {
            // case I.A
            xi = -c_d / (aa[i]*dd[j] + b_c);
            eta = b_d / (aa[j]*dd[i] + b_c);
        }
        else if (std::abs(cc[i]) < ZERO_TOL) {
            // case I.B.a
            xi = dd[i] / bb[i];
            eta = (bb[i]*dd[j]-bb[j]*dd[i])/(aa[j]*dd[i]+bb[i]*cc[j]);
        }
        else {
            // case I.B.b
            quadraticRoots(r1, r2, aa[j]*bb[i], cc[j]*bb[i]-aa[j]*dd[i]-bb[j]*cc[i], dd[j]*cc[i]-cc[j]*dd[i]);
            xi = pickRoot(r1, r2);
            eta = (dd[i]-bb[i])*xi/cc[i];
        }
    }
    else {
        // case II
        if (std::abs(aa[j]) > ZERO_TOL) {
            // case II.A
            if (std::abs(a_b) > ZERO_TOL) {
                // case II.A.a
                if (std::abs(a_c) > ZERO_TOL) {
                    // case II.A.a.1
                    quadraticRoots(r1, r2, a_b, -b_c - a_d, -c_d);
                    xi = pickRoot(r1, r2);
                    eta = (a_d-a_b*xi)/(a_c);
                }
                else {
                    // case II.A.a.2
                    xi = a_d/a_b;
                    eta = -aa[i]*b_d/(cc[i]*a_b + aa[i]*a_d);
                }
            }
            else {
                // case II.A.b
                xi = -aa[i]*c_d/(bb[i]*a_c+aa[i]*a_d);
                eta = a_d/a_c;
            }
        }
        else {
            // case II.B
            if (std::abs(bb[j]) < ZERO_TOL) {
                // case II.B.a
                xi = -c_d/(aa[i]*dd[j]+bb[i]*cc[j]);
                eta = dd[j]/cc[j];
            }
            else {
                // case II.B.b
                quadraticRoots(r1, r2, aa[i]*bb[j], cc[i]*bb[j]-aa[i]*dd[j]-bb[i]*cc[j], dd[i]*cc[j]-cc[i]*dd[j]);
                xi = pickRoot(r1, r2);
                eta = (dd[j]-bb[j]*xi)/cc[j];
            }
        }
    }
    return Aos::Vec<ndim>(xi, eta);
}

template <int ndim> scalar_type CubicGLL<ndim>::quad16interpolation(const Aos::Vec<ndim>& refcrds, const std::vector<scalar_type>& edgevals, const std::vector<scalar_type>& ctrvals) const {
    return  edgevals[0]*basis0(refcrds[0])*basis3(refcrds[1]) +
            edgevals[1]*basis0(refcrds[0])*basis2(refcrds[1]) +
            edgevals[2]*basis0(refcrds[0])*basis1(refcrds[1]) +
            edgevals[3]*basis0(refcrds[0])*basis0(refcrds[1]) +
            edgevals[4]*basis1(refcrds[0])*basis0(refcrds[1]) +
            edgevals[5]*basis2(refcrds[0])*basis0(refcrds[1]) + 
            edgevals[6]*basis3(refcrds[0])*basis0(refcrds[1]) + 
            edgevals[7]*basis3(refcrds[0])*basis1(refcrds[1]) +
            edgevals[8]*basis3(refcrds[0])*basis2(refcrds[1]) +
            edgevals[9]*basis3(refcrds[0])*basis3(refcrds[1]) +
            edgevals[10]*basis2(refcrds[0])*basis3(refcrds[1]) +
            edgevals[11]*basis1(refcrds[0])*basis3(refcrds[1]) +
            ctrvals[0]*basis1(refcrds[0])*basis2(refcrds[1]) + 
            ctrvals[1]*basis1(refcrds[0])*basis1(refcrds[1]) +
            ctrvals[2]*basis2(refcrds[0])*basis1(refcrds[1]) +
            ctrvals[3]*basis2(refcrds[0])*basis2(refcrds[1]);
}

template <int ndim> scalar_type CubicGLL<ndim>::bilinearMapJacobian(const Aos::Vec<ndim>& refcrds, const std::vector<Aos::Vec<ndim>>& corners) const {
    const scalar_type a = 0.25*(-(1+refcrds[1])*corners[0][0] - 
                                 (1-refcrds[1])*corners[1][0] +
                                 (1-refcrds[1])*corners[2][0] + 
                                 (1+refcrds[1])*corners[3][0]);
        
    const scalar_type b = 0.25*(-(1+refcrds[1])*corners[0][1] - 
                                 (1-refcrds[1])*corners[1][1] +
                                 (1-refcrds[1])*corners[2][1] + 
                                 (1+refcrds[1])*corners[3][1]);

    const scalar_type c = 0.25*( (1-refcrds[0])*corners[0][0] - 
                                 (1-refcrds[0])*corners[1][0] -
                                 (1+refcrds[0])*corners[2][0] + 
                                 (1+refcrds[0])*corners[3][0]);
        
    const scalar_type d = 0.25*( (1-refcrds[0])*corners[0][1] - 
                                 (1-refcrds[0])*corners[1][1] -
                                 (1+refcrds[0])*corners[2][1] + 
                                 (1+refcrds[0])*corners[3][1]);
        
//      return std::abs(a*d - b*c);
//     std::cout << "det = " << a*d - b*c << std::endl;
   return twoByTwoDeterminant(a, b, c, d);
}

template struct CubicGLL<2>;
template struct CubicGLL<3>;

}
