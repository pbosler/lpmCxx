#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <memory>
#include "LpmAosTypes.hpp"
#include "LpmAosParticle.hpp"
#include "LpmAosParticleSet.hpp"
#include "LpmAosParticleFactory.hpp"
#include "LpmTypeDefs.h"

using namespace Lpm::Aos;

using Lpm::index_type;
using Lpm::scalar_type;
using Lpm::PI;

struct Input {
    index_type nmax;
    std::string prog_name;
    std::string ofroot; 
    index_type wi;
    
    Input(int argc, char* argv[]);
    std::string statusMsg() const;
};

struct Output {    
    Output() {}
    void writePyData(std::ostream& os, const ParticleSet<1>& active, const ParticleSet<1>& passive) const;
};

inline scalar_type bottomTopo(const scalar_type x) {return 0.1*x*x;}

int main(int argc, char* argv[]) {

    Input input(argc, argv);
    std::cout << input.statusMsg();
    Output output;
    
    std::shared_ptr<ParticleFactory<1>> fac(new SWEParticleFactory<1>());
    
    ParticleSet<1> active(fac, input.nmax);
    ParticleSet<1> passive(fac, input.nmax+1);


    {
        const scalar_type len = 2.0*PI/input.nmax;
        const scalar_type dx = len;
    
        std::vector<scalar_type> asb(input.nmax,0.0);
        std::vector<scalar_type> psb(input.nmax+1, 0.0);
        std::vector<scalar_type> ah(input.nmax, 0.0);
        std::vector<scalar_type> ph(input.nmax+1, 0.0);
        std::vector<scalar_type> as(input.nmax, 0.0);
        std::vector<scalar_type> ps(input.nmax+1, 0.0);
        Vec<1> x0;

        for (index_type i=0; i<input.nmax; ++i) {
            x0.x[0] = -PI + (i+0.5)*dx;
            asb[i] = bottomTopo(x0.x[0]);
            as[i] = (std::abs(x0.x[0]) <= 2.0 ? 0.4 : asb[i]);
            ah[i] = as[i] - asb[i];
            //std::cout << "i = " << i << ": x0 = " << x0 << std::endl;
            active.insert(x0, len);
        }
        
        for (index_type i=0; i<=input.nmax; ++i) {
            x0.x[0] = -PI + i*dx;
            psb[i] = bottomTopo(x0.x[0]);
            ps[i] = (std::abs(x0.x[0]) <= 2.0 ? 0.4 : psb[i]);
            ph[i] = ps[i] - psb[i];
            //std::cout << "i = " << i << ": x0 = " << x0 << std::endl;
            passive.insert(x0);
        }
        active.initScalarFieldFromVector("bottom_height", asb);
        active.initScalarFieldFromVector("depth", ah);
        active.initScalarFieldFromVector("surf_height", as);
        passive.initScalarFieldFromVector("bottom_height", psb);
        passive.initScalarFieldFromVector("depth", ph);
        passive.initScalarFieldFromVector("surf_height",ps);
    }
    
    std::cout << active.infoString() << std::endl;
    std::cout << passive.infoString() << std::endl;
    
    std::ofstream fs("tmp/swe1d_0.py");
    output.writePyData(fs, active, passive);
    fs.close();
return 0;
}

Input::Input(int argc, char* argv[]) {
    prog_name = argv[0];
    nmax = 10;
    ofroot = "null";
    wi = 1;
    for (int i=1; i<argc; ++i) {
        const std::string& token = argv[i];
        if (token == "-nmax" || token == "-n") {
            nmax = std::stoi(argv[++i]);
        }
        else if (token == "-o") {
            ofroot = std::string(argv[++i]);
        }
        else if (token == "-wi") {
            wi = std::stoi(argv[++i]);
        }
    }
}

std::string Input::statusMsg() const {
    std::ostringstream ss;
    ss << prog_name << " input summary: " << std::endl;
    ss << "\tnmax = " << nmax << std::endl;
    ss << "\toutput file root = " << ofroot << std::endl;
    ss << "\twrite interval = " << wi << std::endl;
    return ss.str();
}

void Output::writePyData(std::ostream& os, const ParticleSet<1>& active, const ParticleSet<1>& passive) const {
    os << "import numpy as np" << std::endl;
    {
        os << "activex = np.array([";
        for (index_type i=0; i<active.n(); ++i) {
            const Vec<1> x = active.physCrd(i);
            os << x.x[0] << (i<active.n()-1? "," : "");
        }
        os << "])" << std::endl;
        os << "active_lag = np.array([";
        for (index_type i=0; i<active.n(); ++i) {
            const Vec<1> a = active.lagCrd(i);
            os << a.x[0] << (i<active.n()-1 ? "," : "");
        }
        os << "])" << std::endl;
        os << "active_sb = np.array([";
        const std::vector<scalar_type> sb = active.getScalarFieldValues("bottom_height");
        for (index_type i=0; i<active.n(); ++i) {
            os << sb[i] << (i<active.n()-1 ? "," : "");
        }
        os << "])" << std::endl;
        const std::vector<scalar_type> h = active.getScalarFieldValues("depth");
        const std::vector<scalar_type> s = active.getScalarFieldValues("surf_height");
        os << "active_h = np.array([";
        for (index_type i=0; i<active.n(); ++i) {
            os << h[i] << (i<active.n()-1 ? "," : "");
        }
        os << "])" << std::endl;
        os << "active_s = np.array([";
        for (index_type i=0; i<active.n(); ++i) {
            os << s[i] << (i<active.n()-1 ? "," : "");
        }
        os << "])" << std::endl;
        const std::vector<scalar_type> div = active.getScalarFieldValues("divergence");
        os << "active_div = np.array([";
        for (index_type i=0; i<active.n(); ++i) {
            os << div[i] << (i<active.n()-1 ? "," : "");
        }
        os << "])" << std::endl;
        os << "active_len = np.array([";
        for (index_type i=0; i<active.n(); ++i) {
            os << active.getPtr(i)->weight() << (i<active.n()-1 ? "," : "");
        }
        os << "])" << std::endl;
        const std::vector<scalar_type> u = active.getScalarFieldValues("velocity");
        os << "active_u = np.array([";
        for (index_type i=0; i<active.n(); ++i) {
            os << u[i] << (i<active.n()-1 ? "," : "");
        }
        os << "])" << std::endl;
    }
    {
        os << "passivex = np.array([";
        for (index_type i=0; i<passive.n(); ++i) {
            const Vec<1> x = passive.physCrd(i);
            os << x.x[0] << (i<passive.n()-1 ? "," : "");
        }
        os << "])" << std::endl;
        os << "passive_lag = np.array([";
        for (index_type i=0; i<passive.n(); ++i) {
            const Vec<1> a = passive.lagCrd(i);
            os << a.x[0] << (i<passive.n()-1 ? "," : "");
        }
        os << "])" << std::endl;
        os << "passive_sb = np.array([";
        const std::vector<scalar_type> sb = passive.getScalarFieldValues("bottom_height");
        for (index_type i=0; i<passive.n(); ++i) {
            os << sb[i] << (i<passive.n()-1 ? "," : "");
        }
        os << "])" << std::endl;
        const std::vector<scalar_type> h = passive.getScalarFieldValues("depth");
        const std::vector<scalar_type> s = passive.getScalarFieldValues("surf_height");
        os << "passive_h = np.array([";
        for (index_type i=0; i<passive.n(); ++i) {
            os << h[i] << (i<passive.n()-1 ? "," : "");
        }
        os << "])" << std::endl;
        os << "passive_s = np.array([";
        for (index_type i=0; i<passive.n(); ++i) {
            os << s[i] << (i<passive.n()-1 ? "," : "");
        }
        os << "])" << std::endl;
        const std::vector<scalar_type> div = passive.getScalarFieldValues("divergence");
        os << "passive_div = np.array([";
        for (index_type i=0; i<passive.n(); ++i) {
            os << div[i] << (i<passive.n()-1 ? "," : "");
        }
        os << "])" << std::endl;
        const std::vector<scalar_type> u = passive.getScalarFieldValues("velocity");
        os << "passive_u = np.array([";
        for (index_type i=0; i<passive.n(); ++i) {
            os << u[i] << (i<passive.n()-1 ? "," : "");
        }
        os << "])" << std::endl;
    }
    {        
        os << "import matplotlib.pyplot as plt" << std::endl;
        os << "fig, ax = plt.subplots(2,3, sharex=True, sharey=False)" << std::endl;
        os << "" << std::endl;
        os << "ax[0,0].plot(activex, active_s, '-o',label='$s_a$')" << std::endl;
        os << "ax[0,0].plot(passivex, passive_s, '-+', label='$s_p$')" << std::endl;
        os << "ax[0,0].plot(passivex, passive_sb, '-', label='$s_b$')" << std::endl;
        os << "" << std::endl;
        os << "ax[0,0].legend()" << std::endl;
        os << "" << std::endl;
        os << "ax[0,1].plot(activex, active_h, '-o', label='$h_a$')" << std::endl;
        os << "ax[0,1].plot(passivex, passive_h, '-+', label='$h_p$')" << std::endl;
        os << "ax[0,1].legend()" << std::endl;
        os << "" << std::endl;
        os << "ax[1,0].plot(activex, active_div, '-o', label=r'$\\delta_a$')" << std::endl;
        os << "ax[1,0].plot(passivex, passive_div, '-+', label=r'$\\delta_p$')" << std::endl;
        os << "ax[1,0].legend()" << std::endl;
        os << "" << std::endl;
        os << "ax[1,1].plot(activex, active_len, '-o', label='$l_a$')" << std::endl;
        os << "ax[1,1].legend()" << std::endl;
        os << "" << std::endl;
        os << "ax[0,2].plot(activex, active_u, '-o', label='$u_a$')" << std::endl;
        os << "ax[0,2].plot(passivex, passive_u,'-+', label='$u_p$')" << std::endl;
        os << "ax[0,2].legend()" <<  std::endl;
        os << "plt.show()" << std::endl;
        os << "plt.close(fig)" << std::endl;
    }
}

