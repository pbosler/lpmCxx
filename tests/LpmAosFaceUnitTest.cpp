#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmOutputMessage.h"
#include "LpmLogger.h"
#include "LpmAosFace.hpp"
#include "LpmAosFaceFactory.hpp"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <memory>
#include <vector>

using namespace Lpm;

typedef std::vector<index_type> ind_vec;

int main (int argc, char* argv[]) {

    std::unique_ptr<Logger> log(new Logger(OutputMessage::debugPriority));
    std::stringstream ss;
    const std::string nullstr;
    {
        ss << "Test info: \n \t title: " << "Lpm Aos Face Unit Test" << std::endl;
        ss << "\t objectives: " << std::endl;
        ss << "\t 1. Verify basic functionality of Lpm::Face class." << std::endl;
        OutputMessage introMsg(ss.str(), OutputMessage::tracePriority, "main");
        log->logMessage(introMsg);
        ss.str(nullstr);
    }
    
    QuadFaceFactory<3> quadFactory;
    TriFaceFactory<2> triFactory;
    
    const ind_vec ctr = {0};
    const ind_vec triverts = {0, 1, 2};
    const ind_vec quadverts = {0, 1, 2, 3};
    const ind_vec triedges = {10, 11, 12};
    const ind_vec quadedges = {10, 11, 12, 13};
    const index_type prt = 0;
    const scalar_type ar = 1.0;
    
    
    std::cout << "test start" << std::endl;
    std::cout << "Quad Face: " << std::endl;
    std::unique_ptr<Face<3>> qf = quadFactory.createFace(ctr, quadverts, quadedges, prt, ar);
    std::cout << qf->infoString();

    std::cout << "Tri face: " << std::endl;
    std::unique_ptr<Face<2>> tf = triFactory.createFace(ctr, triverts, triedges, prt, ar);
    std::cout << tf->infoString();

return 0;
}

