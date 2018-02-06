#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmOutputMessage.h"
#include "LpmLogger.h"
#include "LpmXyzVector.h"
#include "LpmBox3d.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <memory>

using namespace Lpm;

int main (int argc, char* argv[]) {

    std::unique_ptr<Logger> log(new Logger(OutputMessage::debugPriority));
    std::stringstream ss;
    const std::string nullstr;
    {
        ss << "Test info: \n \t title: " << "Box3d unit test" << std::endl;
        ss << "\t objectives: " << std::endl;
        ss << "\t 1. Verify Box3d class functions" << std::endl;
        OutputMessage introMsg(ss.str(), OutputMessage::tracePriority, "main");
        log->logMessage(introMsg);
        ss.str(nullstr);
    }

    {   // box 3d unit tests
        
        Box3d unitBox(-1,1,-1,1,-1,1);
        std::cout << "Unit Box : " << unitBox.infoString();
        std::cout << "\tvolume = " << unitBox.volume() << std::endl;
        std::cout << "\tcentroid = " << unitBox.centroid() << std::endl;
        std::cout << "\tface centroids : " << std::endl;
        const std::vector<XyzVector> fctrs = unitBox.faceCentroids();
        std::cout << "\t\t";
        for (int i=0; i<6; ++i) 
            std::cout << fctrs[i] << "    ";
        std::cout << std::endl;
        std::cout << "\tcontains origin ? : " << (unitBox.containsPoint(XyzVector(0.0, 0.0, 0.0)) ? "true" : "false") << std::endl;
        std::cout << "\tlongestEdge = " << unitBox.longestEdge() << std::endl;
        std::cout << "\tshortestEdge = " << unitBox.shortestEdge() << std::endl;
        std::cout << "\taspectRatio = " << unitBox.aspectRatio() << std::endl;
        std::cout << "\tminRadius = " << unitBox.minRadius() << std::endl;
        std::cout << "\tmaxRadius = " << unitBox.maxRadius() << std::endl;
        
        std::vector<Box3d> kids = unitBox.bisectAll();
        for (int i = 0; i < 8; ++i) {
            std::cout << "child " << i << ": " << kids[i].infoString();
        }
        
        Box3d boxCopy(kids[1]);
        std::cout << "BoxCopy " << boxCopy.infoString();
        Box3d boxAssign = boxCopy;
        std::cout << "BoxAssign " << boxAssign.infoString();
        
        Box3d boxR(1.0);
        std::cout << "padded unit box: " << boxR.infoString();
        
        Box3d box1(-2, -1, 3, 4, -1, 2);
        XyzVector queryPt(-1.5, 3.75, 0.0);
        std::cout << "test box: " << box1.infoString();
        std::cout << "\taspectRatio = " << box1.aspectRatio() << std::endl;
        std::cout << "\tminRadius = " << box1.minRadius() << std::endl;
        std::cout << "\tmaxRadius = " << box1.maxRadius() << std::endl;
        std::cout << "\tclosestPoint to origin: " << box1.closestPointInBox() << std::endl;
        std::cout << "\tfarthestPoint from origin: " << box1.farthestPointInBox() << std::endl;    
        std::cout << "\tcontains query point " << queryPt << "?  " << (box1.containsPoint(queryPt) ? "yes" : "no") << std::endl;
        std::cout << "\tclosest point to query: " << box1.closestPointInBox(queryPt) << std::endl;
        std::cout << "\tfarthest point from query: " << box1.farthestPointInBox(queryPt) << std::endl;
        std::cout << "\tintersects unit sphere? " << (box1.intersectsSphere() ? "yes" : "no") << std::endl;
        std::cout << "\tintersects sphere centered at origin with radius = 3.5? " << (box1.intersectsSphere(XyzVector(), 3.5) ? "yes" : "no" ) << std::endl;
    }
return 0;
}

