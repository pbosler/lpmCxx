#include "LpmMeshSeed.h"
#include "LpmOutputMessage.h"
#include "LpmEuclideanCoords.h"
#include "LpmSphericalCoords.h"
#include "LpmTriFaces.h"
#include "LpmQuadFaces.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <exception>
#include <cmath>

namespace Lpm {

std::unique_ptr<Logger> MeshSeed::log(new Logger(OutputMessage::debugPriority, "MeshSeed_log"));

void MeshSeed::initMeshFromSeed(std::shared_ptr<Coords> crds, std::shared_ptr<Edges> edges, std::shared_ptr<Faces> faces,
    const scalar_type domainRadius) {
    readSeedFile();
    bool errorState = false;
    if (crds->n() > 0) {
        OutputMessage errMsg("coordinates have already been initialized", OutputMessage::errorPriority, "MeshSeed::initMeshFromSeed");
        log->logMessage(errMsg);
        errorState = true;
    }
    if (edges->n() > 0) {
        OutputMessage errMsg("edges have already been initialized", OutputMessage::errorPriority, "MeshSeed::initMeshFromSeed");
        log->logMessage(errMsg);
        errorState = true;
    }
    if (faces->n() > 0) {
        OutputMessage errMsg("faces have already been initialized", OutputMessage::errorPriority, "MeshSeed::initMeshFromSeed");
        log->logMessage(errMsg);
        errorState = true;
    }
    if (errorState) {
        throw std::logic_error("initMeshFromSeed expects empty mesh objects");
    }
    if (crds->nMax() < _nCoords) {
        OutputMessage errMsg("not enough memory in coordinates", OutputMessage::errorPriority, "MeshSeed::initMeshFromSeed");
        log->logMessage(errMsg);
        errorState = true;
    }
    if (edges->nMax() < _nEdges) {
        OutputMessage errMsg("not enough memory in edges", OutputMessage::errorPriority, "MeshSeed::initMeshFromSeed");
        log->logMessage(errMsg);
        errorState = true;
    }
    if (faces->nMax() < _nFaces) {
        OutputMessage errMsg("not enough memory in faces", OutputMessage::errorPriority, "MeshSeed::initMeshFromSeed");
        log->logMessage(errMsg);
        errorState = true;
    }
    if (errorState) {
        throw std::out_of_range("mesh objects insufficient to contain MeshSeed");
    }
    
    for (index_type i = 0; i < vertCrds.size(); ++i) {
        vertCrds[i].scale(domainRadius);
        crds->insert(vertCrds[i]);
    }
    for (index_type i = 0; i < edgeOrigs.size(); ++i)
        edges->insert(edgeOrigs[i], edgeDests[i], edgeLefts[i], edgeRights[i]);
    for (index_type i = 0; i < faceEdges.size(); ++i)
        faces->insert(faceEdges[i]);	
    faces->resetAreas();	
}

std::string MeshSeed::infoString() const {
    std::ostringstream ss;
    std::string fullFilename(LPM_MESH_SEED_DIR);
    fullFilename += "/";
    fullFilename += _fname;
    ss << "MeshSeed info:\n";
    ss << "\tseed file = " << _fname << " (" << fullFilename << ")" << std::endl;
    ss << "\tvertCrds: " << _nCoords << " expected, found " << vertCrds.size() << std::endl;
    for (index_type i = 0; i < vertCrds.size(); ++i) 
        ss << "\t\t" << vertCrds[i] << std::endl;

    ss << "\tedges: " << _nEdges << " expected, found " << edgeOrigs.size() << std::endl;
    for (index_type i = 0; i < edgeOrigs.size(); ++i)
        ss << "\t\t" << edgeOrigs[i] << ", " << edgeDests[i] << ", " << edgeLefts[i] << ", " << edgeRights[i] << std::endl;
        
    ss << "\tfaces: " << _nFaces << " expected, found " << faceEdges.size() << std::endl;
    for (index_type i = 0; i < faceEdges.size(); ++i) {
        ss << "\t\t";
        for (index_type j = 0; j < faceEdges[i].size(); ++j)
            ss << faceEdges[i][j] << " ";
        ss << std::endl;
    }
    return ss.str();
}

void MeshSeed::readSeedFile(){
    vertCrds.clear();
    edgeOrigs.clear();
    edgeDests.clear();
    edgeLefts.clear();
    edgeRights.clear();
    faceEdges.clear();
    std::string line;
    std::string fullFilename(LPM_MESH_SEED_DIR);
    fullFilename += "/";
    fullFilename += _fname;
    std::ifstream seedFile(fullFilename);
    if (!seedFile.is_open()) {
        std::ostringstream ss;
        ss << "cannot open seed file: " << fullFilename;
        OutputMessage errMsg(ss.str(), OutputMessage::errorPriority, "MeshSeed:readSeedFile");
        log->logMessage(errMsg);
        throw std::ios_base::failure("file not found");
    }
    index_type lineNumber = 0;
    index_type edgeSectionHeaderLine = -1;
    index_type faceVertHeaderLine = -1;
    index_type faceSectionHeaderLine = -1;
    index_type vertexDegreeHeaderLine = -1;
    while (std::getline(seedFile, line)) {
        lineNumber += 1;
        if (line.find("edgeO") != std::string::npos) {
            edgeSectionHeaderLine = lineNumber;
        }
        if (line.find("faceverts") != std::string::npos) {
            faceVertHeaderLine = lineNumber;
        }
        if (line.find("faceedges") != std::string::npos) {
            faceSectionHeaderLine = lineNumber;
        }
        if (line.find("vertexdegree") != std::string::npos){
            vertexDegreeHeaderLine = lineNumber;
        }
        
        std::istringstream iss(line);
        
        if (lineNumber > 1 && lineNumber < _nCoords + 2) {
            scalar_type x;
            scalar_type y;
            scalar_type z;
            bool crdError = false;
            switch (_nDim) {
                case (2) : {
                    if (! (iss >> x >> y)) {
                        crdError = true;
                    }
                    vertCrds.push_back(XyzVector(x,y));
                    break;
                }
                case (3) : {
                    if (!(iss >> x >> y >> z)) {
                        crdError = true;
                    }
                    vertCrds.push_back(XyzVector(x,y,z));
                    break;
                }
            }
            if (crdError) {
                std::ostringstream oss;
                oss << "cannot read coordinates from line " << lineNumber;
                OutputMessage errMsg(oss.str(), OutputMessage::errorPriority, "MeshSeed::readSeedFile");
                log->logMessage(errMsg);
            }
        }
        else if (edgeSectionHeaderLine > 0 && lineNumber > edgeSectionHeaderLine && lineNumber < edgeSectionHeaderLine + _nEdges + 1) {
            bool edgeError = false;
            index_type origInd;
            index_type destInd;
            index_type leftInd;
            index_type rightInd;
            if (!(iss >> origInd >> destInd >> leftInd >> rightInd)) {
                edgeError = true;
            }
            
            edgeOrigs.push_back(origInd);
            edgeDests.push_back(destInd);
            edgeLefts.push_back(leftInd);
            edgeRights.push_back(rightInd);
            
            if (edgeError) {
                std::ostringstream oss;
                oss << "cannot read edge record from line " << lineNumber;
                OutputMessage errMsg(oss.str(), OutputMessage::errorPriority, "MeshSeed::readSeedFile");
                log->logMessage(errMsg);
            }
        }
        else if (faceSectionHeaderLine > 0 && lineNumber > faceSectionHeaderLine && lineNumber < faceSectionHeaderLine + _nFaces + 1) {
            bool faceError = false;
            index_type edge0Ind;
            index_type edge1Ind;
            index_type edge2Ind;
            index_type edge3Ind;
            std::vector<index_type> edgeList;
            switch (_nEdgesPerFace) {
                case (3) : {
                    if (!(iss >> edge0Ind >> edge1Ind >> edge2Ind)) {
                        faceError = true;
                    }
                    edgeList = {edge0Ind, edge1Ind, edge2Ind};
                    break;
                }
                case (4) : {
                    if (!(iss >> edge0Ind >> edge1Ind >> edge2Ind >> edge3Ind)) {
                        faceError = true;
                    }
                    edgeList = {edge0Ind, edge1Ind, edge2Ind, edge3Ind};
                    break;
                }
            }
            
            faceEdges.push_back(edgeList);
            
            if (faceError) {
                std::ostringstream oss;
                oss << "cannot read face record from line " << lineNumber;
                OutputMessage errMsg(oss.str(), OutputMessage::errorPriority, "MeshSeed::readSeedFile");
                log->logMessage(errMsg);
            }
        }
    }
// #ifdef DEBUG_ALL
//     std::cout << "edgeSectionHeaderLine = " << edgeSectionHeaderLine << std::endl;
//     std::cout << "faceSectionHeaderLine = " << faceSectionHeaderLine << std::endl;
// #endif
}

index_type TriHexSeed::nFaces(const int recursionLevel) const {
    return 6 * std::pow(4, recursionLevel);
}

index_type QuadRectSeed::nFaces(const int recursionLevel) const {
    return 4 * std::pow(4, recursionLevel);
}

index_type IcosTriSphereSeed::nFaces(const int recursionLevel) const {
    return 20 * std::pow(4, recursionLevel);
}

index_type CubedSphereSeed::nFaces(const int recursionLevel) const {
    return 6 * std::pow(4, recursionLevel);
}

index_type TriHexSeed::nVertices(const int recursionLevel) const {
    index_type result = 0;
    for (int i = std::pow(2, recursionLevel) + 1; i <= std::pow(2, recursionLevel +1); ++i)
    {
        result += i;
    }
    result *= 2;
    result += std::pow(2, recursionLevel+1);
    result += 1;
    return result;
}

index_type QuadRectSeed::nVertices(const int recursionLevel) const {
    index_type result = 3;
    for (int i = 1; i <= recursionLevel; ++i) {
        result += std::pow(2, i);
    }
    result *= result;
    return result;
}

index_type IcosTriSphereSeed::nVertices(const int recursionLevel) const {
    return 2 + 10 * std::pow(4, recursionLevel);
}

index_type CubedSphereSeed::nVertices(const int recursionLevel) const {
    return 2 + 6 * std::pow(4, recursionLevel);
}

index_type TriHexSeed::nEdges(const index_type nverts, const index_type nfaces) const {
    return nfaces + nverts - 1;
}

index_type QuadRectSeed::nEdges(const index_type nverts, const index_type nfaces) const {
    return nfaces + nverts - 1;
}

index_type IcosTriSphereSeed::nEdges(const index_type nverts, const index_type nfaces) const {
    return nfaces + nverts - 2;
}

index_type CubedSphereSeed::nEdges(const index_type nverts, const index_type nfaces) const {
    return nfaces + nverts - 2;
}

}
