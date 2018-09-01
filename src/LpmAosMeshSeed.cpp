#include "LpmAosMeshSeed.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <exception>

namespace Lpm {
namespace Aos {

std::string MeshSeed::infoString() const {
    std::ostringstream ss;
    std::cout << "MESH SEED INFO: id = " << this->idString() << std::endl;
    std::cout << "\tseed file = " << _fname << std::endl;
    if (_ndim == 2) {
        std::cout << "\tr2vertexCrds :" << std::endl;
        for (index_type i=0; i<_nCrds; ++i) {
            std::cout  << "\t" << r2vertexCrds[i] << std::endl;
        }
        std::cout << "\tr3vertexCrds.empty() = " << (r3vertexCrds.empty() ? "true" : "false") << std::endl;
    }
    else if (_ndim == 3) {
        std::cout << "\tr2vertexCrds.empty() = " << (r2vertexCrds.empty() ? "true" : "false") << std::endl;
        std::cout << "\tr3vertexCrds :" << std::endl;
        for (index_type i=0; i<_nCrds; ++i) {
            std::cout << "\t" << r3vertexCrds[i] << std::endl;
        }        
    }
    std::cout << "\tedges :" << std::endl;
    std::cout << std::setw(20) << "orig" << std::setw(20) <<  "dest"<< std::setw(20)  << "left"
        << std::setw(20) << "right" << std::endl;
    for (index_type i=0; i<_nEdges; ++i) {
        std::cout << std::setw(20) << edgeOrigs[i] << std::setw(20) << edgeDests[i] << std::setw(20) << edgeLefts[i]
            << std::setw(20) << edgeRights[i] << std::endl;
    }
    std::cout << "\tface vertices:" << std::endl;
    for (index_type i=0; i<_nFaces; ++i) {
        std::cout << "\t";
        for (index_type j=0; j<_nEdgesPerFace; ++j) {
            std::cout << faceVerts[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "\tface edges:" << std::endl;
    for (index_type i=0; i<_nFaces; ++i) {
        std::cout << "\t";
        for (index_type j=0; j<_nEdgesPerFace; ++j) {
            std::cout << faceEdges[i][j] << " ";
        }
        std::cout << std::endl;
    }
    return ss.str();
}

void MeshSeed::determineMaxAllocations(index_type& nv, index_type& nf, index_type& ne, const index_type maxRec) const {
    nv = this->nVerticesAfterUniformRefinement(maxRec);
    nf = 0;
    ne = 0;
    for (int k=0; k<=maxRec; ++k) {
        nf += this->nFacesAfterUniformRefinement(k);
        ne += this->nEdgesAfterUniformRefinement(nVerticesAfterUniformRefinement(k), nFacesAfterUniformRefinement(k));
    }
}

void MeshSeed::initFromFile() {
    // start with clean slate
    r2vertexCrds.clear();
    r3vertexCrds.clear();
    edgeOrigs.clear();
    edgeDests.clear();
    edgeLefts.clear();
    edgeRights.clear();
    faceVerts.clear();
    faceEdges.clear();
    
    // open mesh seed data file
    std::string fullFilename(LPM_MESH_SEED_DIR);
    fullFilename += "/";
    fullFilename += _fname;
    _fname = fullFilename;
    
    
    std::ifstream seedFile(fullFilename);
    if (!seedFile.is_open()) {
        std::ostringstream ss;
        ss << "error: cannot open seed file " << fullFilename;
        throw std::ios_base::failure(ss.str());
    }
    else {
        std::cout << "...reading mesh seed data from file " << fullFilename << std::endl;
    }
    // parse file
    std::string line;
    index_type lineNumber=0;
    index_type edgeSectionHeaderLine = -1;
    index_type faceVertHeaderLine = -1;
    index_type faceSectionHeaderLine = -1;
    index_type vertexDegreeHeaderLine = -1;
    while (std::getline(seedFile, line)) {
        ++lineNumber;
        if (line.find("edgeO") != std::string::npos) {
            edgeSectionHeaderLine = lineNumber;
//std::cout << edgeSectionHeaderLine << std::endl;
        }
        if (line.find("faceverts") != std::string::npos) {
            faceVertHeaderLine = lineNumber;
//             std::cout << faceVertHeaderLine << std::endl;
        }
        if (line.find("faceedges") != std::string::npos) {
            faceSectionHeaderLine = lineNumber;
//             std::cout << faceSectionHeaderLine << std::endl;
        }
        if (line.find("vertexdegree")!= std::string::npos) {
            vertexDegreeHeaderLine = lineNumber;
//            std::cout << vertexDegreeHeaderLine << std::endl;
        }
        
        std::istringstream iss(line);
        if (lineNumber>1 && lineNumber < _nCrds+2) {
            scalar_type x,y,z;
            bool crdErr = false;
            switch (_ndim) {
                case (2) : {
                    if (!(iss >> x >> y)) {
                        crdErr = true;
                    }
                    //std::cout << "(x,y) = (" << x << "," << y << ")" << std::endl;
                    r2vertexCrds.push_back(Vec<2>(x,y));
                    break;
                }
                case (3) : {
                    if (!(iss >> x >> y >> z)) {
                        crdErr = true;
                    }
                    r3vertexCrds.push_back(Vec<3>(x,y,z));
                    break;
                }
                if (crdErr) {
                    std::ostringstream ss;
                    ss << "error: cannot read coordinate from line " << lineNumber << " of file " << fullFilename;
                    throw std::ios_base::failure(ss.str());
                }
            }
        }
        else if (edgeSectionHeaderLine>0 && lineNumber > edgeSectionHeaderLine && 
            lineNumber < edgeSectionHeaderLine + _nEdges + 1) {
            bool edgeErr = false;
            index_type origInd;
            index_type destInd;
            index_type leftInd;
            index_type rightInd;
            if (!(iss >> origInd >> destInd >> leftInd >> rightInd)) {
                edgeErr = true;
            }
            //std::cout << origInd << " " << destInd << " " << leftInd << " " << rightInd << std::endl;
            edgeOrigs.push_back(origInd);
            edgeDests.push_back(destInd);
            edgeLefts.push_back(leftInd);
            edgeRights.push_back(rightInd);
            
            if (edgeErr) {
                std::ostringstream ss;
                ss << "error: cannot read edge from line " << lineNumber << " of file " << fullFilename;
                throw std::ios_base::failure(ss.str());
            }
        }
        else if (faceVertHeaderLine>0 && lineNumber > faceVertHeaderLine && 
            lineNumber < faceVertHeaderLine + _nFaces + 1) {
            bool faceErr = false;
            index_type v0;
            index_type v1;
            index_type v2;
            index_type v3;
            ind_vec_type vertlist;
            switch (_nEdgesPerFace){
                case (3) : {
                    if (!(iss >> v0 >> v1 >> v2)) {
                        faceErr = true;
                    }
                    vertlist = {v0, v1, v2};
                    //std::cout << vertlist[0] << " " << vertlist[1] << " " << vertlist[2] << std::endl;
                    break;
                }
                case (4): {
                    if (!(iss >> v0 >> v1 >> v2 >> v3)) {
                        faceErr = true;
                    }
                    vertlist = {v0, v1, v2, v3};
                    break;
                }
            }
            faceVerts.push_back(vertlist);
            //std::cout << "faceVerts.size() = " << faceVerts.size() << std::endl;
            if (faceErr){
                std::ostringstream ss;
                ss << "error: cannot read face vertices from line " << lineNumber << " of file " << fullFilename;
                throw std::ios_base::failure(ss.str());
            }
       }
       else if (faceSectionHeaderLine>0 && lineNumber > faceSectionHeaderLine && 
        lineNumber < faceSectionHeaderLine + _nFaces +1) {
            bool faceErr = false;
            index_type e0;
            index_type e1;
            index_type e2;
            index_type e3;
            ind_vec_type edgelist;
            switch (_nEdgesPerFace) {
                case (3) : {
                    if (!(iss >> e0 >> e1 >> e2)) {
                        faceErr = true;
                    }
                    edgelist = {e0, e1, e2};
                    break;
                }
                case (4) : {
                    if (!(iss >> e0 >> e1 >> e2 >> e3)) {
                        faceErr = true;
                    }
                    edgelist = {e0, e1, e2, e3};
                    break;
                }
            }
            //std::cout << edgelist[0] << " " << edgelist[1] << " " << edgelist[2] << std::endl;
            faceEdges.push_back(edgelist);
//             std::cout << "faceEdges.size() = " << faceEdges.size() << std::endl;
            if (faceErr) {
                std::ostringstream ss;
                ss << "error: cannot read face edges from line " << lineNumber << " of file " << fullFilename;
                throw std::ios_base::failure(ss.str());
            }
       }
           
    }
}

index_type TriHexSeed::nFacesAfterUniformRefinement(const index_type recursionLevel) const {
    return 6 * std::pow(4, recursionLevel);
}

index_type QuadRectSeed::nFacesAfterUniformRefinement(const index_type recursionLevel) const {
    return 4 * std::pow(4, recursionLevel);
}

index_type IcosTriSphereSeed::nFacesAfterUniformRefinement(const index_type recursionLevel) const {
    return 20 * std::pow(4, recursionLevel);
}

index_type CubedSphereSeed::nFacesAfterUniformRefinement(const index_type recursionLevel) const {
    return 6 * std::pow(4, recursionLevel);
}

index_type TriHexSeed::nVerticesAfterUniformRefinement(const index_type recursionLevel) const {
    index_type result = 0;
    for (index_type i = std::pow(2, recursionLevel) + 1; i <= std::pow(2, recursionLevel +1); ++i) {
        result += i;
    }
    result *= 2;
    result += std::pow(2, recursionLevel+1);
    result += 1;
    return result;
}

index_type QuadRectSeed::nVerticesAfterUniformRefinement(const index_type recursionLevel) const {
    index_type result = 3;
    for (index_type i = 1; i <= recursionLevel; ++i) {
        result += std::pow(2, i);
    }
    result *= result;
    return result;
}

index_type IcosTriSphereSeed::nVerticesAfterUniformRefinement(const index_type recursionLevel) const {
    return 2 + 10 * std::pow(4, recursionLevel);
}

index_type CubedSphereSeed::nVerticesAfterUniformRefinement(const index_type recursionLevel) const {
    return 2 + 6 * std::pow(4, recursionLevel);
}

index_type TriHexSeed::nEdgesAfterUniformRefinement(const index_type nverts, const index_type nfaces) const {
    return nfaces + nverts - 1;
}

index_type QuadRectSeed::nEdgesAfterUniformRefinement(const index_type nverts, const index_type nfaces) const {
    return nfaces + nverts - 1;
}

index_type IcosTriSphereSeed::nEdgesAfterUniformRefinement(const index_type nverts, const index_type nfaces) const {
    return nfaces + nverts - 2;
}

index_type CubedSphereSeed::nEdgesAfterUniformRefinement(const index_type nverts, const index_type nfaces) const {
    return nfaces + nverts - 2;
}

}
}