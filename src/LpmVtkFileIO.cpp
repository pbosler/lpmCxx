#include "LpmVtkFileIO.h"
#include "LpmXyzVector.h"
#include "LpmTriFaces.h"
#include "LpmQuadFaces.h"
#include "LpmOutputMessage.h"
#include <limits>

namespace Lpm {

std::unique_ptr<Logger> log(new Logger(OutputMessage::debugPriority, "VtkFileIO_log"));

void VtkWriter::writeVTKHeader(std::ofstream& os, const std::string& title) const {
    os << "# vtk DataFile Version 2.0" << std::endl;
    os << title << std::endl;
    os << "ASCII" << std::endl;
    os << "DATASET POLYDATA" << std::endl;
}

void VtkWriter::writeCoordsToVTKPoints(std::ofstream& os, const std::shared_ptr<Coords>& crds) const {
    os << "POINTS " << crds->n() << " double" << std::endl;
    for (index_type i = 0; i < crds->n(); ++i) {
        const XyzVector ptVec = crds->getVec(i);
        os << std::setprecision(std::numeric_limits<long double>::digits10 + 1);
        os << ptVec.x << " " << ptVec.y << " " << ptVec.z << std::endl;
    }
}

void VtkWriter::writeFacesToVTKPolygons(std::ofstream& os, const std::shared_ptr<Faces>& faces) const {
    index_type cellListSize = 0;
    TriFaces* tri_ptr = dynamic_cast<TriFaces*>(faces.get());
    QuadFaces* quad_ptr = dynamic_cast<QuadFaces*>(faces.get());
    if (tri_ptr) {
        cellListSize = 4 * faces->nLeaves();
    }
    else if (quad_ptr) {
        cellListSize = 5 * faces->nLeaves();
    }
    else {
        for (index_type i = 0; i < faces->n(); ++i) {
            if (!faces->hasChildren(i)) {
                cellListSize += (1 + faces->nVerticesAtFace(i));
            }
        }
    }
    os << "POLYGONS " << faces->nLeaves() << " " << cellListSize << std::endl;
    for (index_type i = 0; i < faces->n(); ++i) {
        if (!faces->hasChildren(i)) {
            const std::vector<index_type> vertInds = faces->vertexIndices(i);
            os << vertInds.size() << " ";
            for (int j = 0; j < vertInds.size(); ++j)
                os << vertInds[j] << " ";
            os << std::endl;
        }
    }
}

void VtkWriter::writeVTKPointDataHeader(std::ofstream& os, const index_type nPoints) const {
    os << "POINT_DATA " << nPoints << std::endl;
}

void VtkWriter::writeCoordsToVTKPointData(std::ofstream& os, const std::shared_ptr<Coords>& crds, 
    const bool lagrangian) const {
    os << "SCALARS " << (lagrangian ? "lagrangian_coords " : "physical_coords ") << "double 3" << std::endl;
    os << "LOOKUP_TABLE default" << std::endl;
    for (index_type i = 0; i < crds->n(); ++i) {
        os << std::setprecision(std::numeric_limits<long double>::digits10 + 1);
        const XyzVector ptVec = crds->getVec(i);
        os << ptVec.x << " " << ptVec.y << " " << ptVec.z << std::endl;
    }
}

void VtkWriter::writeFieldToVTKData(std::ofstream& os, const std::shared_ptr<Field>& field) const {
    std::string fieldstring(field->name());
    fieldstring += "_";
    fieldstring += field->units();
    os << "SCALARS " << fieldstring << " double " << field->nDim() << std::endl;
    os << "LOOKUP_TABLE default" << std::endl;
    os << std::setprecision(std::numeric_limits<long double>::digits10 + 1);
    switch (field->nDim()) {
        case (1) : {
            for (index_type i = 0; i < field->n(); ++i)
                os << field->getScalar(i) << std::endl;
            break;
        }
        case (2) : {
            for (index_type i = 0; i < field->n(); ++i) {
                const XyzVector fVec = field->get2dVector(i);
                os << fVec.x << " " << fVec.y << std::endl;
            }
            break;
        }
        case (3) : {
            for (index_type i = 0; i < field->n(); ++i) {
                const XyzVector fVec = field->get3dVector(i);
                os << fVec.x << " " << fVec.y << " " << fVec.z << std::endl;
            }
            break;
        }
    }
}

void VtkWriter::writeVTKCellDataHeader(std::ofstream& os, const index_type nCells) const {
    os << "CELL_DATA " << nCells << std::endl;
}

void VtkWriter::writeFaceAreaToVTKCellData(std::ofstream& os, const std::shared_ptr<Faces>& faces) const {
    os << "SCALARS faceArea double 1" << std::endl;
    os << "LOOKUP_TABLE default" << std::endl;
    os << std::setprecision(std::numeric_limits<long double>::digits10 + 1);
    for (index_type i = 0; i < faces->n(); ++i) {
        if (!faces->hasChildren(i)) {
            os << faces->area(i) << std::endl;
        }
    }
}

void VtkWriter::writePointsAndFieldsToCSV(std::ofstream& os, const std::shared_ptr<Coords>& crds, 
            const std::vector<std::shared_ptr<Field>>& fields) const {
    os << "x,y,z";
    const bool hasFields = !fields.empty();
    int lastFieldInd = -1;
    if (hasFields) {
        os << ",";
        lastFieldInd = fields.size()-1;
        for (int i=0; i<lastFieldInd; ++i) {
            switch (fields[i]->nDim()) {
                case (1) : {
                    os << fields[i]->name() << ",";
                break;
                }
                case (2) : {
                    os << fields[i]->name() + "0," + fields[i]->name() + "1,";
                break;
                }
                case (3) : {
                    os << fields[i]->name() + "0," + fields[i]->name() + "1," + fields[i]->name() + "2,";
                break;
                }
            }
        }
        switch (fields[lastFieldInd]->nDim()) {
            case (1) : {
                os << fields[lastFieldInd]->name() << std::endl;
            break;
            }
            case (2) : {
                os << fields[lastFieldInd]->name() + "0," + fields[lastFieldInd]->name() + "1" << std::endl;
            break;
            }
            case (3) : {
                os << fields[lastFieldInd]->name() + "0," + fields[lastFieldInd]->name() + "1," +
                      fields[lastFieldInd]->name() + "2"<< std::endl;
            break;
            }
        }
        
    }
    for (index_type i=0; i<crds->n(); ++i) {
        const XyzVector crdvec = crds->getVec(i);
        os << crdvec.x << "," << crdvec.y << "," << crdvec.z;
        if (hasFields) {
            os << ",";
            for (int j=0; j<lastFieldInd; ++j) {
                switch (fields[j]->nDim()) {
                    case (1) : {
                        os << fields[j]->getScalar(i) << ",";
                        break;
                    }
                    case (2) : {
                        const XyzVector fvec = fields[j]->get2dVector(i);
                        os << fvec.x << "," << fvec.y << ",";
                        break;
                    }
                    case (3) : {
                        const XyzVector fvec = fields[j]->get3dVector(i);
                        os << fvec.x << "," << fvec.y << "," << fvec.z << ",";
                        break;
                    }
                }
            }
            switch (fields[lastFieldInd]->nDim()) {
                case (1) : {
                    os << fields[lastFieldInd]->getScalar(i);
                break;}
                case (2) : {
                    const XyzVector fvec = fields[lastFieldInd]->get2dVector(i);
                    os << fvec.x << "," << fvec.y;
                break;}
                case (3) : {
                    const XyzVector fvec = fields[lastFieldInd]->get2dVector(i);
                    os << fvec.x << "," << fvec.y << "," << fvec.z;
                break;}
            }
            
        }
        os << std::endl;
    }
}

void VtkWriter::writeFaceFieldToVTKCellData(std::ofstream& os, const std::shared_ptr<Faces>& faces, 
            const std::shared_ptr<Field>& field) const {
    std::string fieldstring(field->name());
    fieldstring += "_";
    fieldstring += field->units();
    os << "SCALARS " << fieldstring << " double " << field->nDim() << std::endl;
    os << "LOOKUP_TABLE default" << std::endl;
    os << std::setprecision(std::numeric_limits<long double>::digits10 + 1);
    switch (field->nDim()) {
        case (1) : {
            for (index_type i = 0; i < faces->n(); ++i) {
                if (!faces->hasChildren(i)) {
                    os << field->getScalar(i) << std::endl;
                }
            }
            break;
        }
        case (2) : {
            for (index_type i = 0; i < faces->n(); ++i) {
                if (!faces->hasChildren(i)) {
                    const XyzVector fVec = field->get2dVector(i);
                    os << fVec.x << " " << fVec.y << std::endl;
                }
            }
            break;
        }
        case (3) : {
            for (index_type i = 0; i < faces->n(); ++i) {
                if (!faces->hasChildren(i)) {
                    const XyzVector fVec = field->get3dVector(i);
                    os << fVec.x << " " << fVec.y << " " << fVec.z << std::endl;
                }
            }
            break;
        }
    }
}


}
