#ifndef _LPM_VTK_FILE_IO_H_
#define _LPM_VTK_FILE_IO_H_

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmCoords.h"
#include "LpmField.h"
#include "LpmFaces.h"
#include "LpmLogger.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <memory>

namespace Lpm {

class VtkWriter {
    public:
        void writeVTKHeader(std::ofstream& os, const std::string& title = "") const;
        void writeCoordsToVTKPoints(std::ofstream& os, const std::shared_ptr<Coords>& crds) const;
        void writeFacesToVTKPolygons(std::ofstream& os, const std::shared_ptr<Faces>& faces) const;
        
        void writeVTKPointDataHeader(std::ofstream& os, const index_type nPoints) const;
        void writeCoordsToVTKPointData(std::ofstream& os, const std::shared_ptr<Coords>& crds, 
            const bool lagrangian = true) const;
        void writeFieldToVTKData(std::ofstream& os, const std::shared_ptr<Field>& field) const;
        
        void writeVTKCellDataHeader(std::ofstream& os, const index_type nCells) const;
        void writeFaceAreaToVTKCellData(std::ofstream& os, const std::shared_ptr<Faces>& faces) const;
        void writeFaceFieldToVTKCellData(std::ofstream& os, const std::shared_ptr<Faces>& faces, 
            const std::shared_ptr<Field>& field) const;
    protected:
        static std::unique_ptr<Logger> log;
};

}

#endif
