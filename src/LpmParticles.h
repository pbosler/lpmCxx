#ifndef _LPM_PARTICLES_H_
#define _LPM_PARTICLES_H_

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmOutputMessage.h"
#include "LpmLogger.h"
#include "LpmCoords.h"
#include <memory>
#include <vector>
#include <map>

namespace Lpm {

class Particles {
    public :
    
    protected :
        std::map<std::string, std::unique_ptr<Field>> fieldMap;
        std::unique_ptr<Coords> coords;
        std::unique_ptr<Coords> lagCoords;
    
};

}
#endif