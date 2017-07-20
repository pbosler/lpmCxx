#ifndef _LPM_QUAD_FACES_H_
#define _LPM_QUAD_FACES_H_

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmEdges.h"
#include "LpmCoords.h"
#include "LpmFaces.h"
#include <tuple>
#include <memory>

namespace Lpm {

class QuadFaces : public Faces {
    public:
        QuadFaces(const index_type nMax, const std::shared_ptr<Edges> edge_ptr, const std::shared_ptr<Coords> crd_ptr,
            const std::shared_ptr<Coords> lag_crd_ptr = 0, const bool sim3d = false);
            
        void divide(const index_type i);
         
};

}

#endif
