#ifndef _LPM_TRI_FACES_H_
#define _LPM_TRI_FACES_H_

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmFaces.h"

namespace Lpm {

class TriFaces : public Faces {
    public:
        typedef std::tuple<index_type, index_type, index_type, index_type> quad_index_type;
        
        TriFaces(const index_type nMax, const std::shared_ptr<Edges> edge_ptr, const std::shared_ptr<Coords> crd_ptr,
            const std::shared_ptr<Coords> lag_crd_ptr = 0, const bool sim3d = false);
            
        void divide(const index_type i);
        
};

}

#endif
