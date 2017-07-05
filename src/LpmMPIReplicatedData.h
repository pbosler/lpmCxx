#ifndef _LPM_MPI_REPLICATED_DATA_H_
#define _LPM_MPI_REPLICATED_DATA_H_

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include <vector>
#include <string>

namespace Lpm {

class MPIReplicatedData {
    public:
        MPIReplicatedData(const index_type nItems, const index_type nProcs = 1);
        
        inline index_type n() const {return _nItems;}
        
        inline index_type startIndex(const index_type procRank) const {return _procStartIndex[procRank];}
        inline index_type endIndex(const index_type procRank) const {return _procEndIndex[procRank];}
        inline index_type msgSize(const index_type procRank) const {return _procMsgSize[procRank];}
        
        std::string infoString() const;
        
        inline void setNItems(const index_type nItems) {_nItems = nItems;}
        void loadBalance();
    
    protected:
        std::vector<index_type> _procStartIndex;
        std::vector<index_type> _procEndIndex;
        std::vector<index_type> _procMsgSize;
        
        index_type _nItems;
        index_type _nProcs;
    
};

}

#endif
