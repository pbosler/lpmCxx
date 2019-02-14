#include "LpmMPIReplicatedData.h"
#include <sstream>
#include <algorithm>
#include <exception>
#include <string>
#include <sstream>

namespace Lpm {

MPIReplicatedData::MPIReplicatedData(const index_type nItems, const int rank, const int nProcs) :
     _nItems(nItems), _nProcs(nProcs), _rank(rank), _procStartIndex(nProcs, -1), _procEndIndex(nProcs, -1), _procMsgSize(nProcs, -1) {
     loadBalance();
}

void MPIReplicatedData::loadBalance() {
    if (_nItems < _nProcs) {
        std::ostringstream ss;
        ss << "error: cannot distribute " << _nItems << " items over " << _nProcs << " ranks.";
        throw std::logic_error(ss.str());
    }
    const index_type chunkSize = index_type(_nItems / _nProcs);
    for (index_type i = 0; i < _nProcs; ++i) {
        _procStartIndex[i] = i * chunkSize;
        _procEndIndex[i] = (i + 1) * chunkSize - 1;
    }
    _procEndIndex[_nProcs-1] = _nItems - 1;
    
    for (index_type i = 0; i < _nProcs; ++i) {
        _procMsgSize[i] = _procEndIndex[i] - _procStartIndex[i] + 1;
    }
}

std::string MPIReplicatedData::infoString() const {
    std::stringstream ss;
    ss << "MPIReplicatedData info:" << std::endl;
    ss << "\tdistributed " << _nItems << " items over " << _nProcs << " mpi ranks." << std::endl;
    auto minmaxlocs = std::minmax_element(_procMsgSize.begin(), _procMsgSize.end());
    ss << "\tmin nItems per proc = " << *minmaxlocs.first << ", max nItems per proc = " << *minmaxlocs.second << std::endl;
    for (index_type i = 0; i < _nProcs; ++i) {
        ss << "\t\tproc " << i << " (startInd, endInd, msgSize) = (" << _procStartIndex[i] << ", " << _procEndIndex[i]
           << ", " << _procMsgSize[i] << ")" <<std::endl;
    }
    return ss.str();
}


}
