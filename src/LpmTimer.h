#ifndef _LPM_TIMER_H_
#define _LPM_TIMER_H_

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include <mpi.h>
#include <string>

namespace Lpm {

class Timer {
    public :
        Timer(const std::string& name="") : _name(name), _startTime(0.0), _endTime(0.0) {};
        
        inline void start() {_startTime = MPI_Wtime();}
        inline void end() {_endTime = MPI_Wtime();}
        
        std::string infoString() const;
        
        inline scalar_type elapsed() const {return _endTime - _startTime;}
        
    protected:
        scalar_type _startTime;
        scalar_type _endTime;
        std::string _name;
};

}

#endif
