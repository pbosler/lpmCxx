#include "LpmTimer.h"
#include <sstream>

namespace Lpm {

std::string Timer::infoString(const bool verbose) const {
    std::stringstream ss;
    if (verbose) {
        ss << "Timer(" << _name << ") : (start, end, elapsed) = (" << _startTime << ", " << _endTime << ", "
        << _endTime - _startTime << ")" << std::endl;
    }
    else {
        ss << _name << " elapsed time: " << _endTime - _startTime << std::endl;
    }
    return ss.str();
}

}
