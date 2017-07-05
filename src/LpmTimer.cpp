#include "LpmTimer.h"
#include <sstream>

namespace Lpm {

std::string Timer::infoString() const {
    std::stringstream ss;
    ss << "Timer(" << _name << ") : (start, end, elapsed) = (" << _startTime << ", " << _endTime << ", "
       << _endTime - _startTime << ")" << std::endl;
    return ss.str();
}

}
