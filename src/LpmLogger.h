//
//  Logger.h
//  LPM
//
//  Created by Peter Bosler on 10/31/14.
//  Copyright (c) 2014 Peter Bosler. All rights reserved.
//

#ifndef __LPM__Logger__
#define __LPM__Logger__
/**	@brief Logger header 
	@file
	@author Peter Bosler <pabosle@sandia.gov>
 */
#include "LpmTypeDefs.h"
#include "LpmConfig.h"
#include <iostream>
#include "LpmOutputMessage.h"
#include <vector>
#include <string>

namespace Lpm {

/**
 @class Logger
 @brief Logger class handles console output.  All clients must access Logger via the Logger::Instance member function.
 
 It suppresses output below a given priority level, and consolidates 
 output to process zero for all message below a second (higher) priority level.  
 Error messages are output from all processes.  
 */
class Logger
{
public:
    void logMessage( const OutputMessage msg) const;
    
    void logMessage( const LongMessage longMsg ) const;
    
    Logger( const OutputMessage::priority logLevel = OutputMessage::debugPriority, const std::string logid = "", const int procRank = 0);
    
    inline OutputMessage::priority getPriorityLevel() const {return _baseLevel;}
    inline OutputMessage::priority getAllOutputPriorityLevel() const {return _allOutputLevel;}
    
    inline void setProcRank(const int rank) {_procRank = rank;}
    
    inline void startSection() {_tablevel += 1;}
    void endSection() {_tablevel = std::max(_tablevel - 1, 0);}
private:
    OutputMessage::priority _baseLevel;
    OutputMessage::priority _allOutputLevel;
    std::string _key;
    int _tablevel;
    int _procRank;
};

}
#endif /* defined(__LPM__Logger__) */
