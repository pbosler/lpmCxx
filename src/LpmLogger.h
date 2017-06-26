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
#include <iostream>
#include "OutputMessage.h"
#include <vector>
#include <string>

/**
 @class Logger
 @brief Logger class handles console output.  All clients must access Logger via the Logger::Instance member function.
 
 It suppresses output below a given priority level, and consolidates 
 output to process zero for all message below a second (higher) priority level.  
 Error messages are output from all processes.  
 
 Implemented as a singleton pattern.  All clients must access the Logger via the @ref Logger::Instance() method.
 */
class Logger
{
public:
    static Logger* Instance( const OutputMessage::priority level = OutputMessage::debugPriority,
                             const int procRank = 0, const int numProcs = 1);
    void logMessage( const OutputMessage msg) const;
    
    void logMessage( const LongMessage longMsg ) const;
    
protected:
    Logger( const OutputMessage::priority logLevel = OutputMessage::debugPriority,
           const int procRank = 0, const int numProcs = 1);
    
private:
    OutputMessage::priority _baseLevel;
    OutputMessage::priority _allOutputLevel;
    int _procRank;
    int _numProcs;
    static Logger* _instance;
};

#endif /* defined(__LPM__Logger__) */
