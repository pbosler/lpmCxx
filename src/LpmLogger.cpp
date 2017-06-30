//
//  Logger.cpp
//  LPM
//
//  Created by Peter Bosler on 10/31/14.
//  Copyright (c) 2014 Peter Bosler. All rights reserved.
//
/**	@brief Logger implementation
	@file
	@author Peter Bosler <pabosle@sandia.gov>
 */

#include "LpmLogger.h"

using std::cout;
using std::cerr;
using std::endl;

namespace Lpm {

/** @brief constructor
 @param logLevel
 @param procRank
 @param numProcs
*/
Logger::Logger( const OutputMessage::priority logLevel, const std::string logid, const int procRank, const int numProcs) :
_baseLevel(logLevel), _procRank(procRank), _numProcs(numProcs), _allOutputLevel(OutputMessage::warningPriority), 
_key(logid) {};


/**
 @brief Checks an OutputMessage for its priority level, outputs to console from process 0 if priority level is met.
        Outputs error messages from all processes.
 @param msg
 */
void Logger::logMessage( const OutputMessage msg) const
{
    if ( msg.getPriority() < _allOutputLevel) {
        if ( _procRank == 0 && msg.getPriority() >= _baseLevel)
            cout << "------- start message --------" << std::endl;
            cout << _key << " proc0 (" << msg.priorityString() << ")" << std::endl;
            cout << msg.getMessage() << std::endl;
            cout << "------- end message ----------" << std::endl;
    }
    else {
        cout << "proc " << _procRank << " " << _key << ":" << msg << endl;
    }
};

/**
 @brief Checks a LongMessage for its priority level, outputs to console from process 0 if priority level is met.
        Outputs error messages from all processes.
 @param longMsg
 */
void Logger::logMessage( const LongMessage longMsg ) const
{
	if ( longMsg.getPriority() < _allOutputLevel )
	{
		if ( _procRank == 0 && longMsg.getPriority() >= _baseLevel ) {
			cout << "proc 0 " << _key << ": " << longMsg << std::endl;
		}
	}
	else {
		cout << "proc " << _procRank << _key << ": " << longMsg << endl;
	}
}

}