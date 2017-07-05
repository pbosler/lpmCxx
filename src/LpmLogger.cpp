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
#include <sstream>

using std::cout;
using std::cerr;
using std::endl;

namespace Lpm {

/** @brief constructor
 @param logLevel
 @param procRank
 @param numProcs
*/
Logger::Logger( const OutputMessage::priority logLevel, const std::string logid, const int procRank) :
_baseLevel(logLevel), _allOutputLevel(OutputMessage::warningPriority), _key(logid), _tablevel(0), _procRank(procRank) {};


/**
 @brief Checks an OutputMessage for its priority level, outputs to console from process 0 if priority level is met.
        Outputs error messages from all processes.
 @param msg
 */
void Logger::logMessage( const OutputMessage msg) const
{
    std::stringstream ss;
    for (int i = 0; i < _tablevel; ++i) 
        ss << "\t";
    const std::string tabstr = ss.str();
    std::string formattedMessage(tabstr);
    std::string::size_type found_pos;
    std::string::size_type start_pos = 0;
    const std::string msgstr = msg.getMessage();
    for (index_type i = 0; i < msgstr.size(); ++i) {
        if (msgstr[i] == '\n') {
            formattedMessage += msgstr[i];
            formattedMessage += tabstr;
        }
        else
            formattedMessage += msgstr[i];
    }
    
    if ( msg.getPriority() >= _baseLevel && msg.getPriority() < _allOutputLevel) {
        if ( _procRank == 0 && msg.getPriority() >= _baseLevel) {
            cout << tabstr << "------- start message --------" << std::endl;
            cout << tabstr << _key << " proc_" << _procRank << " (" << msg.priorityString() << ")" << std::endl;
            cout << tabstr << formattedMessage << std::endl;
            cout << tabstr << "------- end message ----------" << std::endl;
        }
    }
    else {
        cout << _key << " proc_" << _procRank << " :" << formattedMessage << endl;
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