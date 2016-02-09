//
//  OutputMessage.cpp
//  LPM
//
//  Created by Peter Bosler on 10/31/14.
//  Copyright (c) 2014 Peter Bosler. All rights reserved.
//
/**	@brief OutputMessage implementation
	@file
	@author Peter Bosler <pabosle@sandia.gov>
 */
#include "OutputMessage.h"

/** @brief Main constructor.
    @param msg
 @param pri
 @param codeOrigin
 */
OutputMessage::OutputMessage( const std::string msg, const priority pri, const std::string codeOrigin)
{
    _msgString = msg;
    _pri = pri;
    _codeOrigin = codeOrigin;
};

/** @brief returns a string that indicates a message's priority level
    @return priString
 */
std::string OutputMessage::priorityString() const
{
    std::string priString;
    switch (OutputMessage::_pri)
    {
        case OutputMessage::debugPriority:
            priString = "debugMessage :";
            break;
        case OutputMessage::tracePriority:
            priString = "traceMessage :";
            break;
        case OutputMessage::remarkPriority:
            priString = "Remark :";
            break;
        case OutputMessage::warningPriority:
            priString = "WARNING : ";
            break;
        case OutputMessage::errorPriority:
            priString = "ERROR : ";
            break;
        default:
            priString = "undefined priority :";
            break;
    }
    return priString;
};

/** @brief Prints an OutputMessage to an output stream
    @param os
 */
void OutputMessage::printMsg(std::ostream& os) const
{
    os << "          msg : " << _msgString << std::endl;
    os << " code Location: " << _codeOrigin << std::endl;
};

/** @brief adds a string to the content of a LongMessage
	@param s
*/
void LongMessage::addStringToMessage( const std::string s)
{
	_msgStrings.push_back(s);
};

/** @brief Prints a LongMessage to an output stream
	@param os
*/
void LongMessage::printMsg( std::ostream& os ) const
{
	os << "      msg title : " << _msgString << std::endl;
	os << "  code location : " << _codeOrigin << std::endl;
	for ( int i = 0; i < _msgStrings.size(); ++i)
	{
		os << _msgStrings[i] << std::endl;
	}
};