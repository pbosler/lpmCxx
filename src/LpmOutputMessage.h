//
//  OutputMessage.h
//  LPM
//
//  Created by Peter Bosler on 10/31/14.
//  Copyright (c) 2014 Peter Bosler. All rights reserved.
//

#ifndef __LPM__OutputMessage__
#define __LPM__OutputMessage__
/** @file OutputMessage.h
	@brief OutputMessage header file. Defines OutputMessage and LongMessage classes.
	@author Peter Bosler <pabosle@sandia.gov>
 */
#include "LpmTypeDefs.h"
#include "LpmConfig.h"
#include <string>
#include <iostream>
#include <vector>

namespace Lpm {

/**
  @class OutputMessage 
  @brief A formatting class for handling console output. @n
  A message is defined as having __content__, a __priority level__, and an __origin__ (the code 
    location where the message originates).
  
  Messages are divided into 5 priority levels : @n
    @ref debugPriority = a debugging message, used for developing new code. @n
    @ref tracePriority = a message whose purpose is to state its location in the code. 
        Typically a trace message with have content like "entering <class>::<function>" to
        help debugging or progress tracking in an application. @n
    @ref  remarkPriority = a message for general output, often to indicate context or state of a program. @n
    @ref warningPriority = a message that indicates a problem that may not require aborting the program. @n
    @ref errorPriority = a message that indicates a problem that likely requires an abort.
 
  OutputMessage is used by Logger class to handle console output in a parallel environment.
*/
class OutputMessage
{
public:
    enum priority { debugPriority,   ///< indicates a message used for debugging
                    tracePriority,   ///< a message used to report location in code
                    remarkPriority,  ///< a general remark
                    warningPriority, ///< a warning message
                    errorPriority    ///< an error message
                    };
    /// constructor
    OutputMessage( const std::string msg = " ", const priority pri = debugPriority, const std::string codeOrigin = " " );
    /// destructor
	virtual ~OutputMessage(){};
    
    std::string priorityString() const;

    inline std::string getMessage() const {return _msgString;}
    inline void resetMsgString(const std::string& newstr) { _msgString = newstr;}
    
    inline priority getPriority() const { return _pri; };
    
    void printMsg(std::ostream& os) const;
    
protected:
    std::string _msgString;
    priority _pri;
    std::string _codeOrigin;
};

std::string priorityString(const OutputMessage::priority);

/// basic output
inline std::ostream& operator << ( std::ostream& os, const OutputMessage& msg)
{
    os << msg.priorityString() << std::endl;
    msg.printMsg(os);
    return os;
};


/** @class LongMessage
	@brief An extension of OutputMessage to enable messages whose content is best represented in multiple strings.
*/
class LongMessage : public OutputMessage
{
public:
	/** @brief Constructor.
		@param msgTitle
		@param pri
		@param codeOrigin
		@param content
	*/
	LongMessage( const std::string msgTitle = " ", const priority pri = debugPriority, 
				 const std::string codeOrigin = " ", const std::vector< std::string > content = std::vector<std::string>() ) :
				 OutputMessage( msgTitle, pri, codeOrigin ), _msgStrings(content) {};
	
	void printMsg( std::ostream & os ) const;
	
	void addStringToMessage( const std::string s );
	
	/// returns the number of strings in a LongMessage object's content _msgStrings
	int numStrings() const { return _msgStrings.size(); }
	
private:
	std::vector< std::string > _msgStrings;	
};

/// basic output
inline std::ostream& operator << ( std::ostream& os, const LongMessage& msg)
{
    os << msg.priorityString() << std::endl;
    msg.printMsg(os);
    return os;
};

} //namespace
#endif /* defined(__LPM__OutputMessage__) */
