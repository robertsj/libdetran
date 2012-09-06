//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   GenException.hh
 * \author Jeremy Roberts
 * \date   04/09/2011
 * \brief  GenException class definition.
 * \note   Modified version of K. Huff's class from cyclus
 */
//---------------------------------------------------------------------------//

#ifndef GENEXCEPTION_HH
#define GENEXCEPTION_HH

#include <iostream>
#include <exception>
#include <string>

namespace detran_utilities
{

/**
 *  A generic mechanism to manually manage exceptions
 */
class GenException: public std::exception
{

protected:

    /// The message associated with this exception.
    std::string myMessage;
    
    /// A string to prepend to all message of this class.
    static std::string prepend;
    
public:
    
    /// Constructs a new GenException with the default message.
    GenException();
    
    /**
     * @brief Constructs a new GenException with a provided message
     *
     * @param line line of code erring
     * @param file file in which error occurs
     * @param msg the message
     */
    GenException(int line, std::string file, std::string msg);
    
    /**
     * Returns the error message associated with this GenException.
     *
     * @return the message
     */
    virtual const char* what() const throw();
    
    /**
     * Destroys this GenException.
     */
    virtual ~GenException() throw();
    
};

/// Easy macro for throwing exceptions.
#define THROW(m) throw detran_utilities::GenException(__LINE__,__FILE__,m);

} // end namespace detran_utilities

#endif // GENEXCEPTION_HH
