//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   GenException.hh
 *  @author Jeremy Roberts
 *  @date   04/09/2011
 *  @brief  GenException class definition.
 *  @note   Modified version of K. Huff's class from cyclus
 */
//---------------------------------------------------------------------------//

#ifndef detran_utilities_GENEXCEPTION_HH
#define detran_utilities_GENEXCEPTION_HH

#include "utilities/utilities_export.hh"
#include "utilities/Definitions.hh"
#include <iostream>
#include <exception>
#include <string>

namespace detran_utilities
{

/**
 *  A generic mechanism to manually manage exceptions
 */
class UTILITIES_EXPORT GenException: public std::exception
{

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
  virtual const char* what() const throw ();

  /**
   * Destroys this GenException.
   */
  virtual ~GenException() throw ();


protected:

  /// The message associated with this exception.
  std::string myMessage;
  /// A string to prepend to all message of this class.
  static std::string prepend;

};

/// Easy macro for throwing exceptions.
#define THROW(m) throw detran_utilities::GenException(__LINE__,__FILE__,m);

} // end namespace detran_utilities

#endif // detran_utilities_GENEXCEPTION_HH
