//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   GenException.cc
 *  @author Jeremy Roberts
 *  @date   04/09/2011
 *  @brief  Member definitions of class GenException
 *  @note   Modified version of K. Huff's class from cyclus
 */
//---------------------------------------------------------------------------//

#include <sstream>
#include "GenException.hh"

namespace detran_utilities
{

//---------------------------------------------------------------------------//
std::string itoa(int i)
{
  std::stringstream out;
  out << i;
  return out.str();
}

//---------------------------------------------------------------------------//
std::string dtoa(double d)
{
  std::stringstream out;
  out << d;
  return out.str();
}

std::string GenException::prepend = "detran exception";

//---------------------------------------------------------------------------//
GenException::GenException()
{
	myMessage = prepend;
}

//---------------------------------------------------------------------------//
GenException::GenException(int line, std::string file, std::string msg)
{
	myMessage = prepend + "\n" 
              + "           on line: " + itoa(line) + "\n"
              + "           in file: " + file + "\n" 
              + "           message: " + msg;
}

//---------------------------------------------------------------------------//
const char* GenException::what() const throw()
{
	return myMessage.c_str();
}

//---------------------------------------------------------------------------//
GenException::~GenException() throw() 
{
  /* ... */
}

} // end namespace detran_utilities
