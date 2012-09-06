//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Warning.hh
 * \author Jeremy Roberts
 * \date   Mar 20, 2012
 * \brief  Warning macros
 */
//---------------------------------------------------------------------------//

#ifndef WARNING_HH_
#define WARNING_HH_

#include "detran_config.h"
#include <iostream>
#include <exception>
#include <string>

namespace detran_utilities
{

enum WARNING
{
  USER_INPUT,
  SOLVER_CONVERGENCE,
  END_WARN
};

#ifdef DETRAN_ENABLE_DEBUG

inline void warning(int type, std::string message)
{
  std::cout << "WARNING::" << type << " " << message << std::endl;
}

#else

inline void warning(int type, std::string message){}

#endif

} // end namespace detran_utilities

#endif /* WARNING_HH_ */
