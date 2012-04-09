/*
 * warning.hh
 *
 *  Created on: Mar 20, 2012
 *      Author: robertsj
 */

#ifndef WARNING_HH_
#define WARNING_HH_

#include <iostream>
#include <exception>
#include <string>
#include "detran_config.h"

namespace detran
{

enum WARNING
{
  USER_INPUT,
  SOLVER_CONVERGENCE,
  END_WARN
};

#ifdef DETRAN_ENABLE_DEBUG

void warning(int type, std::string message);

#else

inline void warning(int type, std::string message){}

#endif

} // end namespace detran

#endif /* WARNING_HH_ */
