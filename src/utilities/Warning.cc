/*
 * Warning.cc
 *
 *  Created on: Mar 20, 2012
 *      Author: robertsj
 */

#include "Warning.hh"

namespace detran_utils
{

#ifdef DETRAN_ENABLE_DEBUG

void warning(int type, std::string message)
{
  std::cout << "WARNING::" << type << " " << message << std::endl;
}

#else


#endif

} // end namespace detran_utils
