//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   StdOutUtils.hh
 *  @brief  Utilities for printing things to standard output
 *  @author Jeremy Roberts
 *  @date   Nov 8, 2012
 */
//---------------------------------------------------------------------------//

#ifndef detran_ioutils_STDOUTUTILS_HH_
#define detran_ioutils_STDOUTUTILS_HH_

#include "utilities/Definitions.hh"
#include <cstdio>
#include <iostream>
#include <string>

namespace detran_ioutils
{

template <class T>
inline void print_vec(const T &v,
                      const std::string name = "")
{
}

template <>
inline void print_vec(const detran_utilities::vec_dbl &v,
                      const std::string name)
{
  std::cout << " vec_dbl: " << name << std::endl;
  for (detran_utilities::size_t i = 0; i < v.size(); ++i)
    printf(" %4i   %16.8f \n", i, v[i]);
}



} // end namespace detran

#endif // detran_ioutils_STDOUTUTILS_HH_

//---------------------------------------------------------------------------//
//              end of file StdOutUtils.hh
//---------------------------------------------------------------------------//
