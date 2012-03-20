//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Definitions.hh
 * \author Jeremy Roberts
 * \brief  Useful type definitions, etc.
 */
//---------------------------------------------------------------------------//

#ifndef DETRAN_UTILS_DEFINITIONS_HH_
#define DETRAN_UTILS_DEFINITIONS_HH_

#include <vector>

namespace detran_utils
{

typedef std::vector<int>      vec_int;
typedef std::vector<double>   vec_dbl;

typedef std::vector<vec_int>  vec2_int;
typedef std::vector<vec_dbl>  vec2_dbl;

typedef std::vector<vec2_int> vec3_int;
typedef std::vector<vec2_dbl> vec3_dbl;

} // end namespace detran_utils

#endif // DETRAN_UTILS_DEFINITIONS_HH_
