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

namespace detran
{

typedef std::vector<int>      vec_int;
typedef std::vector<double>   vec_dbl;

typedef std::vector<std::vector<int> >  vec2_int;
typedef std::vector<std::vector<double> >  vec2_dbl;

typedef std::vector<std::vector<std::vector<int> > > vec3_int;
typedef std::vector<std::vector<std::vector<double> > > vec3_dbl;

typedef unsigned int u_int;

} // end namespace detran


#endif // DETRAN_UTILS_DEFINITIONS_HH_
