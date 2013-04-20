//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   EquationTypeTraits.hh
 *  @brief  EquationTypeTraits
 *  @author Jeremy Roberts
 *  @date   Jan 14, 2013
 */
//---------------------------------------------------------------------------//

#ifndef detran_EQUATIONTYPETRAITS_HH_
#define detran_EQUATIONTYPETRAITS_HH_

#include "DimensionTraits.hh"
#include "Equation_DD_1D.hh"

namespace detran
{

enum EQUATION_TYPES
{
  DD1D, DD2D, DD3D
};

template <int N>
class EquationTypeTraits
{
  /* ... */
};

template <>
class EquationTypeTraits<DD1D>
{
  typedef Equation_DD_1D equation_type;
  typedef _1D            dimension_type;
};


} // end namespace detran

#endif // EQUATIONTYPETRAITS_HH_ 

//---------------------------------------------------------------------------//
//              end of file EquationTypeTraits.hh
//---------------------------------------------------------------------------//
