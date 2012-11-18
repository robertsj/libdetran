//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   SyntheticDiscreteSource.cc
 *  @brief  SyntheticDiscreteSource
 *  @author Jeremy Roberts
 *  @date   Nov 15, 2012
 */
//---------------------------------------------------------------------------//

#include "SyntheticDiscreteSource.hh"

namespace detran
{

//---------------------------------------------------------------------------//
SyntheticDiscreteSource::SyntheticDiscreteSource(const size_t number_groups,
                                                 SP_mesh mesh,
                                                 SP_quadrature quadrature,
                                                 SP_material material)
  : Base(number_groups, mesh, quadrature, material, true)
{
  /* ... */
}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file SyntheticDiscreteSource.cc
//---------------------------------------------------------------------------//
