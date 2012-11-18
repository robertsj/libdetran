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
  , d_source(d_number_groups)
{
  // Preconditions
  Require(d_quadrature);

  // Initialize the source vector
  for (size_t g = 0; g < d_number_groups; ++g)
  {
    d_source[g].resize(d_quadrature->number_angles(),
                       vec_dbl(d_mesh->number_cells(), 0));
  }
}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file SyntheticDiscreteSource.cc
//---------------------------------------------------------------------------//
