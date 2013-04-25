//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   SyntheticSource.cc
 *  @brief  SyntheticSource
 *  @author Jeremy Roberts
 *  @date   Nov 15, 2012
 */
//---------------------------------------------------------------------------//

#include "SyntheticSource.hh"
#include "BDFCoefficients.hh"

namespace detran
{

//---------------------------------------------------------------------------//
SyntheticSource::SyntheticSource(const size_t   number_groups,
                                 SP_mesh        mesh,
                                 SP_quadrature  quadrature,
                                 SP_material    material,
                                 bool           discrete)
  : Base(number_groups, mesh, quadrature, discrete)
  , d_material(material)
  , d_norm(detran_angle::Quadrature::angular_norm(d_mesh->dimension()))
{
  /* ... */
}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file SyntheticSource.cc
//---------------------------------------------------------------------------//
