//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  IsotropicSource.cc
 *  @brief IsotropicSource member definitions
 *  @note  Copyright (C) Jeremy Roberts 2012-2013
 */
//----------------------------------------------------------------------------//

#include "IsotropicSource.hh"

namespace detran_external_source
{

//----------------------------------------------------------------------------//
IsotropicSource::IsotropicSource(size_t number_groups,
                                 SP_mesh mesh,
                                 spectra_type &spectra,
                                 vec_int &map,
                                 SP_quadrature quadrature)
  : ExternalSource(number_groups, mesh, quadrature)
  , d_source_spectra(spectra)
  , d_source_map(map)
  , d_norm(detran_angle::Quadrature::angular_norm(mesh->dimension()))
{
  Require(d_source_spectra.size()    > 0);
  Require(d_source_spectra[0].size() == d_number_groups);
  Require(d_source_map.size()        == d_mesh->number_cells());

  d_mesh->add_mesh_map("ISOTROPICSOURCE", d_source_map);
}

} // end namespace detran_external_source

//----------------------------------------------------------------------------//
//              end of file IsotropicSource.cc
//----------------------------------------------------------------------------//
