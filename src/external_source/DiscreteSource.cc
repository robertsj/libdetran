//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  DiscreteSource.cc
 *  @brief DiscreteSource member definitions
 *  @note  Copyright (C) Jeremy Roberts 2012-2013
 */
//----------------------------------------------------------------------------//

#include "DiscreteSource.hh"

namespace detran_external_source
{

//----------------------------------------------------------------------------//
DiscreteSource::DiscreteSource(size_t         number_groups,
                               SP_mesh        mesh,
                               vec3_dbl       spectra,
                               vec_int        map,
                               SP_quadrature  quadrature)
  : ExternalSource(number_groups, mesh, quadrature, true)
  , d_source_spectra(spectra)
  , d_source_map(map)
{
  Require(d_source_spectra.size()       > 0);
  Require(d_source_spectra[0].size()    == d_number_groups);
  Require(d_source_spectra[0][0].size() == d_number_angles);
  Require(d_source_map.size()           == d_mesh->number_cells());

  d_mesh->add_mesh_map("DISCRETESOURCE", d_source_map);
}

} // end namespace detran_external_source

//----------------------------------------------------------------------------//
//              end of file DiscreteSource.cc
//----------------------------------------------------------------------------//
