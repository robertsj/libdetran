//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   IsotropicSource.hh
 * \author robertsj
 * \date   Apr 4, 2012
 * \brief  IsotropicSource class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef ISOTROPICSOURCE_HH_
#define ISOTROPICSOURCE_HH_

#include "ExternalSource.hh"

namespace detran
{

//===========================================================================//
/*!
 * \class IsotropicSource
 * \brief Isotropic volume source.
 */
//===========================================================================//

class IsotropicSource : public ExternalSource
{

public:

  // Source types
  typedef SP<IsotropicSource>  SP_source;
  typedef vec2_dbl             spectra_type;

  IsotropicSource(SP_mesh mesh,
                  SP_quadrature quadrature,
                  int number_groups,
                  spectra_type &spectra,
                  vec_int &map)
    :  ExternalSource(mesh, quadrature, number_groups)
  {
    set_source(spectra, map);
  }

  virtual double source(int cell, int group)
  {
    Require(cell >= 0);
    Require(cell < d_mesh->number_cells());
    Require(group >= 0);
    Require(group < d_number_groups);
    return d_source_spectra[d_source_map[cell]][group];
  }

  virtual double source(int cell, int group, int angle)
  {
    Require(cell >= 0);
    Require(cell < d_mesh->number_cells());
    Require(group >= 0);
    Require(group < d_number_groups);
    Require(angle >=0);
    Require(angle < d_number_angles);
    return d_source_spectra[d_source_map[cell]][group] * 1.0; // \todo norm
  }

  void set_source(spectra_type &spectra, vec_int &map)
  {
    Require(spectra.size()>0);
    Require(spectra[0].size() == d_number_groups);
    d_source_spectra = spectra;
    d_mesh->add_coarse_mesh_map("ISOTROPICSOURCE", map);
    d_source_map = d_mesh->mesh_map("ISOTROPICSOURCE");
  }

private:

  /// Source spectra
  spectra_type d_source_spectra;

  /// Fine mesh source map
  vec_int d_source_map;

};

} // end namespace detran

#endif /* ISOTROPICSOURCE_HH_ */

//---------------------------------------------------------------------------//
//              end of IsotropicSource.hh
//---------------------------------------------------------------------------//
