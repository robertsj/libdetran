//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  DiscreteSource.hh
 *  @brief DiscreteSource class definition
 *  @note  Copyright (C) Jeremy Roberts 2012-2013
 */
//----------------------------------------------------------------------------//

#ifndef detran_external_source_DISCRETESOURCE_HH_
#define detran_external_source_DISCRETESOURCE_HH_

#include "ExternalSource.hh"

namespace detran_external_source
{

//----------------------------------------------------------------------------//
/**
 *  @class DiscreteSource
 *  @brief Discrete source definition
 */
//----------------------------------------------------------------------------//

class EXTERNAL_SOURCE_EXPORT DiscreteSource: public ExternalSource
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::vec3_dbl            vec3_dbl;
  typedef detran_utilities::vec_int             vec_int;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param number_groups  Number of energy groups
   *  @param mesh           Pointer to mesh
   *  @param spectra        Angle-dependent spectra,
   *                          size = [#spectra][#groups][#angle]
   *  @param map            Map of where spectra are located, size = [#cells]
   *  @param quadrature     Quadrature
   */
  DiscreteSource(size_t         number_groups,
                 SP_mesh        mesh,
                 vec3_dbl       spectra,
                 vec_int        map,
                 SP_quadrature  quadrature);

  static SP_externalsource
  Create(size_t                     number_groups,
         SP_mesh                    mesh,
         detran_utilities::vec3_dbl spectra,
         detran_utilities::vec_int  map,
         SP_quadrature  quadrature)
  {
    SP_externalsource
      p(new DiscreteSource(number_groups, mesh, spectra, map, quadrature));
    return p;
  }

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL EXTERNAL SOURCES MUST IMPLEMENT THESE
  //--------------------------------------------------------------------------//

  double source(const size_t cell, const size_t group)
  {
    Require(cell < d_mesh->number_cells());
    Require(group < d_number_groups);
    double value = 0;
    for (size_t o = 0; o < d_quadrature->number_octants(); ++o)
    {
      for (size_t a = 0; a < d_quadrature->number_angles_octant(); ++a)
      {
        int angle = d_quadrature->index(o, a);
        value += d_source_spectra[d_source_map[cell]][group][angle] *
                 d_quadrature->weight(a);
      }
    }
    return value;
  }

  double source(const size_t cell,
                const size_t group,
                const size_t angle)
  {
    Require(cell < d_mesh->number_cells());
    Require(group < d_number_groups);
    Require(angle < d_number_angles);
    return d_source_spectra[d_source_map[cell]][group][angle];
  }

private:


  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Angle-dependent source spectra
  vec3_dbl d_source_spectra;
  /// Fine mesh source map
  vec_int d_source_map;

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

};

} // end namespace detran_external_source

#endif /* detran_external_source_DISCRETESOURCE_HH_ */

//----------------------------------------------------------------------------//
//              end of DiscreteSource.hh
//----------------------------------------------------------------------------//
