//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  IsotropicSource.hh
 *  @brief IsotropicSource class definition
 *  @note  Copyright (C) Jeremy Roberts 2012-2013
 */
//----------------------------------------------------------------------------//

#ifndef detran_external_source_ISOTROPICSOURCE_HH_
#define detran_external_source_ISOTROPICSOURCE_HH_

#include "ExternalSource.hh"

namespace detran_external_source
{

//----------------------------------------------------------------------------//
/**
 *  @class IsotropicSource
 *  @brief Isotropic volume source.
 */
//----------------------------------------------------------------------------//

class EXTERNAL_SOURCE_EXPORT IsotropicSource : public ExternalSource
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::vec2_dbl             spectra_type;
  typedef detran_utilities::vec_int              vec_int;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param number_groups  Number of energy groups
   *  @param mesh           Pointer to mesh
   *  @param spectra        Vector of spectra, size = [#spectra][#groups].  The
   *                        unit is n/cc-sec.
   *  @param map            Map of where spectra are located, size = [#cells]
   *  @param quadrature     Pointer to quadrature (optional)
   */
  IsotropicSource(size_t number_groups,
                  SP_mesh mesh,
                  spectra_type &spectra,
                  vec_int &map,
                  SP_quadrature quadrature = SP_quadrature(0));

  /// SP constructor
  static SP_externalsource
  Create(detran_utilities::size_t          number_groups,
         detran_geometry::Mesh::SP_mesh    mesh,
         std::vector<std::vector<double> > &spectra,
         std::vector<int>                  &map,
         SP_quadrature quadrature = SP_quadrature(0))
  {
    SP_externalsource
      p(new IsotropicSource(number_groups, mesh, spectra, map, quadrature));
    return p;
  }

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL EXTERNAL SOURCES MUST IMPLEMENT THESE
  //--------------------------------------------------------------------------//

  double source(const size_t cell, const size_t group)
  {
    // Preconditions
    Require(cell < d_mesh->number_cells());
    Require(group < d_number_groups);

    return d_source_spectra[d_source_map[cell]][group];
  }

  double source(const size_t cell,
                const size_t group,
                const size_t angle)
  {
    Require(cell < d_mesh->number_cells());
    Require(group < d_number_groups);
    Require(angle < d_number_angles); // note, #angle = -1 if no quad set

    return d_source_spectra[d_source_map[cell]][group] * d_norm;
  }

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Source spectra
  spectra_type d_source_spectra;
  /// Fine mesh source map
  vec_int d_source_map;
  /// Angular norm
  double d_norm;

};

} // end namespace detran

#endif /* detran_external_source_ISOTROPICSOURCE_HH_ */

//----------------------------------------------------------------------------//
//              end of IsotropicSource.hh
//----------------------------------------------------------------------------//
