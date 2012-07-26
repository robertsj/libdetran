//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   FixedSource.hh
 * \author robertsj
 * \date   Jul 26, 2012
 * \brief  FixedSource class definition.
 */
//---------------------------------------------------------------------------//

#ifndef FIXEDSOURCE_HH_
#define FIXEDSOURCE_HH_

// Detran
#include "Material.hh"
#include "Mesh.hh"

// Utilities
#include "Definitions.hh"
#include "InputDB.hh"

// System
#include "petsc.h"

namespace detran_diffusion
{

/*!
 *  \class FixedSource
 *  \brief Fixed source for diffusion problems.
 *
 *  Fixed sources can be be volumetric or on the boundary
 *  as incident partial currents.
 *
 */
class FixedSource
{

public:

  typedef detran::InputDB::SP_input     SP_input;
  typedef detran::Material::SP_material SP_material;
  typedef detran::Mesh::SP_mesh         SP_mesh;
  typedef detran::vec_dbl               vec_dbl;
  typedef detran::vec2_dbl              vec2_dbl;

  /*!
   *  \brief Constructor
   *  \param input          Input database
   *  \param mesh           Mesh
   *  \param number_groups  Number of energy groups

   */
  FixedSource(SP_input      input
              SP_mesh       mesh
              int           number_groups);

  /*!
   *  \brief Add a volume source
   *
   *  \param spectra        Energy spectra
   *  \param spectra_map    Map of spectra locations
   */
  add_volume_source(vec2_dbl      &spectra,
                    vec_int       &spectra_map)

  /*!
   *  \brief Add a boundary source.
   *
   *  Note, surface sources are defined in terms of
   *  partial currents.
   *
   *  \param number_groups  Number of energy groups
   *  \param spectra        Energy spectra
   *  \param spectra_map    Map of spectra locations
   */
  add_boundary_source(vec2_dbl      &source,
                      int           side);

private:

  /// \name Private Data
  /// \{

  /// Right hand side vector
  Vec d_b;

  /// Energy groups
  int d_number_groups;

  /// One group spatial size
  int d_group_size;

  /// \}

  /// \name Implementation
  /// \{

  void construct();

  /// \}

};

} // end namespace detran_diffusion


#endif /* FIXEDSOURCE_HH_ */
