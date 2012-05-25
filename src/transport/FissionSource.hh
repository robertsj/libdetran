//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   FissionSource.hh
 * \author robertsj
 * \date   Apr 4, 2012
 * \brief  FissionSource class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef FISSIONSOURCE_HH_
#define FISSIONSOURCE_HH_

// Detran
#include "Material.hh"
#include "Mesh.hh"
#include "State.hh"

// Utilities
#include "DBC.hh"
#include "Definitions.hh"
#include "SP.hh"

namespace detran
{

//===========================================================================//
/*!
 * \class FissionSource
 * \brief 
 */
//===========================================================================//

class FissionSource: public Object
{

public:

  typedef SP<FissionSource>                 SP_source;
  typedef State::SP_state                   SP_state;
  typedef Mesh::SP_mesh                     SP_mesh;
  typedef Material::SP_material             SP_material;

  /*!
   *  \brief Constructor
   *
   *  \param state      State vector
   *  \param mesh       Cartesian mesh
   *  \param material   Materials
   *
   */
  FissionSource(SP_state state, SP_mesh mesh, SP_material material);

  /// SP Constructor
  static SP<FissionSource>
  Create(SP<State> state, SP<Mesh> mesh, SP<Material> material)
  {
    SP_source p;
    p = new FissionSource(state, mesh, material);
    return p;
  }

  /// Initialize to thermal fission cross-section, normalized.
  void initialize();

  /// Update the fission density.
  void update();

  /*!
   *   \brief Setup the fission source for an outer iteration.
   *
   *   This sets a new scaling factor \f$ k \f$ and precomputes the
   *   quantity \f$ v = fd(k)^{-1} \f$.
   *
   *   \param scale     Scaling factor (typically 1/keff)
   */
  void setup_outer(double scale);

  /*!
   *   \brief Return the fission source in a group.
   *
   *   The group fission source is just that component of the density
   *   released in  a particular group.  Mathematically, this is just
   *
   *   \f[
   *     q_{f,g} = \frac{\chi_g}{4\pi k} \sum_g \nu\Sigma_{f,g} \phi_g \, .
   *   \f]
   *
   *   Note, the scaling factor is actually arbitrary.  For 2-D and 3-D, it
   *   is \f$ 4\pi \f$, possibly with the eigenvalue \f$ k \f$.  The client
   *   sets this in \ref update.
   *
   *   Note also that this returns a moments source, so the client must
   *   apply the moments-to-discrete operator.
   *
   *   @param   g   Group of the source.
   *   @return      Source vector.
   */
  const State::moments_type& source(int g);

  /*!
   *   \brief Return the fission density.
   *
   *   \f[
   *     fd = \sum_g \nu\Sigma_{f,g} \phi_g \, .
   *   \f]
   *
   *   @return      Fission density vector.
   */
  const State::moments_type& density();

  /*!
   *   \brief Set the fission density.
   *   \param   f   User-defined density.
   */
  void set_density(State::moments_type& f)
  {
    d_density = f;
  }

  /// Unimplemented DBC function
  bool is_valid() const
  {
    return true;
  }

private:

  /// State vector
  SP_state d_state;
  /// Mesh
  SP_mesh d_mesh;
  /// Materials
  SP_material d_material;

  /// \f$ q_{fg} = norm \times \chi_g \sum_g \nu \Sigma_{fg} \phi_g \f$ .
  State::moments_type d_source;
  /// \f$ d = \sum_g \nu \Sigma_{fg} \phi_g \f$ .
  State::moments_type d_density;
  /// Scaling factor
  double d_scale;
  /// Number of groups.
  int d_number_groups;

};

} // end namespace detran

#endif /* FISSIONSOURCE_HH_ */

//---------------------------------------------------------------------------//
//              end of FissionSource.hh
//---------------------------------------------------------------------------//
