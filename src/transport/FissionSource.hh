//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  FissionSource.hh
 *  @brief FissionSource class definition
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_FISSIONSOURCE_HH_
#define detran_FISSIONSOURCE_HH_

#include "transport/transport_export.hh"
#include "transport/State.hh"
#include "geometry/Mesh.hh"
#include "material/Material.hh"
#include "utilities/DBC.hh"
#include "utilities/Definitions.hh"
#include "utilities/SP.hh"

namespace detran
{

/**
 *  @class FissionSource
 *  @brief Defines the isotropic source from fission reactions.
 */
class TRANSPORT_EXPORT FissionSource
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<FissionSource>       SP_fissionsource;
  typedef State::SP_state                           SP_state;
  typedef detran_geometry::Mesh::SP_mesh            SP_mesh;
  typedef detran_material::Material::SP_material    SP_material;
  typedef detran_utilities::vec_int                 vec_int;
  typedef detran_utilities::size_t                  size_t;
  typedef State::moments_type                       moments_type;
  typedef State::vec_moments_type                   vec_moments_type;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *
   *  @param state      State vector
   *  @param mesh       Cartesian mesh
   *  @param material   Materials
   *
   */
  FissionSource(SP_state state, SP_mesh mesh, SP_material material);

  /// SP Constructor
  static SP_fissionsource Create(SP_state    state,
                                 SP_mesh     mesh,
                                 SP_material material);

  //--------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //--------------------------------------------------------------------------//

  // METHODS FOR OUTER ITERATIONS (FISSION *FIXED*)

  /// Initialize to sum of cell nu*fission cross-section, normalized.
  void initialize();

  /// Update the fission density.
  void update();

  /**
   *   @brief Setup the fission source for an outer iteration.
   *
   *   This sets a new scaling factor @f$ C @f$ and precomputes the
   *   quantity @f$ v = C \times fd @f$.
   *
   *   @param scale     Scaling factor (typically 1/keff)
   */
  void setup_outer(const double scale = 1.0);

  /**
   *   @brief Set the scaling factor only
   *   @param scale     Scaling factor (typically 1/keff)
   */
  void set_scale(const double scale = 1.0);

  /**
   *   @brief Return the fission source in a group.
   *
   *   The group fission source is just that component of the density
   *   released in  a particular group.  Mathematically, this is just
   *
   *   @f[
   *       q_{f,g} = @frac{\chi_g}{4\pi k} \sum_g \nu\Sigma_{f,g} \phi_g \, .
   *   @f]
   *
   *   Note, the scaling factor is actually arbitrary.  For 2-D and 3-D, it
   *   is @f$ 4\pi @f$, possibly with the eigenvalue @f$ k @f$.  The client
   *   sets this in \ref update.
   *
   *   Note also that this returns a moments source, so the client must
   *   apply the moments-to-discrete operator.
   *
   *   @param   g   Group of the source.
   *   @return      Source vector.
   */
  const moments_type& source(const size_t g) const;

  /**
   *   @brief Return the fission density.
   *
   *   @f[
   *       fd = \sum_g \nu\Sigma_{f,g} \phi_g \, .
   *   @f]
   *
   *   @return      Fission density vector.
   */
  const moments_type& density()const ;
  moments_type& density();

  /**
   *   @brief Set the fission density.
   *   @param   f   User-defined density.
   */
  void set_density(moments_type& f)
  {
    d_density = f;
  }

  // METHODS TO TREAT FISSION LIKE SCATTER

  /**
   *  @brief Build the within group fission source.
   *
   *  This constructs
   *  @f[
   *      q_g = \chi_g \nu \Sigma_{fg} \phi_g \, .
   *  @f]
   *
   *  @param   g        Group for this problem
   *  @param   phi      Const reference to group flux.
   *  @param   source   Mutable reference to moments source.
   *
   */
  void build_within_group_source(const size_t        g,
                                 const moments_type &phi,
                                 moments_type       &source);

  /**
   *  \brief Build the in-fission source.
   *
   *  This constructs
   *  @f[
   *      q_g = \chi_g \sum^G_{g',g\ne g'} \Sigma_{fg'} \phi_{g'} \, .
   *  @f]
   *
   *  This *assumes* the state is up-to-date.
   *
   *  @param   g        Group for this problem
   *  @param   source   Mutable reference to moments source.
   *
   */
  void build_in_fission_source(const size_t  g,
                               moments_type &source);

  /**
   *   @brief Fill a group source vector with a fission source given
   *          a client-defined flux vector
   *
   *   For Krylov methods, we can bring the entire flux-dependent
   *   terms to the left hand side, thus treating the flux implicitly.
   *   This function allows the client to build the total fission
   *   source into a given group based on an arbitrary, client-defined
   *   multigroup flux vector.
   *
   *   @param   g       group of source being constructed
   *   @param   phi     multigroup fluxes
   *   @param   source  moment vector of group source to contribute to
   */
  void build_total_group_source(const size_t                   g,
                                const State::vec_moments_type &phi,
                                State::moments_type           &source);

  /// Get the state
  SP_state state() {return d_state;}

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// State vector
  SP_state d_state;
  /// Mesh
  SP_mesh d_mesh;
  /// Materials
  SP_material d_material;
  /// Material map
  vec_int d_mat_map;
  /// @f$ q_{fg} = norm \times \chi_g \sum_g \nu \Sigma_{fg} \phi_g @f$ .
  vec_moments_type d_source;
  /// @f$ d = \sum_g \nu \Sigma_{fg} \phi_g @f$ .
  moments_type d_density;
  /// Scaling factor
  double d_scale;
  /// Number of groups.
  size_t d_number_groups;
  /// Adjoint flag
  bool d_adjoint;

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

  /// The "from" group
  size_t g_from(const size_t g, const size_t gp) const;
  /// the "to" group
  size_t g_to(const size_t g, const size_t gp) const;

};

TRANSPORT_TEMPLATE_EXPORT(detran_utilities::SP<FissionSource>)

} // end namespace detran

#include "FissionSource.i.hh"

#endif /* detran_FISSIONSOURCE_HH_ */

//----------------------------------------------------------------------------//
//              end of FissionSource.hh
//----------------------------------------------------------------------------//
