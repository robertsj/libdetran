//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  ScatterSource.hh
 *  @brief ScatterSource class definition
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_SCATTERSOURCE_HH_
#define detran_SCATTERSOURCE_HH_

#include "transport/transport_export.hh"
#include "transport/State.hh"
#include "material/Material.hh"
#include "geometry/Mesh.hh"
#include "utilities/DBC.hh"
#include "utilities/MathUtilities.hh"
#include "utilities/SP.hh"
#include <iostream>

namespace detran
{

//----------------------------------------------------------------------------//
/**
 *  @class ScatterSource
 *  @brief Methods for constructing various scattering sources.
 *
 *  See the individual methods for detailed information.
 */
//----------------------------------------------------------------------------//

class TRANSPORT_EXPORT ScatterSource
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<ScatterSource>       SP_scattersource;
  typedef detran_geometry::Mesh::SP_mesh            SP_mesh;
  typedef detran_material::Material::SP_material    SP_material;
  typedef State::SP_state                           SP_state;
  typedef State::moments_type                       moments_type;
  typedef detran_utilities::vec_int                 vec_int;
  typedef detran_utilities::vec_size_t              vec_size_t;
  typedef detran_utilities::size_t                  size_t;
  typedef vec_size_t                                groups_t;
  typedef groups_t::iterator                        groups_iter;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param mesh           Pointer to mesh
   *  @param material       Pointer to material
   *  @param state          Pointer to state vector
   */
  ScatterSource(SP_mesh mesh, SP_material material, SP_state state);

  //--------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //--------------------------------------------------------------------------//

  /**
   *  @brief Build the within group scattering source.
   *
   *  This constructs
   *  @f[
   *      q_g = \mathbf{S}_{gg} \phi_g \, .
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
   *  @brief Build the in-scatter source.
   *
   *  This constructs
   *  @f[
   *      q_g = \sum^G_{g',g\ne g'} \mathbf{S}_{gg'}\phi_{g'}
   *  @f]
   *  or for adjoint problems
   *  @f[
   *      q_g = \sum^G_{g',g\ne g'} \mathbf{S}_{g'g}\phi_{g'} \, .
   *  @f]
   *
   *  This \e assumes the state is up-to-date.
   *
   *  @param   g    Group for this problem
   *  @param   s    Mutable reference to moments source.
   *
   */
  void build_in_scatter_source(const size_t  g,
                               moments_type &s);
  /**
   *  @brief Build the downscatter source.
   *
   *  This constructs
   *  @f[
   *      q_g = \sum^{g_{cutoff}_{g'} \mathbf{S}_{gg'}\phi_{g'} \, .
   *  @f]
   *
   *  This \e assumes the state is up-to-date.
   *
   *  This is useful when creating the fixed source for Krylov
   *  multigroup solves when Gauss-Seidel has been used for the
   *  downscatter block.
   *
   *  @param   g        Group for this problem
   *  @param   g_cutoff Highest group to contribute to downscatter
   *  @param   s        Mutable reference to moments source.
   */
  void build_downscatter_source(const size_t  g,
                                const size_t  g_cutoff,
                                moments_type &s);

  /**
   *  @brief Build the total scatter source.
   *
   *  In some cases, including all scattering is required, as is the case
   *  when performing multigroup Krylov solves.  This constructs
   *
   *  @f[
   *      q_g = \sum^G_{g'} \mathbf{S}_{gg'}\phi_{g'} \, .
   *  @f]
   *
   *  Because Gauss-Seidel can be used to solve downscatter blocks,
   *  a cutoff group is passed to exclude the solved portion of
   *  the problem.
   *
   *  @param   g        Group for this problem.
   *  @param   g_cutoff Highest group to contribute to downscatter.
   *  @param   phi      Const reference to multigroup flux moments.
   *  @param   s        Mutable reference to moments source.
   *
   */
  void build_total_group_source(const size_t                   g,
                                const size_t                   g_cutoff,
                                const State::vec_moments_type &phi,
                                moments_type                  &s);

protected:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Mesh
  SP_mesh d_mesh;
  /// Material
  SP_material d_material;
  /// State
  SP_state d_state;
  /// Material map
  vec_int d_mat_map;
  /// Adjoint
  bool d_adjoint;
  /// Vector of group indices for in-scatter, down-scatter, and total-scatter
  groups_t d_groups;

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

  /// Lower bound for a group loop
  size_t lower(const size_t g) const;
  /// Upper bound for a group loop
  size_t upper(const size_t g) const;
  /// The "from" group
  size_t g_from(const size_t g, const size_t gp) const;
  /// the "to" group
  size_t g_to(const size_t g, const size_t gp) const;
  /// Set group index vector
  groups_iter groups(const size_t g_start, const size_t g_finish, bool inc);

};

TRANSPORT_TEMPLATE_EXPORT(detran_utilities::SP<ScatterSource>)

} // end namespace detran

//----------------------------------------------------------------------------//
// INLINE MEMBER DEFINITIONS
//----------------------------------------------------------------------------//

#include "ScatterSource.i.hh"

#endif /* detran_SCATTERSOURCE_HH_ */

//----------------------------------------------------------------------------//
//              end of ScatterSource.hh
//----------------------------------------------------------------------------//
