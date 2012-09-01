//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ScatterSource.hh
 * \author robertsj
 * \date   Apr 4, 2012
 * \brief  ScatterSource class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef SCATTERSOURCE_HH_
#define SCATTERSOURCE_HH_

// Detran
#include "Material.hh"
#include "Mesh.hh"
#include "State.hh"

// Utilities
#include "DBC.hh"
#include "SP.hh"

#include <iostream>

namespace detran
{

//---------------------------------------------------------------------------//
/*!
 *  \class ScatterSource
 *  \brief Methods for constructing various scattering sources.
 */
//---------------------------------------------------------------------------//

class ScatterSource: public Object
{

public:

  typedef SP<ScatterSource>         SP_source;
  //
  typedef Mesh::SP_mesh             SP_mesh;
  typedef Material::SP_material     SP_material;
  typedef State::SP_state           SP_state;
  //
  typedef State::moments_type       moments_type;

  ScatterSource(SP_mesh mesh, SP_material material, SP_state state);

  /*!
   *  \brief Build the within group scattering source.
   *
   *  This constructs
   *  \f[
   *      q_g = \mathbf{S}_{gg} \phi_g \, .
   *  \f]
   *
   *  \param   g        Group for this problem
   *  \param   phi      Const reference to group flux.
   *  \param   source   Mutable reference to moments source.
   *
   */
  void build_within_group_source(int g,
                                 const moments_type &phi,
                                 moments_type &source)
  {
    // Preconditions
    Require( (g >= 0) and (g < d_material->number_groups()) );
    Require( phi.size() == source.size() );

    for (int cell = 0; cell < d_mesh->number_cells(); cell++)
    {
      source[cell] += phi[cell] * d_material->sigma_s(d_mat_map[cell], g, g);
    }
  }

  /*!
   *  \brief Build the in-scatter source.
   *
   *  This constructs
   *  \f[
   *      q_g = \sum^G_{g',g\ne g'} \mathbf{S}_{gg'}\phi_{g'} \, .
   *  \f]
   *
   *  This *assumes* the state is up-to-date.
   *
   *  \param   g        Group for this problem
   *  \param   source   Mutable reference to moments source.
   *  \param   g_down   Highest group to contribute to downscatter
   *  \param   g_up     Lowest group to contribute to upscatter
   *
   */
  void build_in_scatter_source(int g,
                               moments_type &source)
  {
    // Preconditions.
    Require( (g >= 0) and (g < d_material->number_groups()) );

    // Add downscatter.
    for (int gp = d_material->lower(g); gp < g; gp++) //
    {
      const moments_type &phi = d_state->phi(gp);
      for (int cell = 0; cell < d_mesh->number_cells(); cell++)
      {
        source[cell] += phi[cell] * d_material->sigma_s(d_mat_map[cell], g, gp);
      }
    }
    // Add upscatter.
    for (int gp = g + 1; gp <= d_material->upper(g); gp++)
    {
      const moments_type &phi = d_state->phi(gp);
      for (int cell = 0; cell < d_mesh->number_cells(); cell++)
      {
        source[cell] += phi[cell] * d_material->sigma_s(d_mat_map[cell], g, gp);
      }
    }
  }

  /*!
   *  \brief Build the downscatter source.
   *
   *  This constructs
   *  \f[
   *      q_g = \sum^{g_{cutoff}_{g'} \mathbf{S}_{gg'}\phi_{g'} \, .
   *  \f]
   *
   *  This \e assumes the state is up-to-date.
   *
   *  This is useful when creating the fixed source for Krylov
   *  multigroup solves when Gauss-Seidel has been used for the
   *  downscatter block.
   *
   *  \param   g        Group for this problem
   *  \param   g_down   Highest group to contribute to downscatter
   *  \param   source   Mutable reference to moments source.
   */
  void build_downscatter_source(int g,
                                int g_cutoff,
                                moments_type &source)
  {
    // Preconditions.
    Require( (g >= 0) and (g < d_material->number_groups()) );
    Require( (g_cutoff >= 0) and (g_cutoff <= d_material->number_groups()) );

    // Add downscatter.
    for (int gp = d_material->lower(g); gp < g_cutoff; gp++) //
    {
      moments_type phi = d_state->phi(gp);
      for (int cell = 0; cell < d_mesh->number_cells(); cell++)
      {
        source[cell] += phi[cell] * d_material->sigma_s(d_mat_map[cell], g, gp);
      }
    }

  }

  /*!
   *  \brief Build the total scatter source.
   *
   *  In some cases, including all scattering is required, as is the case
   *  when performing multigroup Krylov solves.  This constructs
   *
   *  \f[
   *      q_g = \sum^G_{g'} \mathbf{S}_{gg'}\phi_{g'} \, .
   *  \f]
   *
   *  Because Gauss-Seidell can be used to solve downscatter blocks,
   *  a cutoff group is passed to exclude the solved portion of
   *  the problem.
   *
   *  \param   g        Group for this problem.
   *  \param   g_cutoff Highest group to contribute to downscatter.
   *  \param   phi      Const reference to multigroup flux moments.
   *  \param   source   Mutable reference to moments source.
   *
   */
  void build_total_group_source(int g,
                                int g_cutoff,
                                const State::vec_moments_type &phi,
                                moments_type &source)
  {
    Require( (g >= 0) and (g < d_material->number_groups()) );


    for (int gp = g_cutoff; gp <= d_material->upper(g); gp++)
    {
      for (int cell = 0; cell < d_mesh->number_cells(); cell++)
      {
        source[cell] += phi[gp][cell] *
            d_material->sigma_s(d_mat_map[cell], g, gp);
      }
    }
  }

  bool is_valid() const
  {
    return true;
  }

protected:

  /// Mesh
  SP_mesh d_mesh;

  /// Material
  SP_material d_material;

  /// State
  SP_state d_state;

  /// Material map
  vec_int d_mat_map;

};

} // namespace detran

#endif /* SCATTERSOURCE_HH_ */

//---------------------------------------------------------------------------//
//              end of ScatterSource.hh
//---------------------------------------------------------------------------//
