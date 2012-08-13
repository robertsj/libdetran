//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   CoarseMesh.hh
 * \brief  CoarseMesh 
 * \author Jeremy Roberts
 * \date   Aug 8, 2012
 */
//---------------------------------------------------------------------------//

#ifndef COARSEMESH_HH_
#define COARSEMESH_HH_

// Detran
#include "Mesh.hh"

// Utilities
#include "DBC.hh"
#include "SP.hh"

// System
#include <iostream>

namespace detran
{

class CoarseMesh
{

public:

  typedef SP<CoarseMesh>    SP_coarsemesh;
  typedef Mesh::SP_mesh     SP_mesh;

  /*!
   *  \class CoarseMesh
   *  \brief Create a coarse mesh for acceleration
   *
   *  This is used for nonlinear acceleration schemes, and,
   *  in principle, could be used to develop multigrid
   *  preconditioners for Krylov solvers.
   *
   *  This is a very limited interface.  A possible second
   *  approach (for the user) would be to assign
   *  a coarse mesh map (independent of the coarse mesh
   *  used for material assignment)
   */
  CoarseMesh(SP_mesh fine_mesh, u_int level);

  /// Get the original fine mesh
  SP_mesh get_fine_mesh() const
  {
    Require(d_fine_mesh);
    return d_fine_mesh;
  }

  /// Get the coarse mesh
  SP_mesh get_coarse_mesh() const
  {
    Require(d_coarse_mesh);
    return d_coarse_mesh;
  }

  /*!
   *  \brief Get the coarse mesh index for a fine mesh
   *  \param  ijk fine mesh index
   *  \param  dim dimension of index
   *  \return     coarse mesh index
   */
  int fine_to_coarse(u_int ijk, u_int dim) const
  {
    //Require(dim < d_fine_mesh->dimension());
    Require(ijk < d_fine_mesh->number_cells(dim));
    return d_fine_to_coarse[dim][ijk];
  }

  /*!
   *  \brief Determine whether a fine mesh edge is on a coarse mesh edge
   *  \param  ijk fine mesh edge index
   *  \param  dim dimension of index
   *  \return     nonnegative index of coarse edge (otherwise -1)
   */
  int coarse_edge_flag(u_int ijk, u_int dim) const
  {
    Require(dim < d_fine_mesh->dimension());
    Require(ijk < d_fine_mesh->number_cells(dim) + 1);
    return d_coarse_edge_flag[dim][ijk];
  }

private:

  /// Fine mesh
  SP_mesh d_fine_mesh;

  /// Coarse mesh
  SP_mesh d_coarse_mesh;

  /// Level
  const u_int d_level;

  /// Fine to coarse map
  vec2_int d_fine_to_coarse;

  /// Coarse mesh edge flag
  vec2_int d_coarse_edge_flag;

};

} // end namespace detran

#endif // COARSEMESH_HH_ 

//---------------------------------------------------------------------------//
//              end of file CoarseMesh.hh
//---------------------------------------------------------------------------//
