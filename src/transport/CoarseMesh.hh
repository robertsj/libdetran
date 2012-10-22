//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   CoarseMesh.hh
 *  @brief  CoarseMesh
 *  @author Jeremy Roberts
 *  @date   Aug 8, 2012
 */
//---------------------------------------------------------------------------//

#ifndef detran_COARSEMESH_HH_
#define detran_COARSEMESH_HH_

#include "geometry/Mesh.hh"
#include "utilities/DBC.hh"
#include "utilities/SP.hh"
#include <iostream>

namespace detran
{

class CoarseMesh
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<CoarseMesh>  SP_coarsemesh;
  typedef detran_geometry::Mesh::SP_mesh    SP_mesh;
  typedef detran_utilities::size_t          size_t;
  typedef detran_utilities::vec_int         vec_int;
  typedef detran_utilities::vec2_int        vec2_int;
  typedef detran_utilities::vec_dbl         vec_dbl;
  typedef detran_utilities::vec2_dbl        vec2_dbl;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @class CoarseMesh
   *  @brief Create a coarse mesh for acceleration
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
  CoarseMesh(SP_mesh fine_mesh, const size_t level);

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

  /**
   *  @brief Get the coarse mesh index for a fine mesh
   *  @param  ijk fine mesh index
   *  @param  dim dimension of index
   *  @return     coarse mesh index
   */
  size_t fine_to_coarse(const size_t ijk, const size_t dim) const
  {
    //Require(dim < d_fine_mesh->dimension());
    Require(ijk < d_fine_mesh->number_cells(dim));
    return d_fine_to_coarse[dim][ijk];
  }

  /**
   *  @brief Determine whether a fine mesh edge is on a coarse mesh edge
   *  @param  ijk fine mesh edge index
   *  @param  dim dimension of index
   *  @return     nonnegative index of coarse edge (otherwise -1)
   */
  size_t coarse_edge_flag(const size_t ijk, const size_t dim) const
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
  const size_t d_level;
  /// Fine to coarse map
  vec2_int d_fine_to_coarse;
  /// Coarse mesh edge flag
  vec2_int d_coarse_edge_flag;

};

} // end namespace detran

#endif // detran_COARSEMESH_HH_

//---------------------------------------------------------------------------//
//              end of file CoarseMesh.hh
//---------------------------------------------------------------------------//
