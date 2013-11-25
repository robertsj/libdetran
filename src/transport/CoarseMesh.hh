//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  CoarseMesh.hh
 *  @brief CoarseMesh class definition
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 *
 *  \todo Rename to CoarseMesher
 */
//----------------------------------------------------------------------------//

#ifndef detran_COARSEMESH_HH_
#define detran_COARSEMESH_HH_

#include "transport/transport_export.hh"
#include "geometry/Mesh.hh"
#include "utilities/DBC.hh"
#include "utilities/SP.hh"
#include <iostream>

namespace detran
{

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
/**
 *  @example transport/test/test_CoarseMesh.cc
 *  @brief   Test of CoarseMesh class
 */
class TRANSPORT_EXPORT CoarseMesh
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<CoarseMesh>  SP_coarsemesh;
  typedef detran_geometry::Mesh::SP_mesh    SP_mesh;
  typedef detran_utilities::size_t          size_t;
  typedef detran_utilities::vec_int         vec_int;
  typedef detran_utilities::vec2_int        vec2_int;
  typedef detran_utilities::vec_dbl         vec_dbl;
  typedef detran_utilities::vec2_dbl        vec2_dbl;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param fine_mesh    The original fine mesh
   *  @param level        Number of fine meshes desired per coarse mesh
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
   *  @brief Get the coarse mesh index for a fine mesh in one dimension
   *
   *  This functionality is useful when sweeping through the fine mesh. Any
   *  cell (i, j, k) can have up to three surfaces (in x, y, and z) coincident
   *  with a coarse cell (I, J, K).  For tallying surface integrals, the
   *  surface area contribution of (i, j, k) is needed, and this method
   *  allows the user to query the fine mesh cell and hence dimensions in a
   *  general fashion.
   *
   *  @param  ijk fine mesh index
   *  @param  dim dimension of index
   *  @return     coarse mesh index along the dimension of interest
   */
  size_t fine_to_coarse(const size_t ijk, const size_t dim) const
  {
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

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

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

//----------------------------------------------------------------------------//
//              end of file CoarseMesh.hh
//----------------------------------------------------------------------------//
