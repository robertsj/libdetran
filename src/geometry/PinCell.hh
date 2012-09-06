//----------------------------------*-C++-*----------------------------------//
/*!
 *  \file   PinCell.hh
 *  \author Jeremy Roberts
 *  \brief  PinCell class definition
 *  \date   Mar 23, 2012
 */
//---------------------------------------------------------------------------//

#ifndef PinCell_HH_
#define PinCell_HH_

#include "Mesh2D.hh"

namespace detran_geometry
{

//---------------------------------------------------------------------------//
/*!
 *  \class PinCell
 *  \brief Two-dimensional Cartesian mesh.
 *
 *  This is mostly a convenience interface for producing a 2D pin cell
 *  mesh.  Two discretization options are provided.
 *
 *  The first takes the
 *  actual geometry (i.e. a square with possible concentric cylinders
 *  at its center) and defines a uniformly-spaced mesh.  The materials
 *  assigned for each mesh are those at the mesh cell center.  This, of
 *  course, is not a nonconservative approximation, but it is
 *  simple and lends itself to cases required a uniform mesh.
 *
 *  The second approximates the cylinder using a single step staircase
 *  approximation.  The cylinder is a union of a square region of
 *  side length \f$ L \f$.  Each face has a rectangular extension
 *  centered at the face center, of width \f$ L/2 \f$, and thickness
 *  $\f \delta \f$.  A brute-force integration couple with nonlinear
 *  optimization was used to find the optimal relation between
 *  the pin radius \f$ r \f$, the box width \f$ L \f$, and the thickness
 *  \f$ \delta \f$ such that volume is conserved and the pin area
 *  not covered by the approximate pin is minimized.
 *  These were
 *  \f [
 *       L = 1.60496875 r \, ,
 *  \f ]
 *  and
 *  \f [
 *       \delta = 0.176223981031790 r \, .
 *  \f ]
 *  Note, the initial decision to use a single step of half the face
 *  width was arbitrary.
 */
//---------------------------------------------------------------------------//
class PinCell
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<PinCell>   SP_pincell;
  typedef Mesh::SP_mesh                   SP_mesh;
  typedef Mesh::vec_dbl                   vec_dbl;
  typedef Mesh::vec_int                   vec_int;

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Constructor.
   *
   *  \param    pitch       Pin cell pitch (assumed square)
   *  \param    radii       Vector of fuel pin radii (can be zero length)
   *  \param    mat_map     Region material map (cell-center outward)
   *  \param    fuel_flag   Designates this as a fuel pin
   */
  PinCell(double pitch, vec_dbl radii, vec_int mat_map, bool fuel_flag = true);

  /// SP Constructor
  static SP_pincell
  Create(double pitch, vec_dbl radii, vec_int mat_map)
  {
    SP_pincell p(new PinCell(pitch, radii, mat_map));
    return p;
  }

  /// Return the smart pointer to my mesh.
  Mesh::SP_mesh mesh()
  {
    return d_mesh;
  }

  /*!
   *  \brief Mesh the pin cell.
   *  \param number_meshes  Number of meshes per direction.
   *  \param flag           Switch between uniform and "best fit"
   */
  void meshify(int number_meshes, bool flag = false);

  bool is_fuel() const
  {
    return d_fuel_flag;
  }

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Underlying meshed object
  Mesh2D::SP_mesh d_mesh;

  /// Pin cell pitch.
  double d_pitch;

  /// Pin shell radii.
  vec_dbl d_radii;

  // Region material map.
  vec_int d_mat_map;

  /// Fuel flag (true if fuel).
  bool d_fuel_flag;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  /*!
   * \brief Determine in what region a mesh center resides.
   *
   *
   * \param  i  Horizontal index.
   * \param  j  Vertical index.
   * \return      Region index.
   */
  int find_region(int i, int j, double width);

};

} // end namespace detran_geometry

#endif /* PinCell_HH_ */

//---------------------------------------------------------------------------//
//              end of PinCell.hh
//---------------------------------------------------------------------------//
