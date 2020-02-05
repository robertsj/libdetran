//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  PinCell.hh
 *  @brief PinCell class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts  
 */
//----------------------------------------------------------------------------//

#ifndef detran_geometry_PinCell_HH_
#define detran_geometry_PinCell_HH_

#include "Mesh2D.hh"

namespace detran_geometry
{

//----------------------------------------------------------------------------//
/**
 *  @class PinCell
 *  @brief Represents a pin cell
 */
//----------------------------------------------------------------------------//
class GEOMETRY_EXPORT PinCell
{

public:

  /**
   *  @brief Division schemes identifiers
   *
   *  The division schemes are illustrated below.  The dividing planes
   *  pass through the pin center.  The regions created are indexed in
   *  order from the pin center outward, and from the first octant counter
   *  clockwise.
   *
   *  @code
   *   __ __     __ __     _____     __ __
   *  |     |   |  |  |   |\   /|   |\ | /|
   *  |     |   |__|__|   | \ / |   |_\|/_|
   *  |     |   |  |  |   | / \ |   | /|\ |
   *  |__ __|   |__|__|   |/___\|   |/_|_\|
   *
   *  @endcode
   *
   */
  enum DIVISION_SCHEMES
  {
    DIVISION_NONE,
    DIVISION_HV,
    DIVISION_DIAG,
    DIVISION_HV_DIAG,
    END_DIVISION_SCHEMES
  };

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<PinCell>   SP_pincell;
  typedef Mesh::SP_mesh                   SP_mesh;
  typedef Mesh::vec_dbl                   vec_dbl;
  typedef Mesh::vec_int                   vec_int;

  //--------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor.
   *
   *  @param    pitch       Pin cell pitch (assumed square)
   *  @param    radii       Vector of fuel pin radii (can be zero length)
   *  @param    mat_map     Region material map (cell-center outward)
   *  @param    fuel_flag   Designates this as a fuel pin
   */
  PinCell(const Point      &pitch,
          const vec_int    &mat_map,
          const vec_dbl    &radii = vec_dbl(0),
          const size_t      division = DIVISION_NONE,
          const Point      &pincenter = Point(0));

  /// SP Constructor
  static SP_pincell
  Create(const Point      &pitch,
         const vec_int    &mat_map,
         const vec_dbl    &radii = vec_dbl(0),
         const size_t      division = DIVISION_NONE,
         const Point      &pincenter = Point(0));

  /// Return the smart pointer to my mesh.
  Mesh::SP_mesh mesh()
  {
    return d_mesh;
  }

  /**
   *  @brief Mesh the pin cell
   *
   *  There are two discretization options. The first takes the
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
   *
   *  @param number_meshes  Number of meshes per direction.
   *  @param flag           Switch between uniform and "best fit"
   */
  void meshify(int number_meshes, bool flag = false);

  const Point& pitch() const {return d_pitch;}
  const Point& pin_center() const {return d_pin_center;}
  const vec_dbl& radii() const {return d_radii;}
  const vec_int& mat_map() const {return d_mat_map;}

  size_t division() const {return d_division;}

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Pin cell pitch.
  Point d_pitch;
  /// Region material map.
  vec_int d_mat_map;
  /// Pin shell radii.
  vec_dbl d_radii;
  /// Division scheme
  size_t d_division;
  /// Pin center
  Point d_pin_center;
  /// Meshed pin cell
  Mesh2D::SP_mesh d_mesh;

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

  /**
   *  @brief Determine in what region a mesh center resides.
   *
   *
   *  @param  i  Horizontal index.
   *  @param  j  Vertical index.
   *  @return      Region index.
   */
  int find_region(int i, int j, double width);

};

GEOMETRY_TEMPLATE_EXPORT(detran_utilities::SP<PinCell>)
GEOMETRY_TEMPLATE_EXPORT(std::vector<detran_utilities::SP<PinCell> >)

} // end namespace detran_geometry

#endif /* detran_geometry_PinCell_HH_ */

//----------------------------------------------------------------------------//
//              end of PinCell.hh
//----------------------------------------------------------------------------//
