//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  Macrobody.hh
 *  @brief Macrobody class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#ifndef detran_geometry_MACROBODY_HH_
#define detran_geometry_MACROBODY_HH_

#include "CSG.hh"
#include "QuadraticSurfaceFactory.hh"

namespace detran_geometry
{

/**
 *  @class Macrobody
 *  @brief Base class for macrobody nodes
 *
 *  A "macrobody" is just a simple solid region.  In most CSG applications,
 *  primitives consist of macrobodies rather than the surfaces that make
 *  up the body.  In some cases, using macrobodies may yield more efficient
 *  tracking since the intersections returned are guaranteed to be with
 *  the macrobody and not with the surfaces.
 */
class Macrobody: public CSG_Node
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef Surface::vec_surface  vec_surface;

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /// Constructor with optional origin
  Macrobody(c_Point &origin = Point(0, 0, 0));
  /// Does the node contain the point?
  bool contains(c_Point &r) const;
  /// Where does the node intersect the ray?
  virtual vec_point intersections(const Ray &r, c_dbl t_max) = 0;

protected:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Surfaces that define the macrobody
  vec_surface d_surfaces;
  /// Origin
  Point d_origin;

};

/**
 *  @class RightParallelpiped
 *  @brief Orthogonal cuboid (or right parallelpiped)
 *
 *  Equivalent to MCNP's RPP
 */
class RightParallelpiped: public Macrobody
{

public:

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *
   *  To obtain a 2-D parallelpiped, set the depth to zero.
   *
   *  @param  W       width (x)
   *  @param  L       length (y)
   *  @param  D       depth (z)
   *  @param  origin  bottom left corner of the cuboid
   */
  RightParallelpiped(c_dbl W, c_dbl L, c_dbl D, c_Point &origin);

  /// Does the node contain the point?
  bool contains(c_Point &r) const;
  /// Where does the node intersect the ray?
  vec_point intersections(const Ray &r, c_dbl t_max);

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  //@{
  ///  Dimensions
  double d_length;
  double d_width;
  double d_depth;
  //@}

};


/**
 *  @class RightCircularCylinder
 *  @brief RightCircularCylinder aligned with z-axis
 *
 *  Equivalent to MCNP's RCC
 */
class RightCircularCylinder: public Macrobody
{

public:

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *
   *  To obtain a 2-D cylinder, set the height to zero.
   *
   *  @param  R       radius
   *  @param  H       height
   *  @param  origin  bottom center of the cylinder
   */
  RightCircularCylinder(c_dbl    R,
                        c_dbl    H,
                        c_Point &origin = Point(0, 0, 0));

  /// Does the node contain the point?
  bool contains(c_Point &r) const;
  /// Where does the node intersect the ray?
  vec_point intersections(const Ray &r, c_dbl t_max);

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  //@{
  /// Dimensions
  double d_radius;
  double d_height;
  //@}

};

/**
 *  @class RightHexagonalPrism
 *  @brief HexagonalPrism aligned with z-axis
 *
 *  Equivalent to MCNP's RHP or HEX
 */
class RightHexagonalPrism: public Macrobody
{

public:

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *
   *  To obtain a 2-D prism, set the height to zero.
   *
   *  @param  L       length of a side
   *  @param  H       height
   *  @param  rotate  western point (false) or face (true)
   *  @param  origin  bottom center of the prism
   */
  RightHexagonalPrism(c_dbl    L,
                      c_dbl    H,
                      bool     rotate = false,
                      c_Point &origin = Point(0, 0, 0));

  /// Does the node contain the point?
  bool contains(c_Point &r) const;
  /// Where does the node intersect the ray?
  vec_point intersections(const Ray &r, c_dbl t_max);

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  //@{
  /// Dimensions
  double d_length;
  double d_height;
  //@}

};

} // end namespace detran_geometry

#endif /* detran_geometry_MACROBODY_HH_ */

//----------------------------------------------------------------------------//
//              end of Macrobody.hh
//----------------------------------------------------------------------------//
