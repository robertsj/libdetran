//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Surface.hh
 *  @brief Surface class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_geometry_SURFACE_HH_
#define detran_geometry_SURFACE_HH_

#include "geometry/Ray.hh"
#include "utilities/SP.hh"

namespace detran_geometry
{

/**
 *  @class Surface
 *  @brief Represents arbitrary three-dimensional surface
 */
class Surface
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<Surface>   SP_surface;
  typedef std::vector<SP_surface>         vec_surface;
  typedef std::vector<Point>              vec_point;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR AND DESTRUCTOR
  //--------------------------------------------------------------------------//

  /// Constructor
  Surface()
  {
    /* ... */
  }

  /// Virtual destructor
  virtual ~Surface()
  {
    /* ... */
  }

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /// Determine whether a point is "outside" (positive) or "inside" the surface
  bool sense(const Point &r)
  {
    return f(r) >= 0.0;
  }

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE --- ALL SURFACE TYPES MUST IMPLEMENT THESE
  //--------------------------------------------------------------------------//

  /// Residual function that implicitly defines the surface, i.e. F(r) = 0
  virtual double f(const Point &r) = 0;

  /**
   *  @brief Compute intersection points made with a ray in order
   *  @param    ray     Ray being cast
   *  @param    t_max   Optional maximum ray length to consider
   */
  virtual vec_point intersections(const Ray    &ray,
                                  const double  t_max = -1) = 0;

};


} // namespace detran_geometry

#endif /* detran_geometry_SURFACE_HH_ */

//----------------------------------------------------------------------------//
//              end of file Surface.hh
//----------------------------------------------------------------------------//
