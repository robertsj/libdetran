//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  QuadraticSurfaceFactory.hh
 *  @brief QuadraticSurfaceFactory
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//


#ifndef detran_geometry_QUADRATICSURFACEFACTORY_HH_
#define detran_geometry_QUADRATICSURFACEFACTORY_HH_

#include "QuadraticSurface.hh"

namespace detran_geometry
{

/**
 *  @struct QuadraticSurfaceFactory
 *  @brief  Convenience functions for creating simple quadratic surfaces
 */
struct QuadraticSurfaceFactory
{

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef QuadraticSurface      QS;
  typedef QS::SP_surface        SP_surface;
  typedef const double          c_double;

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  static SP_surface PlaneX(c_double x)
  {
    return QS::Create(0, 0, 0,   0, 0, 0,   1, 0, 0,  -x);
  }

  static SP_surface PlaneY(c_double y)
  {
    return QS::Create(0, 0, 0,   0, 0, 0,   0, 1, 0,  -y);
  }

  static SP_surface PlaneZ(c_double z)
  {
    return QS::Create(0, 0, 0,   0, 0, 0,   0, 0, 1,  -z);
  }

  /// a*x + b*y + c*z = d
  static SP_surface Plane(c_double a, c_double b, c_double c, c_double d)
  {
    return QS::Create(0, 0, 0,   0, 0, 0,   a, b, c,  -d);
  }

  /// (y-y0)^2 + (z-z0)^2 = r^2
  static SP_surface CylinderX(c_double y0, c_double z0, c_double r)
  {
    return QS::Create(0, 1, 1,   0, 0, 0,   0, -2*y0, -2*z0,   y0*y0+z0*z0-r*r);
  }

  /// (x-x0)^2 + (z-z0)^2 = r^2
  static SP_surface CylinderY(c_double x0, c_double z0, c_double r)
  {
    return QS::Create(1, 0, 1,   0, 0, 0,   -2*x0, 0, -2*z0,   x0*x0+z0*z0-r*r);
  }

  /// (x-x0)^2 + (y-y0)^2 = r^2
  static SP_surface CylinderZ(c_double x0, c_double y0, c_double r)
  {
    return QS::Create(1, 1, 0,   0, 0, 0,   -2*x0, -2*y0, 0,   x0*x0+y0*y0-r*r);
  }

};

} // end namespace detran_geometry

#endif /* QUADRATICSURFACEFACTORY_HH_ */

//----------------------------------------------------------------------------//
//              end of file QuadraticSurfaceFactory.hh
//----------------------------------------------------------------------------//
