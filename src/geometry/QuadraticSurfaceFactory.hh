//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  QuadraticSurfaceFactory.hh
 *  @brief QuadraticSurfaceFactory struct definition
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
  typedef const double          c_dbl;

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  static SP_surface CreatePlaneX(c_dbl x)
  {
    return QS::Create(0, 0, 0,   0, 0, 0,   1, 0, 0,  -x);
  }

  static SP_surface CreatePlaneY(c_dbl y)
  {
    return QS::Create(0, 0, 0,   0, 0, 0,   0, 1, 0,  -y);
  }

  static SP_surface CreatePlaneZ(c_dbl z)
  {
    return QS::Create(0, 0, 0,   0, 0, 0,   0, 0, 1,  -z);
  }

  /// a*x + b*y + c*z = d
  static SP_surface CreatePlane(c_dbl a, c_dbl b, c_dbl c, c_dbl d)
  {
    return QS::Create(0, 0, 0,   0, 0, 0,   a, b, c,  -d);
  }

  /// (y-y0)^2 + (z-z0)^2 = r^2
  static SP_surface CreateCylinderX(c_dbl y0, c_dbl z0, c_dbl r)
  {
    return QS::Create(0, 1, 1,   0, 0, 0,   0, -2*y0, -2*z0,   y0*y0+z0*z0-r*r);
  }

  /// (x-x0)^2 + (z-z0)^2 = r^2
  static SP_surface CreateCylinderY(c_dbl x0, c_dbl z0, c_dbl r)
  {
    return QS::Create(1, 0, 1,   0, 0, 0,   -2*x0, 0, -2*z0,   x0*x0+z0*z0-r*r);
  }

  /// (x-x0)^2 + (y-y0)^2 = r^2
  static SP_surface CreateCylinderZ(c_dbl x0, c_dbl y0, c_dbl r)
  {
    return QS::Create(1, 1, 0,  0, 0, 0,  -2*x0, -2*y0, 0,  x0*x0+y0*y0-r*r);
  }

  /// (x-x0)^2 + (y-y0)^2 + (z-z0)^2 = r^2
  static SP_surface CreateSphere(c_dbl x0, c_dbl y0, c_dbl z0, c_dbl r)
  {
    return QS::Create(1,1,1, 0,0,0, -2*x0,-2*y0,-2*z0, x0*x0+y0*y0+z0*z0-r*r);
  }

};

} // end namespace detran_geometry

#endif /* QUADRATICSURFACEFACTORY_HH_ */

//----------------------------------------------------------------------------//
//              end of file QuadraticSurfaceFactory.hh
//----------------------------------------------------------------------------//
