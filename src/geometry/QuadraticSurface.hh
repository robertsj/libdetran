//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  QuadraticSurface.hh
 *  @brief QuadraticSurface class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_geometry_QUADRATICSURFACE_HH_
#define detran_geometry_QUADRATICSURFACE_HH_

#include "Surface.hh"
#include "callow/matrix/MatrixDense.hh"

namespace detran_geometry
{

/**
 *  @class QuadraticSurface
 *  @brief Defines surfaces that are second order functions of x, y, and z.
 *
 *  These surfaces are defined by
 *  @f[
 *      f(x,y,z) = Ax^2 + By^2 + Cz^2 + Dxy + Exz + Fyz + Gx + Hy + Iz + J = 0
 *  @f]
 */
class QuadraticSurface: public Surface
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef const double  c_double;
  typedef const Point   c_Point;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR AND DESTRUCTOR
  //--------------------------------------------------------------------------//

  /// Constructor
  QuadraticSurface(c_double A, c_double B, c_double C, c_double D, c_double E,
                   c_double F, c_double G, c_double H, c_double I, c_double J);

  /// Virtual destructor
  virtual ~QuadraticSurface(){}

  /// SP constructor
  static SP_surface Create(c_double A, c_double B, c_double C, c_double D,
                           c_double E, c_double F, c_double G, c_double H,
                           c_double I, c_double J);

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE --- ALL SURFACE TYPES MUST IMPLEMENT THESE
  //--------------------------------------------------------------------------//

  /// Residual function that implicitly defines the surface, i.e. f(r) = 0
  virtual double f(c_Point &r);

  /**
   *  @brief Default implementation for quadratic surfaces
   *
   *  Quadratic surfaces by definition can have at most two intersections
   *  with a given ray.  We solve this as follows.
   *
   *  Let
   *  @f[
   *      \vec{r} = [x,y,z,1]^T
   *  @f]
   *  and
   *  @f[
   *      \mathbf{M} =
   *      \left (
   *        \begin{array}{cccc}
   *         2A &  D &  E &  G \\
   *          D & 2B &  F &  H \\
   *          E &  F & 2C &  I \\
   *          G &  H &  I & 2J \\
   *        \end{array}
   *      \right ) \, .
   *  @f]
   *  Then we can express the residual function as
   *  @f[
   *      f(\vec{r}) = \frac{1}{2} \vec{r}^T \mathbf{M} \vec{r} = 0 \, .
   *  @f]
   *  Letting @f$ \vec{r} = \vec{r}_0 + t \vec{d} @f$, we have
   *  @f[
   *      2f(\vec{r}_0 + t \vec{d})
   *        = \overbrace{(\vec{d}^T   \mathbf{M} \vec{d})}^{a} t^2
   *        + \overbrace{(\vec{r}_0^T \mathbf{M} \vec{d} +
   *                      \vec{d}^T   \mathbf{M} \vec{r}_0)}^{b} t
   *        + \overbrace{\vec{r}_0^T  \mathbf{M} \vec{r}_0}^{c}
   *        = 0
   *  @f]
   *  which is just a quadratic equation in @f$ t @f$ that yields 0, 1,
   *  or two roots (intersections).
   *
   *  @param    ray     Ray
   *  @param    t_max   Optional maximum ray length to consider
   *  @return           Vector of intersection points
   */
  virtual vec_point intersections(const Ray &ray,
                                  c_double   t_max = -1);


protected:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  double d_A;
  double d_B;
  double d_C;
  double d_D;
  double d_E;
  double d_F;
  double d_G;
  double d_H;
  double d_I;
  double d_J;

  /// Coefficients in matrix form
  callow::MatrixDense d_M;

};

} // end namespace detran_geometry

#endif /* detran_geometry_QUADRATICSURFACE_HH_ */

//----------------------------------------------------------------------------//
//              end of file QuadraticSurface.hh
//----------------------------------------------------------------------------//
