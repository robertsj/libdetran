//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  QuadraticSurface.cc
 *  @brief QuadraticSurface member definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "QuadraticSurface.hh"
#include "utilities/SoftEquivalence.hh"

namespace detran_geometry
{

//----------------------------------------------------------------------------//
QuadraticSurface::
QuadraticSurface(c_double A, c_double B, c_double C, c_double D, c_double E,
                 c_double F, c_double G, c_double H, c_double I, c_double J)
  : d_A(A), d_B(B), d_C(C),
    d_D(D), d_E(E), d_F(F),
    d_G(G), d_H(H), d_I(I),
    d_J(J),
    d_M(4, 4, 0.0)
{
  double m[] = {2*d_A,   d_D,   d_E,   d_G,
                  d_D, 2*d_B,   d_F,   d_H,
                  d_E,   d_F, 2*d_C,   d_I,
                  d_G,   d_H,   d_I, 2*d_J};
  for (size_t i = 0; i < 16; ++i) d_M[i] = m[i];
  d_M.assemble();
}

//----------------------------------------------------------------------------//
QuadraticSurface::SP_surface
QuadraticSurface::
Create(c_double A, c_double B, c_double C, c_double D, c_double E,
       c_double F, c_double G, c_double H, c_double I, c_double J)
{
  SP_surface p(new QuadraticSurface(A, B, C, D, E, F, G, H, I, J));
  return p;
}

//----------------------------------------------------------------------------//
double QuadraticSurface::f(const Point &r)
{
  double x = r.x(), y = r.y(), z = r.z();
  return d_A*x*x + d_B*y*y + d_C*z*z + d_D*x*y + d_E*x*z + d_F*y*z +
         d_G*x + d_H*y + d_I*z + d_J;
}

//----------------------------------------------------------------------------//
QuadraticSurface::vec_point
QuadraticSurface::intersections(const Ray     &r,
                                const double   max_t)
{
  double R_a[4] = {r.origin.x(),    r.origin.y(),    r.origin.z(),    1.0};
  double D_a[4] = {r.direction.x(), r.direction.y(), r.direction.z(), 0.0};
  callow::Vector R(4, R_a), D(4, D_a), M_R(4, 0.0), M_D(4, 0.0);
  d_M.multiply(R, M_R);
  d_M.multiply(D, M_D);
  double a = D.dot(M_D), b = D.dot(M_R) + R.dot(M_D), c = R.dot(M_R);
  double t[] = {0.0, 0.0};

  if (std::abs(a) < 1.0e-13)
  {
    // for a ~ 0, reduces to linear
    if (std::abs(b) > 1.0e-13) t[0] = -c / b;
  }
  else if (b*b > 4.0*a*c)
  {
    t[0] = (-b + std::sqrt(b*b - 4.0*a*c)) / (2.0*a);
    t[1] = (-b - std::sqrt(b*b - 4.0*a*c)) / (2.0*a);
  }

  vec_point points;
  bool idx[] = {t[1] <= t[0], t[1] > t[0]};
  for (size_t i = 0; i < 2; ++i)
  {
    double tt = t[idx[i]];
    if (tt <= 0.0) continue; // neglect intersections behind us
    if (max_t < 0 || tt <= max_t) points.push_back(r.origin + tt*r.direction);
  }
  return points;
}

} // end namespace detran_geometry

//----------------------------------------------------------------------------//
//              end of file QuadraticSurface.hh
//----------------------------------------------------------------------------//
