//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  BaseGL.hh
 *  @brief BaseGL (and BaseDGL) class definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_angle_BASEGL_HH_
#define detran_angle_BASEGL_HH_

#include "angle/BaseQuadrature.hh"
#include "utilities/Constants.hh"

namespace detran_angle
{

/**
 *  @class BaseGL
 *  @brief Base Gauss-Legendre generator
 *
 *  The Gauss-Legendre quadrature approximates
 *  @f[
 *      \int^{1}_{-1} f(x) dx
 *        \approx \sum^m_{i=1} W_m f(x_i) \, .
 *  @f]
 *
 *  The abscissa @f$ x_i @f$ of the m-point quadrature
 *  turn out to be the zeroes of the
 *  Legendre polynomial of degree @f$ m + 1 @f$.  The weights
 *  can be found numerically, as done in this implementation,
 *  which is a modified version of the function
 *  gauleg from <B>Numerical Recipes in C</B>
 *  by W.H. Press, S.A. Teukolsky, W.T. Vetterling, and
 *  & B.P. Flannery, (Cambridge Univ. Press)
 */
class BaseGL: public BaseQuadrature
{

public:

  static std::string name() {return "gl";}

protected:

  void generate_gl_parameters(vec_dbl &x, vec_dbl &w)
  {
    size_t j, i, m = x.size();
    double z1, z, pp, p3, p2, p1;
    // The roots are symmetric, so we only find half of them.
    for (i = 1; i <= m; ++i)
    {
      // Loop over the desired roots.
      z = std::cos(detran_utilities::pi * (i - 0.25) / (m + 0.5));
      // Starting with the above approximation to the ith root, we enter
      // the main loop of refinement by Newon's method.
      int count = 0;
      do
      {
        p1 = 1.0;
        p2 = 0.0;
        // Recurrence to get Legendre polynomial.
        for (j = 1; j <= m; j++)
        {
          p3 = p2;
          p2 = p1;
          p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
        }
        // p1 is now the desired Legendre polynomial. We next compute
        // pp, its derivative, by a standard relation involving also
        // p2, the polynomial of one lower order.
        pp = m * (z * p1 - p2) / (z * z - 1.0);
        z1 = z;
        z = z1 - p1 / pp; // <-- Newton's method
      } while (std::abs(z - z1) > 1e-15 && ++count < 200);
      // Store the root and compute the weight.
      x[i - 1] = -z; // negative puts it in increasing order
      w[i - 1] = 2.0 / ((1.0 - z * z) * pp * pp);
    }
  }

private:

  void build_impl(c_dbl a, c_dbl b)
  {
    // Generate base parameters
    generate_gl_parameters(d_x, d_w);
    // Scale and shift
    scale_and_shift(-1, 1, a, b, d_x, d_w);
  }

};

/**
 *  @class BaseDGL
 *  @brief Base double Gauss-Legendre generator
 *
 *  The double Gauss-Legendre quadrature breaks the domain into two halves in
 *  each of which an independent set of G-L parameters is used.
 */
class BaseDGL: public BaseGL
{

public:

  static std::string name() {return "dgl";}

private:

  void build_impl(c_dbl a, c_dbl b)
  {
    Insist(d_m % 2 == 0, "Double G-L requires an even number of points.");
    vec_dbl xh(d_m/2, 0.0);
    vec_dbl wh(d_m/2, 0.0);
    generate_gl_parameters(xh, wh);
    scale_and_shift(-1.0, 1.0, 0.0, 1.0, xh, wh);
    for (size_t i = 0; i < d_m/2; ++i)
    {
      size_t j = i + d_m/2;
      size_t k = d_m/2 - i - 1;
      d_x[j] =  xh[i];
      d_w[j] =  wh[i];
      d_x[k] = -xh[i];
      d_w[k] =  wh[i];
    }
    scale_and_shift(-1.0, 1.0, a, b, d_x, d_w);
  }

};

} // end namespace detran_angle

#endif /* detran_angle_BASEQUADRATURE_HH_ */

//----------------------------------------------------------------------------//
//              end of BaseGL.hh
//----------------------------------------------------------------------------//
