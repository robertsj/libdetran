//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  BaseGC.hh
 *  @brief BaseGC (and BaseDGC) class definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_angle_BASEGC_HH_
#define detran_angle_BASEGC_HH_

#include "angle/BaseQuadrature.hh"

namespace detran_angle
{

/**
 *  @class BaseGC
 *  @brief Base Gauss-Chebyshev generator
 *
 *  The Gauss-Chebyshev quadrature approximates
 *  @f[
 *      \int^{1}_{-1} \frac{f(x)}{\sqrt{1-x^2}} dx
 *        \approx \sum^m_{i=1} W_i f(x_i) \, .
 *  @f]
 *  Now, the integral we want is not actually weighted, and so we
 *  set the true weights to be
 *  @f[
 *     w_i = W_i \sqrt{1-x_i^2} \, .
 *  @f]
 *  Note that these weights are not normalized exactly, meaning the constant
 *  function cannot be integrated exactly.  Where appropriate, the client
 *  should normalize as needed.
 *
 *  The abscissa of the m-point quadrature are the zeros of the
 *  Chebyshev polynomial of degree m + 1.
 *
 */
class ANGLE_EXPORT BaseGC: public BaseQuadrature
{

public:

  static std::string name() {return "cl";}

protected:

  void generate_gc_parameters(vec_dbl &x, vec_dbl &w)
  {
    const double pi = detran_utilities::pi;
    size_t m = x.size();

    // Base weight (constant for a particular order)
    double W = pi / (double)m;

    // Generate abscissa
    for (size_t i = 1; i <= m; ++i)
    {
      x[i-1] = -std::cos((i-0.5)*pi/m); // negative puts into increasing order
    }

    // Generate weights
    double w_tot = 0.0;
    for (size_t i = 0; i < m; ++i)
    {
      w[i] = W * std::sqrt(1.0 - x[i]*x[i]);
      w_tot += w[i];
    }
  }

private:

  void build_impl(c_dbl a, c_dbl b)
  {
    // Generate base parameters
    generate_gc_parameters(d_x, d_w);
    // Scale and shift
    scale_and_shift(-1.0, 1.0, a, b, d_x, d_w);
  }

};

/**
 *  @class BaseDGC
 *  @brief Base double Gauss-Chebyshev generator
 *
 *  The double Gauss-Chebyshev quadrature breaks the domain into two halves in
 *  each of which an independent set of G-C parameters is used.
 */
class ANGLE_EXPORT BaseDGC: public BaseGC
{

public:

  static std::string name() {return "cl";}

private:

  void build_impl(c_dbl a, c_dbl b)
  {
    Insist(d_m % 2 == 0, "Double G-C requires an even number of points.");
    vec_dbl xh(d_m/2, 0.0);
    vec_dbl wh(d_m/2, 0.0);
    generate_gc_parameters(xh, wh);
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

#endif /* detran_angle_BASEGC_HH_ */

//----------------------------------------------------------------------------//
//              end of BaseGC.hh
//----------------------------------------------------------------------------//
