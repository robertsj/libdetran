//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  BaseUniform.hh
 *  @brief Several classes based on uniformly-spaced abscissa
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#ifndef detran_angle_BASEUNIFORM_HH_
#define detran_angle_BASEUNIFORM_HH_

#include "angle/BaseQuadrature.hh"
#include "utilities/MathUtilities.hh"

namespace detran_angle
{

/**
 *  @class BaseUniform
 *  @brief Provides a quadrature with evenly-spaced points and equal weights
 */
class ANGLE_EXPORT BaseUniform: public BaseQuadrature
{

public:

  static std::string name() {return "uniform";}

private:

  void build_impl(c_dbl a, c_dbl b)
  {
    d_x = detran_utilities::linspace_center(a, b, d_x.size());
    double f = (b - a) / d_m;
    for (size_t i = 0; i < d_w.size(); ++i) d_w[i] = f;
  }

};

/**
 *  @class BaseUniformCosine
 *  @brief Provides a quadrature with evenly-spaced cosines and equal weights
 *
 *  This quadrature really only makes sense if the input range spans a
 *  quadrant (or two).  Here, we map from [-1, 1] to [0, pi]
 */
class ANGLE_EXPORT BaseUniformCosine: public BaseQuadrature
{

public:

  static std::string name() {return "uniformcosine";}

private:

  void build_impl(c_dbl a, c_dbl b)
  {
    using detran_utilities::pi;

    d_x = detran_utilities::linspace_center(-1.0, 1.0, d_x.size());
    for (size_t i = 0; i < d_x.size(); ++i)
      d_x[i] = std::acos(d_x[i]);
    std::sort(d_x.begin(), d_x.end());
    double f = pi / d_m;
    for (size_t i = 0; i < d_w.size(); ++i)
      d_w[i] = f;
    scale_and_shift(0.0, pi, a, b, d_x, d_w);
  }

};

/**
 *  @class BaseSimpson
 *  @brief Simpson's rule
 *
 *  This quadrature is best used for a multiple of three
 */
class ANGLE_EXPORT BaseSimpson: public BaseQuadrature
{

public:

  static std::string name() {return "simpson";}

private:

  void build_impl(c_dbl a, c_dbl b)
  {
    d_x = detran_utilities::linspace_center(0.0, 1.0, d_m);
    double w0 = 1.125 / d_m;
    double w1 = 0.750 / d_m;
    for (size_t i = 0; i < d_m - d_m % 3; ++i)
      d_w[i] = (i % 3) % 2 == 0 ? w0 : w1;
    if (d_m % 3 == 2)
      d_w[d_m-2] = 1.0 / d_m;
    else if (d_m % 3)
      d_w[d_m-1] = 1.0 / d_m;
    scale_and_shift(0.0, 1.0, a, b, d_x, d_w);
  }

};

} // end namespace detran_angle

#endif /* BASEUNIFORM_HH_ */
