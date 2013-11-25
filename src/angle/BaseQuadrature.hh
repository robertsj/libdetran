//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  BaseQuadrature.hh
 *  @brief BaseQuadrature class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_angle_BASEQUADRATURE_HH_
#define detran_angle_BASEQUADRATURE_HH_

#include "angle/angle_export.hh"
#include "utilities/DBC.hh"
#include "utilities/Definitions.hh"
#include "utilities/MathUtilities.hh"
#include "utilities/SP.hh"
#include <cstdio>

namespace detran_angle
{

/**
 *  @class BaseQuadrature
 *  @brief Base quadrature interface for 1-D integration
 */
class ANGLE_EXPORT BaseQuadrature
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<BaseQuadrature>    SP_basequadrature;
  typedef detran_utilities::vec_dbl               vec_dbl;
  typedef detran_utilities::size_t                size_t;
  typedef const double                            c_dbl;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  virtual ~BaseQuadrature(){}

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /**
   *  @brief Build the quadrature
   *  @param    a     lower bound
   *  @param    b     upper bound
   *  @param    m     number of abscissa
   */
  void build(c_dbl a, c_dbl b, const size_t m)
  {
    Require(m > 0);
    Require(b > a);

    d_m = m;
    d_x.resize(d_m, 0.0);
    d_w.resize(d_m, 0.0);

    build_impl(a, b);

    Ensure(d_x.size() == d_w.size());
    Ensure(d_x.size() == d_m);
    for (size_t i = 1; i < d_m; ++i)
    {
      Ensure(d_x[i] > d_x[i - 1]);
    }
  }

  /// Get size of quadrature
  size_t size() const { return d_m; }

  /// Get the abscissa
  vec_dbl& get_x() { return d_x; }
  const vec_dbl& get_x() const { return d_x; }

  /// Get the weights
  vec_dbl& get_w() { return d_w; }
  const vec_dbl& get_w() const { return d_w; }

  /// Pretty print
  void display() const
  {
    using std::cout;
    using std::endl;
    using std::printf;
    cout << "   n             x                   w       " << endl;
    cout << "  ---   ------------------  -----------------" << endl;
    for (size_t n = 0; n < d_m; ++n)
      printf("%4i    %16.13f   %16.13f   \n", n, d_x[n], d_w[n] );
    printf("Total weight:  %16.13f\n", detran_utilities::vec_sum(d_w));
  }

  /// Return the quadrature name.  Subclasses should override.
  static std::string name() {return "base";}

protected:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Number of abscissa
  size_t  d_m;
  /// Vector of abscissa
  vec_dbl d_x;
  /// Vector of weights
  vec_dbl d_w;

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

  /**
   *  @brief Scale and/or shift abscissa and weights
   *  @param    a0    internal lower bound
   *  @param    b0    internal upper bound
   *  @param    a1    user-given lower bound
   *  @param    b1    user-given upper bound
   *  @param    x     abscissa to adjust
   *  @param    w     weights to adjust
   */
  void scale_and_shift(c_dbl a0, c_dbl b0, c_dbl a1, c_dbl b1,
                       vec_dbl &x, vec_dbl &w)
  {
    Require(b0 > a0);
    Require(b1 > a1);
    double scale = (b1 - a1) / (b0 - a0);
    double shift = a1 - a0 * scale;
    detran_utilities::vec_scale(x, scale);
    for (size_t i = 0; i < x.size(); ++i) x[i] += shift;
    detran_utilities::vec_scale(w, scale);
  }

private:

  /// Implementation of the build
  virtual void build_impl(c_dbl a, c_dbl b) = 0;

};

} // end namespace detran_angle

#endif /* detran_angle_BASEQUADRATURE_HH_ */

//----------------------------------------------------------------------------//
//              end of BaseQuadrature.hh
//----------------------------------------------------------------------------//
