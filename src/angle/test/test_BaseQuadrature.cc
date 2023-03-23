//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_BaseQuadrature.cc
 *  @brief Test of various base quadrature classes
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 *
 *  Ideally, any 1-D based quadratures will be tested here.
 */
//----------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "angle/QuadratureFactory.hh"

using namespace detran_angle;
using namespace detran_utilities;
using namespace std;

/*    integrand          range    value
 *  ----------------------------------------------------------------------
 *     0.5              (-1, 1)    1.00
 *     H(x)             (-1, 1)    1.00
 *     x^3              ( 0, 1)    0.25
 *     exp(x^2)*sin(x)  ( 0, 1)    0.77874516052672400466
 */
int number_tests = 4;
double ref_val[]   = {1.0, 1.0, 0.25, 0.77874516052672400466};
double lower[] = {-1.0,-1.0, 0.0, 0.0};
double upper[] = { 1.0, 1.0, 1.0, 1.0};
double integrand(const double x, const int id)
{
  double y = 0;
  if (id == 0)
    y = 0.5;
  else if (id == 1)
    y = x > 0 ? 1.0 : 0.0;
  else if (id == 2)
    y = x * x * x;
  else if (id == 3)
    y = std::exp(x*x) * std::sin(x);
  return y;
}

int number_types = 7;
const char* types[] = {"gl", "dgl", "gc", "dgc",
                       "uniform", "uniformcosine", "simpson"};

//----------------------------------------------------------------------------//
TEST(BaseQuadrature, Regression)
{
  // All available 1-D quadratures are tested on three simple
  // problems to ensure the expected value for each is computed.
  vec2_dbl val(number_types, vec_dbl(number_tests, 0.0));
  vec2_dbl err(number_types, vec_dbl(number_tests, 0.0));
  QuadratureFactory::SP_basequadrature Q;
  for (int i = 0; i < number_types; ++i)
  {
    std::cout << "testing " << types[i] << std::endl;
    for (int j = 0; j < number_tests; ++j)
    {
//      std::cout << "  problem " << j;
      Q = QuadratureFactory::build_base(types[i], 8, lower[j], upper[j]);
      for (int k = 0; k < Q->size(); ++k)
        val[i][j] += Q->get_w()[k] * integrand(Q->get_x()[k], j);
      err[i][j] = val[i][j] - ref_val[j];
//      std::cout << "  value = " << val[i][j]
//                << "  error = " << err[i][j]
//                << "  sum w = " << vec_sum(Q->get_w()) << std::endl;
      std::printf(" %4i %16.13f  %16.13f  %16.13f  %16.13f \n",
                   j, ref_val[j], val[i][j], err[i][j], vec_sum(Q->get_w()));
    }
  }
}

TEST(BaseQuadrature, GL)
{
  {
    QuadratureFactory::SP_basequadrature Q;
    Q = QuadratureFactory::build_base("gl", 6, 0.0, 1.0);
    Q->display();
//    double v_0 = 0.0, v_1 = 0.0;
//    for (int i = 0; i < Q.size(); ++i)
//    {
//      double x = Q.get_x()[i];
//      double w = Q.get_w()[i];
//      std::cout << x << " " << w << std::endl;
//      v_0 += w * x * x * x;
//      v_1 += w * std::exp(x * x) * std::sin(x);
//    }
//    std::cout << " v_0=" << v_0 << " v_1=" << v_1 << std::endl;
//  }
  }
}

TEST(BaseQuadrature, CL)
{
}

TEST(BaseQuadrature, Uniform)
{
  using detran_utilities::pi;
  QuadratureFactory::SP_basequadrature Q;
  Q = QuadratureFactory::build_base("uniformcosine", 6,
                                    0.0, 1.0);
  Q->display();
}

//---------------------------------------------------------------------------//
//              end of test_BaseQuadrature.cc                                //
//---------------------------------------------------------------------------//
