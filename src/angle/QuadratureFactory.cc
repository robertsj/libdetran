//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  QuadratureFactory.cc
 *  @brief QuadratureFactory member definitions
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "QuadratureFactory.hh"
// Base
#include "BaseGL.hh"
#include "BaseGC.hh"
#include "BaseUniform.hh"
#include "TabuchiYamamoto.hh"
// 1D
#include "PolarQuadrature.hh"
// 2D/3D
#include "ProductQuadratureAP.hh"
#include "LevelSymmetric.hh"
#include "UniformEqual.hh"
//
#include "UserQuadrature.hh"

namespace detran_angle
{

//----------------------------------------------------------------------------//
QuadratureFactory::SP_quadrature QuadratureFactory::
build(SP_input input, const int dimension)
{
  using std::string;

  // Set the quadrature type.
  string quad_type;
  if (!input->check("quad_type"))
  {
    quad_type = "u-dgl";
    if (dimension == 1) quad_type = "gl";
  }
  else
  {
    quad_type = input->get<string>("quad_type");
  }

  // All quadratures are parameterized by number of polar and/or azimuth angles
  int np = 1;
  int na = 1;

  // For quadratures with Chebyshev components, the default is
  // not to normalize the weights.
  bool cheb_norm = false;
  if (input->check("quad_chebyshev_normalize"))
    cheb_norm = (0 != input->get<int>("quad_chebyshev_normalize"));

  // Get the number of polar and/or azimuths
  if (input->check("quad_number_polar_octant"))
    np = input->get<int>("quad_number_polar_octant");
  if (input->check("quad_number_azimuth_octant"))
    na = input->get<int>("quad_number_azimuth_octant");

  SP_quadrature q;
  if (dimension == 1)
  {
    //------------------------------------------------------------------------//
    // POLAR QUADRATURES
    //------------------------------------------------------------------------//

    if      (quad_type == "gl")
      q = new PolarGL(np);
    else if (quad_type == "dgl")
      q = new PolarDGL(np);
    else if (quad_type == "gc")
      q = new PolarGC(np, cheb_norm);
    else if (quad_type == "dgc")
      q = new PolarDGC(np, cheb_norm);
    else if (quad_type == "uniform")
      q = new PolarU(np);
    else if (quad_type == "uniformcosine")
      q = new PolarUC(np);
    else if (quad_type == "ty")
      q = new PolarTY(np);
    else if (quad_type == "asdr")
      q = new PolarASDR(np);
    else if (quad_type == "tg")
      q = new PolarTG(np);
  }
  else
  {
    //------------------------------------------------------------------------//
    // PRODUCT QUADRATURES
    //------------------------------------------------------------------------//

    if      (quad_type == "gl-gl")
      q = new Product_GL_GL(dimension, na, np);
    else if (quad_type == "gl-dgl")
      q = new Product_GL_DGL(dimension, na, np);
    else if (quad_type == "gl-u")
      q = new Product_GL_U(dimension, na, np);
    else if (quad_type == "gl-ty")
      q = new Product_GL_TY(dimension, na, np);
    //
    else if (quad_type == "dgl-gl")
      q = new Product_DGL_GL(dimension, na, np);
    else if (quad_type == "dgl-dgl")
      q = new Product_DGL_DGL(dimension, na, np);
    else if (quad_type == "dgl-u")
      q = new Product_DGL_U(dimension, na, np);
    else if (quad_type == "dgl-ty")
      q = new Product_DGL_TY(dimension, na, np);
    //
    else if (quad_type == "u-gl")
      q = new Product_U_GL(dimension, na, np);
    else if (quad_type == "u-dgl")
      q = new Product_U_DGL(dimension, na, np);
    else if (quad_type == "u-gc")
      q = new Product_U_GC(dimension, na, np);
    else if (quad_type == "u-dgc")
      q = new Product_U_DGC(dimension, na, np);
    else if (quad_type == "u-u")
      q = new Product_U_U(dimension, na, np);
    else if (quad_type == "u-ty")
      q = new Product_U_TY(dimension, na, np);
    else if (quad_type == "u-asdr")
      q = new Product_U_ASDR(dimension, na, np);
    else if (quad_type == "u-tg")
      q = new Product_U_TG(dimension, na, np);
    //
    else if (quad_type == "s-gl")
      q = new Product_S_GL(dimension, na, np);
    else if (quad_type == "s-dgl")
      q = new Product_S_DGL(dimension, na, np);
    else if (quad_type == "s-u")
      q = new Product_S_U(dimension, na, np);
    else if (quad_type == "s-ty")
      q = new Product_S_TY(dimension, na, np);
    //
    else if (quad_type == "asqr-asdr")
      q = new Product_ASQR_ASDR(dimension, na, np);
    else if (quad_type == "asqr-gl")
      q = new Product_ASQR_GL(dimension, na, np);
    else if (quad_type == "asqr-dgl")
      q = new Product_ASQR_DGL(dimension, na, np);
    else if (quad_type == "asqr-gc")
      q = new Product_ASQR_GC(dimension, na, np);
    else if (quad_type == "asqr-dgc")
      q = new Product_ASQR_DGC(dimension, na, np);
    else if (quad_type == "asqr-tg")
      q = new Product_ASQR_TG(dimension, na, np);

    //------------------------------------------------------------------------//
    // OTHER QUADRATURES
    //------------------------------------------------------------------------//

    if      (quad_type == "levelsymmetric")
      q = new LevelSymmetric(np, dimension);
    else if (quad_type == "uniformequal")
      q = new UniformEqual(np, dimension);

  }

  if (quad_type == "user")
  {
    Quadrature::vec_dbl mu_v, eta_v, wt_v;
    Insist(input->check("quad_mu"), "quad_mu required for user quadrature");
    Insist(input->check("quad_wt"), "quad_wt required for user quadrature");
    mu_v = input->get<Quadrature::vec_dbl>("quad_mu");
    if (input->check("quad_eta"))
      eta_v = input->get<Quadrature::vec_dbl>("quad_eta");
    wt_v = input->get<Quadrature::vec_dbl>("quad_wt");
    q = new UserQuadrature(dimension, mu_v, eta_v, wt_v);
  }

  Insist(q, "Unsupported quadrature type: " + quad_type);
  return q;
}

//----------------------------------------------------------------------------//
QuadratureFactory::SP_basequadrature
QuadratureFactory::build_base(const std::string &type,
                              const size_t       n,
                              const double       a,
                              const double       b)
{
  Require(n > 0);
  SP_basequadrature q;
  if (type == "gl")
    q = new BaseGL();
  else if (type == "dgl")
    q = new BaseDGL();
  else if (type == "gc")
    q = new BaseGC();
  else if (type == "dgc")
    q = new BaseDGC();
  else if (type == "uniform")
    q = new BaseUniform();
  else if (type == "uniformcosine")
    q = new BaseUniformCosine();
  else if (type == "simpson")
    q = new BaseSimpson();
  else if (type == "ty")
    q = new TabuchiYamamoto();
  else
    THROW("Unsupported base quadrature type: " + type);
  q->build(a, b, n);
  return q;
}

} // end namespace detran_angle

//----------------------------------------------------------------------------//
//              end of QuadratureFactory.cc
//----------------------------------------------------------------------------//
