//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  QuadratureFactory.cc
 *  @brief QuadratureFactory member definitions
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// Detran
#include "QuadratureFactory.hh"
// 1D
#include "PolarQuadrature.hh"
#include "GaussLegendre.hh"
#include "GaussChebyshev.hh"
#include "DPN.hh"
#include "DTN.hh"
// 2D/3D
#include "LevelSymmetric.hh"
#include "QuadrupleRange.hh"
#include "UniformEqual.hh"
// Product
#include "ChebyshevLegendre.hh"
#include "ChebyshevDPN.hh"
#include "LegendreDTN.hh"
// MOC
//#include "Collocated.hh"
//#include "Uniform.hh"

#include "BaseGL.hh"
#include "BaseCL.hh"
#include "BaseUniform.hh"
#include "TabuchiYamamoto.hh"

namespace detran_angle
{

//----------------------------------------------------------------------------//
QuadratureFactory::SP_quadrature QuadratureFactory::
build(SP_input input, const int dimension)
{
  using std::string;

  SP_quadrature q;

  // Set the quadrature type.
  string quad_type = "notset";
  if (!input->check("quad_type"))
  {
    quad_type = "chebyshevdpn";
    if (dimension == 1) quad_type = "gausslegendre";
  }
  else
  {
    quad_type = input->get<string>("quad_type");
  }

  // Set default values.
  bool moc = false;
  int azimuths_octant = 2;
  int polar_octant = 1;
  string polar_type = "TY";

  // All quadratures are parameterized by number of polar and/or
  // number of azimuths
  int np = 1;
  int na = 1;

  // For quadratures with Chebyshev components, the default is
  // not to normalize the weights.
  bool cheb_norm = false;
  if (input->check("quad_chebyshev_normalize"))
    cheb_norm = (0 != input->get<int>("quad_chebyshev_normalize"));

  // First, check whether this is an MOC quadrature or not.  If it
  // is, we need different parameters.
  if (quad_type == "uniform" || quad_type == "collocated")
  {

    // We're using an MOC quadrature (but that doesn't mean this
    // is an MOC problem.
    moc = true;

    // Set the number of azimuths per octant
    if (input->check("quad_azimuths_octant"))
      azimuths_octant = input->get<int>("quad_azimuths_octant");

    // Set the number of polar angles per octant
    if (input->check("quad_polar_octant"))
      polar_octant = input->get<int>("quad_polar_octant");

    // Set the polar type
    if (input->check("polar_type"))
      polar_type = input->get<string>("quad_polar_type");

  }
  else
  {
    // Get the number of polar and/or azimuths
    if (input->check("quad_number_polar_octant"))
      np = input->get<int>("quad_number_polar_octant");
    if (input->check("quad_number_azimuth_octant"))
      na = input->get<int>("quad_number_azimuth_octant");
  }

  //-------------------------------------------------------------------------//
  // 1D (i.e. POLAR) QUADRATURES
  //-------------------------------------------------------------------------//

  if (quad_type == "gl")
  {
    Insist(dimension == 1, "GaussLegendre only for 1D.");
    q = new PolarQuadrature<BaseGL>(np);
  }
  else if (quad_type == "gc")
  {
    Insist(dimension == 1, "GaussLegendre only for 1D.");
    q = new GaussChebyshev(np, cheb_norm);
  }
  else if (quad_type == "dgl")
  {
    Insist(dimension == 1, "DPN only for 1D.");
    q = new DPN(np);
  }
  else if (quad_type == "dgc")
  {
    Insist(dimension == 1, "DTN only for 1D.");
    q = new DTN(np, cheb_norm);
  }

  //-------------------------------------------------------------------------//
  // 2D/3D QUADRATURES
  //-------------------------------------------------------------------------//

  else if (quad_type == "levelsymmetric")
  {
    Insist(dimension > 1, "LevelSymmetric only for 2D or 3D.");
    q = new LevelSymmetric(np, dimension);
  }
  else if (quad_type == "uniformequal")
  {
    Insist(dimension > 1, "UniformEqual only for 2D or 3D.");
    q = new UniformEqual(np, dimension);
  }

  else if (quad_type == "chebyshevlegendre")
  {
    Insist(dimension > 1, "ChebyshevLegendre only for 2D or 3D.");
    q = new ChebyshevLegendre(dimension, na, np);
  }
  else if (quad_type == "chebyshevdpn")
  {
    Insist(dimension > 1, "ChebyshevDPN only for 2D or 3D.");
    q = new ChebyshevDPN(dimension, na, np);
  }
  else if (quad_type == "legendredtn")
  {
    Insist(dimension > 1, "LegendreDTN only for 2D or 3D.");
    q = new LegendreDTN(dimension, na, np);
  }
  else if (quad_type == "quadruplerange")
  {
    Insist(dimension > 1, "QuadrupleRange only for 2D or 3D.");
    q = new QuadrupleRange(dimension, na, np);
  }

  //-------------------------------------------------------------------------//
  // MOC QUADRATURES
  //-------------------------------------------------------------------------//
//
//  else if (quad_type == "collocated")
//  {
//    Insist(dimension > 1, "Collocated only for 2D or 3D.");
//    q = new Collocated(dimension, azimuths_octant, 1, polar_octant, polar_type);
//  }
//  else if (quad_type == "uniform")
//  {
//    Insist(dimension > 1, "Uniform only for 2D or 3D.");
//    int num_space = 10;
//    if (input->check("quad_uniform_number_space"))
//      num_space = input->get<int>("quad_uniform_number_space");
//    q = new Uniform(dimension, azimuths_octant, num_space, polar_octant, polar_type);
//  }
  else
  {
    THROW("Unsupported quadrature selected: " + quad_type);
  }

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
    q = new BaseCL();
  else if (type == "dgc")
    q = new BaseDCL();
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
