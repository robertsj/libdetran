//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   QuadratureFactory.cc
 *  @author robertsj
 *  @date   Apr 11, 2012
 *  @brief  QuadratureFactory member definitions.
 */
//---------------------------------------------------------------------------//

// Detran
#include "QuadratureFactory.hh"
// 1D
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
#include "Collocated.hh"
#include "Uniform.hh"

namespace detran_angle
{

// Build a quadrature
QuadratureFactory::SP_quadrature QuadratureFactory::
build(SP_input input, const int dimension)
{
  SP_quadrature q;
  build(q, input, dimension);
  return q;
}
// Build a quadrature
void QuadratureFactory::
build(SP_quadrature &q, SP_input input, const int dimension)
{
  using std::string;

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
    cheb_norm = input->get<int>("quad_chebyshev_normalize");

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
  // 1D QUADRATURES
  //-------------------------------------------------------------------------//

  if (quad_type == "gausslegendre")
  {
    Insist(dimension == 1, "GaussLegendre only for 1D.");
    q = new GaussLegendre(np);
  }
  else if (quad_type == "gausschebyshev")
  {
    Insist(dimension == 1, "GaussLegendre only for 1D.");
    q = new GaussChebyshev(np, cheb_norm);
  }
  else if (quad_type == "dpn")
  {
    Insist(dimension == 1, "DPN only for 1D.");
    q = new DPN(np);
  }
  else if (quad_type == "dtn")
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

  else if (quad_type == "collocated")
  {
    Insist(dimension > 1, "Collocated only for 2D or 3D.");
    q = new Collocated(dimension, azimuths_octant, 1, polar_octant, polar_type);
  }
  else if (quad_type == "uniform")
  {
    Insist(dimension > 1, "Uniform only for 2D or 3D.");
    int num_space = 10;
    if (input->check("quad_uniform_number_space"))
      num_space = input->get<int>("quad_uniform_number_space");
    q = new Uniform(dimension, azimuths_octant, num_space, polar_octant, polar_type);
  }
  else
  {
    THROW("Unsupported quadrature selected: " + quad_type);
  }

}

void QuadratureFactory::help() const
{

  std::cout << "Available quadratures are: " << std::endl;
  std::cout << "  gausslegendre (1D)" << std::endl;
  std::cout << "  gausschebyshev (1D)" << std::endl;
  std::cout << "  dpn (1D)" << std::endl;
  std::cout << "  dtn (1D)" << std::endl;
  std::cout << "  levelsymmetric" << std::endl;
  std::cout << "  quadruplerange" << std::endl;
  std::cout << "  uniformequal" << std::endl;
  std::cout << "  chebyshevlegendre" << std::endl;
  std::cout << "  chebyshevdpn" << std::endl;

}

} // end namespace detran_angle

//---------------------------------------------------------------------------//
//              end of QuadratureFactory.cc
//---------------------------------------------------------------------------//
