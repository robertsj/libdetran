//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   QuadratureFactory.cc
 * \author robertsj
 * \date   Apr 11, 2012
 * \brief  QuadratureFactory member definitions.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// Detran
#include "QuadratureFactory.hh"
// 1D
#include "GaussLegendre.hh"
// 2D/3D
#include "QuadrupleRange.hh"
#include "LevelSymmetric.hh"
#include "UniformEqual.hh"
// MOC
#include "Collocated.hh"
#include "Uniform.hh"

namespace detran
{

// Build a quadrature
void QuadratureFactory::
build(SP_quadrature &q, SP_input input, int dimension)
{
  using std::string;

  // Set the quadrature type.
  string quad_type;
  if (!input->check("quad_type"))
  {
    if (dimension == 1) quad_type = "gausslegendre";
    if (dimension == 2) quad_type = "quadruplerange";
    if (dimension == 3) quad_type = "levelsymmetric";
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
  int quad_order = 2;

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

    // Set the quadrature order.  Note, this means different
    // things for each quadrature.  See the documentation.
    if (input->check("quad_order"))
      quad_order = input->get<int>("quad_order");

  }

  // SN quadratures
  if (quad_type == "gausslegendre")
  {
    Insist(dimension == 1, "GaussLegendre only for 1D.");
    q = new GaussLegendre(quad_order);
  }
  else if (quad_type == "quadruplerange")
  {
    Insist(dimension > 1, "QuadrupleRange only for 2D or 3D.");
    q = new QuadrupleRange(quad_order, dimension);
  }
  else if (quad_type == "levelsymmetric")
  {
    Insist(dimension > 1, "LevelSymmetric only for 2D or 3D.");
    q = new LevelSymmetric(quad_order, dimension);
  }
  else if (quad_type == "uniformequal")
  {
    Insist(dimension > 1, "UniformEqual only for 2D or 3D.");
    q = new UniformEqual(quad_order, dimension);
  }
  // MOC quadratures
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
    THROW("Unsupported quadrature selected.")
  }

}

// Build a quadrature
void QuadratureFactory::
build(SP_quadrature &q, std::string type, int order, int dimension)
{
  if (type == "gausslegendre")
  {
    Insist(dimension == 1, "GaussLegendre only for 1D.");
    q = new GaussLegendre(order);
  }
  else if (type == "quadruplerange")
  {
    Insist(dimension > 1, "QuadrupleRange only for 2D or 3D.");
    q = new QuadrupleRange(order, dimension);
  }
  else if (type == "levelsymmetric")
  {
    Insist(dimension > 1, "LevelSymmetric only for 2D or 3D.");
    q = new LevelSymmetric(order, dimension);
  }
  else if (type == "uniformequal")
  {
    Insist(dimension > 1, "UniformEqual only for 2D or 3D.");
    q = new UniformEqual(order, dimension);
  }
  else
  {
    THROW("Unsupported quadrature selected.")
  }
}



} // end namespace detran

//---------------------------------------------------------------------------//
//              end of QuadratureFactory.cc
//---------------------------------------------------------------------------//
