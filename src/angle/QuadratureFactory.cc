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

namespace detran
{

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
