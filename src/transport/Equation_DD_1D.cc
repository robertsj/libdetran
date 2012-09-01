//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Equation_DD_1D.cc
 * \author Jeremy Roberts
 * \date   Mar 31, 2012
 * \brief  Equation_DD_1D member definitions.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

#include "Equation_DD_1D.hh"

namespace detran
{

Equation_DD_1D::Equation_DD_1D(SP_mesh mesh,
                               SP_material material,
                               SP_quadrature quadrature,
                               bool update_psi)
  :  Equation<_1D>(mesh, material, quadrature, update_psi)
  ,  d_coef_x(mesh->number_cells_x())
{
  /* ... */
}

void Equation_DD_1D::setup_group(int g)
{
  Require(g >= 0);
  Require(g < d_material->number_groups());
  d_g = g;
}

void Equation_DD_1D::setup_octant(int octant)
{
  Require(octant >= 0);
  Require(octant < 2);
  d_octant = octant;
}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of Equation_DD_1D.cc
//---------------------------------------------------------------------------//
