//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Equation_DD_2D.cc
 * \author Jeremy Roberts
 * \date   Mar 31, 2012
 * \brief  Equation_DD_2D member definitions.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#include "Equation_DD_2D.hh"

namespace detran
{

Equation_DD_2D::Equation_DD_2D(SP_mesh mesh,
                               SP_material material,
                               SP_quadrature quadrature,
                               bool store_psi)
  :  Equation(mesh, material, quadrature, store_psi)
  ,  d_coef_x(mesh->number_cells_x())
  ,  d_coef_y(mesh->number_cells_y())
{
  /* ... */
}

void Equation_DD_2D::setup_group(int g)
{
  Require(g >= 0);
  Require(g < d_material->number_groups());
}

void Equation_DD_2D::setup_octant(int octant)
{
  Require(octant >= 0);
  Require(octant < 4);
  d_octant = octant;
}

void Equation_DD_2D::setup_angle(int angle)
{
  // Currently, only the 1st octant values should be in use.
  Require(angle >= 0);
  Require(angle < d_quadrature->number_angles_octant());
  double mu  = d_quadrature->mu(0, angle);
  double eta = d_quadrature->eta(0, angle);
  for (int i = 0; i < d_mesh->number_cells_x(); i++)
  {
    d_coef_x[i] = 2.0 * mu / d_mesh->dx(i);
  }
  for (int j = 0; j < d_mesh->number_cells_y(); j++)
  {
    d_coef_y[j] = 2.0 * eta / d_mesh->dy(j);
  }
}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of Equation_DD_2D.cc
//---------------------------------------------------------------------------//
