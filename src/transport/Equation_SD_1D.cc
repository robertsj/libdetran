//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Equation_SD_1D.cc
 * \author robertsj
 * \date   Jun 9, 2012
 * \brief  Equation_SD_1D class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// Detran
#include "Equation_SD_1D.hh"

namespace detran
{

Equation_SD_1D::Equation_SD_1D(SP_mesh mesh,
                               SP_material material,
                               SP_quadrature quadrature,
                               bool update_psi)
  :  Equation(mesh, material, quadrature, update_psi)
  ,  d_coef_x(mesh->number_cells_x())
{
  /* ... */
}

void Equation_SD_1D::setup_group(int g)
{
  Require(g >= 0);
  Require(g < d_material->number_groups());
  d_g = g;
}

void Equation_SD_1D::setup_octant(int octant)
{
  Require(octant >= 0);
  Require(octant < 2);
  d_octant = octant;
}

void Equation_SD_1D::setup_angle(int angle)
{
  // Currently, only the 1st octant values should be in use.
  Require(angle >= 0);
  Require(angle < d_quadrature->number_angles_octant());
  double mu  = d_quadrature->mu(0, angle);
  for (int i = 0; i < d_mesh->number_cells_x(); i++)
  {
    d_coef_x[i] = mu / d_mesh->dx(i);
  }
  d_angle = angle;
}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of Equation_SD_1D.cc
//---------------------------------------------------------------------------//



