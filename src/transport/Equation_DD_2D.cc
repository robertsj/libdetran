//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Equation_DD_2D.cc
 *  @author Jeremy Roberts
 *  @date   Mar 31, 2012
 *  @brief  Equation_DD_2D member definitions.
 */
//---------------------------------------------------------------------------//

#include "Equation_DD_2D.hh"

namespace detran
{

//---------------------------------------------------------------------------//
Equation_DD_2D::Equation_DD_2D(SP_mesh       mesh,
                               SP_material   material,
                               SP_quadrature quadrature,
                               const bool    update_psi)
  :  Equation<_2D>(mesh, material, quadrature, update_psi)
  ,  d_coef_x(mesh->number_cells_x())
  ,  d_coef_y(mesh->number_cells_y())
{
  /* ... */
}

//---------------------------------------------------------------------------//
void Equation_DD_2D::setup_group(const size_t g)
{
  Require(g < d_material->number_groups());
  d_g = g;
}

//---------------------------------------------------------------------------//
void Equation_DD_2D::setup_octant(const size_t octant)
{
  Require(octant < 4);
  d_octant = octant;
}

//---------------------------------------------------------------------------//
void Equation_DD_2D::setup_angle(const size_t angle)
{
  Require(angle < d_quadrature->number_angles_octant());
  d_angle = angle;
  double mu  = d_quadrature->mu(0, d_angle);
  double eta = d_quadrature->eta(0, d_angle);
  for (size_t i = 0; i < d_mesh->number_cells_x(); ++i)
  {
    d_coef_x[i] = 2.0 * mu / d_mesh->dx(i);
  }
  for (size_t j = 0; j < d_mesh->number_cells_y(); ++j)
  {
    d_coef_y[j] = 2.0 * eta / d_mesh->dy(j);
  }
}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of Equation_DD_2D.cc
//---------------------------------------------------------------------------//
