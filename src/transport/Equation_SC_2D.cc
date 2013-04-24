//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Equation_SC_2D.cc
 *  @author Jeremy Roberts
 *  @date   Mar 31, 2012
 *  @brief  Equation_SC_2D member definitions.
 */
//---------------------------------------------------------------------------//

#include "Equation_SC_2D.hh"

namespace detran
{

//---------------------------------------------------------------------------//
Equation_SC_2D::Equation_SC_2D(SP_mesh       mesh,
                               SP_material   material,
                               SP_quadrature quadrature,
                               const bool    update_psi)
  :  Equation<_2D>(mesh, material, quadrature, update_psi)
  ,  d_alpha(mesh->number_cells_x())
  ,  d_beta(mesh->number_cells_y())
{
  /* ... */
}

//---------------------------------------------------------------------------//
void Equation_SC_2D::setup_group(const size_t g)
{
  Require(g >= 0);
  Require(g < d_material->number_groups());
  d_g = g;
}

//---------------------------------------------------------------------------//
void Equation_SC_2D::setup_octant(const size_t octant)
{
  Require(octant >= 0);
  Require(octant < 4);
  d_octant = octant;
}

//---------------------------------------------------------------------------//
void Equation_SC_2D::setup_angle(const size_t angle)
{
  Require(angle < d_quadrature->number_angles_octant());
  d_angle = angle;
  double mu  = d_quadrature->mu(0, d_angle);
  double eta = d_quadrature->eta(0, d_angle);
  for (size_t i = 0; i < d_mesh->number_cells_x(); ++i)
  {
    d_alpha[i] = d_mesh->dx(i) / mu;
  }
  for (size_t j = 0; j < d_mesh->number_cells_y(); ++j)
  {
    d_beta[j] = d_mesh->dy(j) / eta;
  }
}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of Equation_SC_2D.cc
//---------------------------------------------------------------------------//
