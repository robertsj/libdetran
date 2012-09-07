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
  :  Equation<_1D>(mesh, material, quadrature, update_psi)
  ,  d_coef_x(mesh->number_cells_x())
{
  /* ... */
}

void Equation_SD_1D::setup_group(const size_t g)
{
  Require(g < d_material->number_groups());
  d_g = g;
}

void Equation_SD_1D::setup_octant(const size_t octant)
{
  Require(octant < 2);
  d_octant = octant;
}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of Equation_SD_1D.cc
//---------------------------------------------------------------------------//



