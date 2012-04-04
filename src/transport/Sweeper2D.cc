//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Sweeper2D.cc
 * \author Jeremy Roberts
 * \date   Mar 24, 2012
 * \brief  Sweeper2D member definitions.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

// Detran
#include "Sweeper2D.hh"

// System
#include <string>

namespace detran
{

Sweeper2D::Sweeper2D(SP_input input,
                     SP_mesh mesh,
                     SP_material material,
                     SP_quadrature quadrature,
                     SP_state state)
  : SweeperBase(input, mesh, material, quadrature, state)
  , d_psi_h(mesh->number_cells_x(), 0.0)
  , d_psi_v(mesh->number_cells_y(), 0.0)
{

  // Check to see if input has equation type. If not, default.
  std::string equation;
  if (input->check("equation"))
  {
    equation = input->get<std::string>("equation");
  }
  else
  {
    equation = "dd";
  }

  // Set the equation.
  if (equation == "dd")
  {
    Equation_DD_2D::SP_equation poo;
    d_equation = new Equation_DD_2D(mesh, material, quadrature, d_update_psi);
  }
  else
  {
    THROW("Unsupported equation type");
  }


}

} // end namespace detran
