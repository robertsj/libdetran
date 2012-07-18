//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Boundary.cc
 * \author Jeremy Roberts
 * \date   Mar 25, 2012
 * \brief  Boundary member definitions.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// Detran
#include "Boundary.hh"
#include "Reflective.hh"
#include "Vacuum.hh"

// System
#include <string>

namespace detran
{

// Constructor.
template<class D>
Boundary<D>::Boundary(SP_input        input,
                      SP_mesh         mesh,
                      SP_quadrature   quadrature)
  : Base(input, mesh)
  , d_quadrature(quadrature)
  , d_boundary_flux(2*D::dimension, // number of sides
                    vec2_boundary_flux(d_number_groups))
  , d_bc(2*D::dimension)
{
  Require(d_quadrature);

  // Allocate the boundary flux container.
  initialize();

  std::vector<std::string> names(6);
  names[Mesh::LEFT]   = "bc_left";
  names[Mesh::RIGHT]  = "bc_right";
  names[Mesh::BOTTOM] = "bc_bottom";
  names[Mesh::TOP]    = "bc_top";
  names[Mesh::SOUTH]  = "bc_south";
  names[Mesh::NORTH]  = "bc_north";

  // Create SP of me.

  // Assign boundary conditions.
  for(int side = 0; side < 2*D::dimension; side++)
  {
    // Vacuum is default.
    std::string type = "vacuum";
    if (input->check(names[side]))
    {
      type = input->get<std::string>(names[side]);
    }
    if (type == "vacuum")
    {
      d_bc[side] = new Vacuum<D>((*this), side, input, mesh, quadrature);
    }
    else if (type == "reflect")
    {
      d_bc[side] = new Reflective<D>((*this), side, input, mesh, quadrature);
      d_is_reflective[side] = true;
      d_has_reflective = true;
    }
    else
    {
      type.append(" is not a supported bc type.");
      THROW(type);
      break;
    }
  }

}

//---------------------------------------------------------------------------//
// Implementation
//---------------------------------------------------------------------------//

template <class D>
void Boundary<D>::initialize()
{
  int nx = d_mesh->number_cells_x();
  int ny = d_mesh->number_cells_y();
  int nz = d_mesh->number_cells_z();
  int na = d_quadrature->number_angles();
  for (int g = 0; g < d_number_groups; g++)
  {
    for (int a = 0; a < d_quadrature->number_angles(); a++)
    {
        // yz planes
        d_boundary_flux[Mesh::LEFT][g].resize(na,
          boundary_flux_type(nz, vec_dbl(ny, 0.0)));
        d_boundary_flux[Mesh::RIGHT][g].resize(na,
          boundary_flux_type(nz, vec_dbl(ny, 0.0)));
        // xz planes
        d_boundary_flux[Mesh::BOTTOM][g].resize(na,
          boundary_flux_type(nz, vec_dbl(nx, 0.0)));
        d_boundary_flux[Mesh::TOP][g].resize(na,
          boundary_flux_type(nz, vec_dbl(nx, 0.0)));
        // xy planes
        d_boundary_flux[Mesh::SOUTH][g].resize(na,
          boundary_flux_type(ny, vec_dbl(nx, 0.0)));
        d_boundary_flux[Mesh::NORTH][g].resize(na,
          boundary_flux_type(ny, vec_dbl(nx, 0.0)));
    }
  }
  d_boundary_flux_size[Mesh::LEFT]   = na * ny * nz;
  d_boundary_flux_size[Mesh::RIGHT]  = na * ny * nz;
  d_boundary_flux_size[Mesh::BOTTOM] = na * nx * nz;
  d_boundary_flux_size[Mesh::TOP]    = na * nx * nz;
  d_boundary_flux_size[Mesh::SOUTH]  = na * nx * ny;
  d_boundary_flux_size[Mesh::NORTH]  = na * nx * ny;
}

template <>
void Boundary<_2D>::initialize()
{
  int nx = d_mesh->number_cells_x();
  int ny = d_mesh->number_cells_y();
  int na = d_quadrature->number_angles();
  for (int g = 0; g < d_number_groups; g++)
  {
    for (int a = 0; a < d_quadrature->number_angles(); a++)
    {
        // vertical sides
        d_boundary_flux[Mesh::LEFT][g].resize(na,
          boundary_flux_type(ny, 0.0));
        d_boundary_flux[Mesh::RIGHT][g].resize(na,
          boundary_flux_type(ny, 0.0));
        // horizontal sides
        d_boundary_flux[Mesh::BOTTOM][g].resize(na,
          boundary_flux_type(nx, 0.0));
        d_boundary_flux[Mesh::TOP][g].resize(na,
          boundary_flux_type(nx, 0.0));
    }
  }
  d_boundary_flux_size[Mesh::LEFT]   = na * ny;
  d_boundary_flux_size[Mesh::RIGHT]  = na * ny;
  d_boundary_flux_size[Mesh::BOTTOM] = na * nx;
  d_boundary_flux_size[Mesh::TOP]    = na * nx;
}

template <>
void Boundary<_1D>::initialize()
{
  int na = d_quadrature->number_angles();
  for (int g = 0; g < d_number_groups; g++)
  {
    for (int a = 0; a < d_quadrature->number_angles(); a++)
    {
        // vertical sides
        d_boundary_flux[Mesh::LEFT][g].resize(na, 0.0);
        d_boundary_flux[Mesh::RIGHT][g].resize(na, 0.0);
    }
  }
  d_boundary_flux_size[Mesh::LEFT]   = na;
  d_boundary_flux_size[Mesh::RIGHT]  = na;
}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of Boundary.cc
//---------------------------------------------------------------------------//
