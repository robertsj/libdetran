//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Boundary.cc
 * \author Jeremy Roberts
 * \date   Mar 25, 2012
 * \brief  Boundary member definitions.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#include "Boundary.hh"

namespace detran
{

// Constructor.
template<class D>
Boundary<D>::Boundary(SP_input        input,
                      SP_mesh         mesh,
                      SP_quadrature   quadrature)
  : d_input(input)
  , d_mesh(mesh)
  , d_quadrature(quadrature)
  , d_number_groups(input->get<int>("number_groups"))
  , d_boundary_flux(2*D::dimension, // number of sides
                    vec2_boundary_flux(d_number_groups))
{
  Require(input);
  Require(mesh);
  Require(quadrature);
  Require(d_number_groups > 0);

  // Allocate the boundary flux container.
  initialize();

}

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
}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of Boundary.cc
//---------------------------------------------------------------------------//
