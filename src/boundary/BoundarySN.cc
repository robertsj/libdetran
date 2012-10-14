//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   BoundarySN.cc
 * \author Jeremy Roberts
 * \date   Mar 25, 2012
 * \brief  Boundary member definitions.
 */
//---------------------------------------------------------------------------//

#include "BoundarySN.hh"
#include "Reflective.hh"
#include "Vacuum.hh"
#include <string>

namespace detran
{

// Constructor.
template<class D>
BoundarySN<D>::BoundarySN(SP_input        input,
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
  names[Mesh::WEST]   = "bc_west";
  names[Mesh::EAST]   = "bc_east";
  names[Mesh::SOUTH]  = "bc_south";
  names[Mesh::NORTH]  = "bc_north";
  names[Mesh::BOTTOM] = "bc_bottom";
  names[Mesh::TOP]    = "bc_top";

  // Assign boundary conditions.
  for(int side = 0; side < 2*D::dimension; side++)
  {
    // Vacuum is default.
    std::string type = "vacuum";
    if (input->check(names[side]))
    {
      // Ugh, this "template" addition is a sneaky thing!
      type = d_input->template get<std::string>(names[side]);
    }
    if (type == "vacuum")
    {
      d_bc[side] = new Vacuum<D>((*this), side, d_input, d_mesh, d_quadrature);
    }
    else if (type == "reflect")
    {
      d_bc[side] = new Reflective<D>((*this), side, d_input, d_mesh, d_quadrature);
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
void BoundarySN<D>::initialize()
{
  THROW("NOT IMPLEMENTED");
}

template <>
void BoundarySN<_3D>::initialize()
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
        d_boundary_flux[Mesh::WEST][g].resize(na,
          boundary_flux_type(nz, vec_dbl(ny, 0.0)));
        d_boundary_flux[Mesh::EAST][g].resize(na,
          boundary_flux_type(nz, vec_dbl(ny, 0.0)));
        // zx planes
        d_boundary_flux[Mesh::SOUTH][g].resize(na,
          boundary_flux_type(nz, vec_dbl(nx, 0.0)));
        d_boundary_flux[Mesh::NORTH][g].resize(na,
          boundary_flux_type(nz, vec_dbl(nx, 0.0)));
        // xy planes
        d_boundary_flux[Mesh::BOTTOM][g].resize(na,
          boundary_flux_type(ny, vec_dbl(nx, 0.0)));
        d_boundary_flux[Mesh::TOP][g].resize(na,
          boundary_flux_type(ny, vec_dbl(nx, 0.0)));
    }
  }
  d_boundary_flux_size[Mesh::WEST]   = na * nz * ny;
  d_boundary_flux_size[Mesh::EAST]   = na * nz * ny;
  d_boundary_flux_size[Mesh::SOUTH]  = na * nz * nx;
  d_boundary_flux_size[Mesh::NORTH]  = na * nz * nx;
  d_boundary_flux_size[Mesh::BOTTOM] = na * ny * nx;
  d_boundary_flux_size[Mesh::TOP]    = na * ny * nx;
}


template <>
void BoundarySN<_2D>::initialize()
{
  int nx = d_mesh->number_cells_x();
  int ny = d_mesh->number_cells_y();
  int na = d_quadrature->number_angles();
  for (int g = 0; g < d_number_groups; g++)
  {
    for (int a = 0; a < d_quadrature->number_angles(); a++)
    {
        // vertical sides
        d_boundary_flux[Mesh::WEST][g].resize(na,
          boundary_flux_type(ny, 0.0));
        d_boundary_flux[Mesh::EAST][g].resize(na,
          boundary_flux_type(ny, 0.0));
        // horizontal sides
        d_boundary_flux[Mesh::SOUTH][g].resize(na,
          boundary_flux_type(nx, 0.0));
        d_boundary_flux[Mesh::NORTH][g].resize(na,
          boundary_flux_type(nx, 0.0));
    }
  }
  d_boundary_flux_size[Mesh::WEST]  = na * ny;
  d_boundary_flux_size[Mesh::EAST]  = na * ny;
  d_boundary_flux_size[Mesh::SOUTH] = na * nx;
  d_boundary_flux_size[Mesh::NORTH] = na * nx;
}

template <>
void BoundarySN<_1D>::initialize()
{
  int na = d_quadrature->number_angles();
  for (int g = 0; g < d_number_groups; g++)
  {
    for (int a = 0; a < d_quadrature->number_angles(); a++)
    {
        // vertical sides
        d_boundary_flux[Mesh::WEST][g].resize(na, 0.0);
        d_boundary_flux[Mesh::EAST][g].resize(na, 0.0);
    }
  }
  d_boundary_flux_size[Mesh::WEST] = na;
  d_boundary_flux_size[Mesh::EAST] = na;
}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of Boundary.cc
//---------------------------------------------------------------------------//
