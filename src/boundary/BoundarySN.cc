//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   BoundarySN.cc
 *  @author Jeremy Roberts
 *  @date   Mar 25, 2012
 *  @brief  Boundary member definitions.
 */
//---------------------------------------------------------------------------//

#include "BoundarySN.hh"
#include "Reflective.hh"
#include "Vacuum.hh"
#include "FixedBoundary.hh"
#include "BoundaryTraits.hh"
#include <string>

namespace detran
{

//---------------------------------------------------------------------------//
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
      type = input->template get<std::string>(names[side]);
    if (type == "reflect")
    {
      d_is_reflective[side] = true;
      d_has_reflective = true;
    }
    else
    {
      d_has_vacuum = true;
    }
  }
}

//---------------------------------------------------------------------------//
template<class D>
void BoundarySN<D>::display(bool inout) const
{
  typedef BoundaryValue<D> BV_T;
  // For a given dimension, provide remaining dimensions
  int remdims[3][2] = {{1,2}, {0,2}, {0,1}};
  // Cell indices
  size_t ijk[3] = {0, 0, 0};
  // Loop over all dimensions
  for (int dim0 = 0; dim0 < D::dimension; ++dim0)
  {
    // Bounding cell indices for this dimension
    size_t bound[2] = {0, d_mesh->number_cells(dim0)-1};
    // Other dimensions
    size_t dim1 = remdims[dim0][0];
    size_t dim2 = remdims[dim0][1];
    // Loop over directions - and +
    for (size_t dir = 0; dir < 2; ++dir)
    {
      // Surface index
      size_t surface = 2 * dim0 + dir;
      std::cout << " SURFACE = " << surface << std::endl;
      // Index and width along this direction
      ijk[dim0] = bound[dir];
      for (size_t g = 0; g < d_number_groups; g++)
      {
        std::cout << "   GROUP = " << g << std::endl;
        for (int oo = 0; oo < std::pow((float)2, (int)(D::dimension-1)); ++oo)
        {
          int o = d_quadrature->outgoing_octant(surface)[oo];
          if (inout)
            o = d_quadrature->incident_octant(surface)[oo];
          std::cout << "     OCTANT = " << o << std::endl;
          for (size_t a = 0; a < d_quadrature->number_angles_octant(); ++a)
          {
            std::cout << "       ANGLE = " << a << std::endl;
            for (ijk[dim1] = 0; ijk[dim1] < d_mesh->number_cells(dim1); ++ijk[dim1])
            {
              //std::cout << "dim1=" << dim1 << " val=" << ijk[dim1] << std::endl;
              if (d_mesh->dimension() == 3)
              {
                for (ijk[dim2] = 0; ijk[dim2] < d_mesh->number_cells(dim2); ++ijk[dim2])
                {
                  std::cout << " dim2="<< dim2 << " val=" << ijk[dim2] << std::endl;
                  std::cout << BV_T::value((*this)(surface, o, a, g), ijk[dim1], ijk[dim2]) << " " << std::endl;
                }
                std::cout << std::endl;
              }
              else
              {
                std::cout << BV_T::value((*this)(surface, o, a, g), ijk[dim1], 0) << " ";
              }
            }
            std::cout << std::endl;
          }
        }
      }
    }
  }
}

//---------------------------------------------------------------------------//
// IMPLEMENTATION
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
template <class D>
void BoundarySN<D>::initialize()
{
  THROW("NOT IMPLEMENTED");
}

//---------------------------------------------------------------------------//
template <>
void BoundarySN<_3D>::initialize()
{
  size_t nx = d_mesh->number_cells_x();
  size_t ny = d_mesh->number_cells_y();
  size_t nz = d_mesh->number_cells_z();
  size_t na = d_quadrature->number_angles();
  for (size_t g = 0; g < d_number_groups; ++g)
  {
    for (size_t a = 0; a < d_quadrature->number_angles(); ++a)
    {
      // yz planes
      d_boundary_flux[Mesh::WEST][g].resize(na, bf_type(nz, vec_dbl(ny, 0.0)));
      d_boundary_flux[Mesh::EAST][g].resize(na, bf_type(nz, vec_dbl(ny, 0.0)));
      // zx planes
      d_boundary_flux[Mesh::SOUTH][g].resize(na, bf_type(nz, vec_dbl(nx, 0.0)));
      d_boundary_flux[Mesh::NORTH][g].resize(na, bf_type(nz, vec_dbl(nx, 0.0)));
      // xy planes
      d_boundary_flux[Mesh::BOTTOM][g].resize(na, bf_type(ny, vec_dbl(nx, 0.0)));
      d_boundary_flux[Mesh::TOP][g].resize(na, bf_type(ny, vec_dbl(nx, 0.0)));
    }
  }
  d_boundary_flux_size[Mesh::WEST]   = na * nz * ny;
  d_boundary_flux_size[Mesh::EAST]   = na * nz * ny;
  d_boundary_flux_size[Mesh::SOUTH]  = na * nz * nx;
  d_boundary_flux_size[Mesh::NORTH]  = na * nz * nx;
  d_boundary_flux_size[Mesh::BOTTOM] = na * ny * nx;
  d_boundary_flux_size[Mesh::TOP]    = na * ny * nx;
}

//---------------------------------------------------------------------------//
template <>
void BoundarySN<_2D>::initialize()
{
  size_t nx = d_mesh->number_cells_x();
  size_t ny = d_mesh->number_cells_y();
  size_t na = d_quadrature->number_angles();
  for (size_t g = 0; g < d_number_groups; ++g)
  {
    for (size_t a = 0; a < d_quadrature->number_angles(); ++a)
    {
      // vertical sides
      d_boundary_flux[Mesh::WEST][g].resize(na, bf_type(ny, 0.0));
      d_boundary_flux[Mesh::EAST][g].resize(na, bf_type(ny, 0.0));
      // horizontal sides
      d_boundary_flux[Mesh::SOUTH][g].resize(na, bf_type(nx, 0.0));
      d_boundary_flux[Mesh::NORTH][g].resize(na, bf_type(nx, 0.0));
    }
  }
  d_boundary_flux_size[Mesh::WEST]  = na * ny;
  d_boundary_flux_size[Mesh::EAST]  = na * ny;
  d_boundary_flux_size[Mesh::SOUTH] = na * nx;
  d_boundary_flux_size[Mesh::NORTH] = na * nx;
}

//---------------------------------------------------------------------------//
template <>
void BoundarySN<_1D>::initialize()
{
  size_t na = d_quadrature->number_angles();
  for (size_t g = 0; g < d_number_groups; ++g)
  {
    for (size_t a = 0; a < d_quadrature->number_angles(); ++a)
    {
      // vertical sides
      d_boundary_flux[Mesh::WEST][g].resize(na, 0.0);
      d_boundary_flux[Mesh::EAST][g].resize(na, 0.0);
    }
  }
  d_boundary_flux_size[Mesh::WEST] = na;
  d_boundary_flux_size[Mesh::EAST] = na;
}

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//---------------------------------------------------------------------------//

BOUNDARY_INSTANTIATE_EXPORT(BoundarySN<_1D>)
BOUNDARY_INSTANTIATE_EXPORT(BoundarySN<_2D>)
BOUNDARY_INSTANTIATE_EXPORT(BoundarySN<_3D>)
BOUNDARY_TEMPLATE_EXPORT(detran_utilities::SP<BoundarySN<_1D> >)
BOUNDARY_TEMPLATE_EXPORT(detran_utilities::SP<BoundarySN<_2D> >)
BOUNDARY_TEMPLATE_EXPORT(detran_utilities::SP<BoundarySN<_3D> >)

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of BoundarySN.cc
//---------------------------------------------------------------------------//
