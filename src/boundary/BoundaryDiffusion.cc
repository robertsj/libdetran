//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   BoundaryDiffusion.cc
 *  @author robertsj
 *  @date   Sep 11, 2012
 *  @brief  BoundaryDiffusion class definition.
 */
//---------------------------------------------------------------------------//

#include "BoundaryDiffusion.hh"

namespace detran
{

//---------------------------------------------------------------------------/
template <class D>
BoundaryDiffusion<D>::BoundaryDiffusion(SP_input        input,
                                        SP_mesh         mesh)
  : Base(input, mesh)
  , d_boundary_flux(2, vec2_boundary_flux(2*D::dimension))
{
  initialize();
}

//---------------------------------------------------------------------------/
template <class D>
typename BoundaryDiffusion<D>::SP_base
BoundaryDiffusion<D>::Create(SP_input         input,
                             SP_mesh          mesh)
{
  SP_boundary p(new BoundaryDiffusion(input, mesh));
  return p;
}

//---------------------------------------------------------------------------//
// Implementation
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------/
template <class D>
void BoundaryDiffusion<D>::initialize()
{
  D::not_implemented();
}

//---------------------------------------------------------------------------/
template <>
void BoundaryDiffusion<_3D>::initialize()
{
  int nx = d_mesh->number_cells_x();
  int ny = d_mesh->number_cells_y();
  int nz = d_mesh->number_cells_z();
  int ng = d_number_groups;
  for (int inout = 0; inout < 2; inout++)
  {
    // yz planes
    d_boundary_flux[inout][Mesh::WEST].resize(
      ng, bf_type(nz, vec_dbl(ny, 0.0)));
    d_boundary_flux[inout][Mesh::EAST].resize(
      ng, bf_type(nz, vec_dbl(ny, 0.0)));
    // xz planes
    d_boundary_flux[inout][Mesh::SOUTH].resize(
      ng, bf_type(nz, vec_dbl(nx, 0.0)));
    d_boundary_flux[inout][Mesh::NORTH].resize(
      ng, bf_type(nz, vec_dbl(nx, 0.0)));
    // xy planes
    d_boundary_flux[inout][Mesh::BOTTOM].resize(
      ng, bf_type(ny, vec_dbl(nx, 0.0)));
    d_boundary_flux[inout][Mesh::TOP].resize(
      ng, bf_type(ny, vec_dbl(nx, 0.0)));
  }
  d_boundary_flux_size[Mesh::WEST]   = ny * nz;
  d_boundary_flux_size[Mesh::EAST]   = ny * nz;
  d_boundary_flux_size[Mesh::SOUTH]  = nx * nz;
  d_boundary_flux_size[Mesh::NORTH]  = nx * nz;
  d_boundary_flux_size[Mesh::BOTTOM] = nx * ny;
  d_boundary_flux_size[Mesh::TOP]    = nx * ny;
}

//---------------------------------------------------------------------------/
template <>
void BoundaryDiffusion<_2D>::initialize()
{
  int nx = d_mesh->number_cells_x();
  int ny = d_mesh->number_cells_y();
  int ng = d_number_groups;
  for (int inout = 0; inout < 2; inout++)
  {
    // vertical sides
    d_boundary_flux[inout][Mesh::WEST].resize(ng, bf_type(ny, 0.0));
    d_boundary_flux[inout][Mesh::EAST].resize(ng, bf_type(ny, 0.0));
    // horizontal sides
    d_boundary_flux[inout][Mesh::SOUTH].resize(ng, bf_type(nx, 0.0));
    d_boundary_flux[inout][Mesh::NORTH].resize(ng, bf_type(nx, 0.0));
  }
  d_boundary_flux_size[Mesh::WEST]  = ny;
  d_boundary_flux_size[Mesh::EAST]  = ny;
  d_boundary_flux_size[Mesh::SOUTH] = nx;
  d_boundary_flux_size[Mesh::NORTH] = nx;
}

//---------------------------------------------------------------------------//
template <>
void BoundaryDiffusion<_1D>::initialize()
{
  int ng = d_number_groups;
  for (int inout = 0; inout < 2; inout++)
  {
    // vertical sides
    d_boundary_flux[inout][Mesh::WEST].resize(ng, 0.0);
    d_boundary_flux[inout][Mesh::EAST].resize(ng, 0.0);
  }
  d_boundary_flux_size[Mesh::WEST]  = 1;
  d_boundary_flux_size[Mesh::EAST]  = 1;
}

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//---------------------------------------------------------------------------//

BOUNDARY_INSTANTIATE_EXPORT(BoundaryDiffusion<_1D>)
BOUNDARY_INSTANTIATE_EXPORT(BoundaryDiffusion<_2D>)
BOUNDARY_INSTANTIATE_EXPORT(BoundaryDiffusion<_3D>)
BOUNDARY_TEMPLATE_EXPORT(detran_utilities::SP<BoundaryDiffusion<_1D> >)
BOUNDARY_TEMPLATE_EXPORT(detran_utilities::SP<BoundaryDiffusion<_2D> >)
BOUNDARY_TEMPLATE_EXPORT(detran_utilities::SP<BoundaryDiffusion<_3D> >)

} // end namespace detran


