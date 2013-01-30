//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   BoundaryFactory.t.hh
 *  @author robertsj
 *  @date   Jan 30, 2013
 *  @brief  BoundaryFactory.t class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_BOUNDARYFACTORY_T_HH_
#define detran_BOUNDARYFACTORY_T_HH_

#include "BoundaryDiffusion.hh"
#include "BoundarySN.hh"
#include "BoundaryMOC.hh"

namespace detran
{

//---------------------------------------------------------------------------//
// DIFFUSION
//---------------------------------------------------------------------------//

template <class D>
typename BoundaryFactory<D, BoundaryDiffusion>::SP_boundary
build(detran_utilities::InputDB::SP_input       input,
      detran_geometry::Mesh::SP_mesh            mesh,
      detran_angle::Quadrature::SP_quadrature   quad)
{
  typename BoundaryDiffusion<D>::SP_boundary
    b(new BoundaryDiffusion<D>(input, mesh));
  return b;
}

//---------------------------------------------------------------------------//
// SN
//---------------------------------------------------------------------------//

template <class D>
typename BoundaryFactory<D, BoundarySN>::SP_boundary
build(detran_utilities::InputDB::SP_input       input,
      detran_geometry::Mesh::SP_mesh            mesh,
      detran_angle::Quadrature::SP_quadrature   quad)
{
  typename BoundarySN<D>::SP_boundary
    b(new BoundarySN<D>(input, mesh, quad));
  return b;
}

//---------------------------------------------------------------------------//
// MOC
//---------------------------------------------------------------------------//

template <class D>
typename BoundaryFactory<D, BoundaryMOC>::SP_boundary
build(detran_utilities::InputDB::SP_input       input,
      detran_geometry::Mesh::SP_mesh            mesh,
      detran_angle::Quadrature::SP_quadrature   quad)
{
  typename BoundaryMOC<D>::SP_boundary
    b(new BoundaryMOC<D>(input, mesh));
  return b;
}

} // end namespace detran

#endif /* detran_BOUNDARYFACTORY_T_HH_ */
