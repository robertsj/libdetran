//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Sweeper3D.cc
 *  @author robertsj
 *  @date   Nov 7, 2012
 *  @brief  Sweeper3D class definition.
 */
//---------------------------------------------------------------------------//

#include "transport/Sweeper3D.hh"
#include "transport/Equation_DD_3D.hh"
//#include "discretization/Equation_SD_3D.hh"
//#include "discretization/Equation_SC_3D.hh"

namespace detran
{

//---------------------------------------------------------------------------//
template <class EQ>
Sweeper3D<EQ>::Sweeper3D(SP_input input,
                         SP_mesh mesh,
                         SP_material material,
                         SP_quadrature quadrature,
                         SP_state state,
                          BoundaryBase<_3D>::SP_boundary boundary,
                         SP_sweepsource sweepsource)
  : Base(input, mesh, material, quadrature, state, boundary, sweepsource)
  , d_boundary(std::dynamic_pointer_cast<Boundary_T>(boundary))
{
    // Preconditions
    Require(d_boundary);
}

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//---------------------------------------------------------------------------//

TRANSPORT_INSTANTIATE_EXPORT(Sweeper3D<Equation_DD_3D>)
TRANSPORT_TEMPLATE_EXPORT(detran_utilities::SP<Sweeper3D<Equation_DD_3D> >)


} // end namespace detran


