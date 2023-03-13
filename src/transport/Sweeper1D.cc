//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Sweeper1D.cc
 *  @author robertsj
 *  @date   Nov 7, 2012
 *  @brief  Sweeper1D class definition.
 */
//---------------------------------------------------------------------------//

#include "Sweeper1D.hh"
#include "Equation_DD_1D.hh"
#include "Equation_SD_1D.hh"
#include "Equation_SC_1D.hh"

namespace detran
{

//---------------------------------------------------------------------------//
template <class EQ>
Sweeper1D<EQ>::Sweeper1D(SP_input input,
                         SP_mesh mesh,
                         SP_material material,
                         SP_quadrature quadrature,
                         SP_state state,
                         BoundaryBase<_1D>::SP_boundary boundary,
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

TRANSPORT_INSTANTIATE_EXPORT(Sweeper1D<Equation_SD_1D>)
TRANSPORT_INSTANTIATE_EXPORT(Sweeper1D<Equation_DD_1D>)
TRANSPORT_INSTANTIATE_EXPORT(Sweeper1D<Equation_SC_1D>)
TRANSPORT_TEMPLATE_EXPORT(detran_utilities::SP<Sweeper1D<Equation_SD_1D> >)
TRANSPORT_TEMPLATE_EXPORT(detran_utilities::SP<Sweeper1D<Equation_DD_1D> >)
TRANSPORT_TEMPLATE_EXPORT(detran_utilities::SP<Sweeper1D<Equation_SC_1D> >)

} // end namespace detran
