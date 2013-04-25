//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Sweeper2D.cc
 *  @author robertsj
 *  @date   Nov 7, 2012
 *  @brief  Sweeper2D class definition.
 */
//---------------------------------------------------------------------------//

#include "transport/Sweeper2D.hh"
#include "transport/Equation_DD_2D.hh"
#include "transport/Equation_SD_2D.hh"
#include "transport/Equation_SC_2D.hh"

namespace detran
{

//---------------------------------------------------------------------------//
template <class EQ>
Sweeper2D<EQ>::Sweeper2D(SP_input input,
                         SP_mesh mesh,
                         SP_material material,
                         SP_quadrature quadrature,
                         SP_state state,
                         SP_boundary boundary,
                         SP_sweepsource sweepsource)
  : Base(input, mesh, material, quadrature, state, boundary, sweepsource)
  , d_boundary(boundary)
{
    // Preconditions
    Require(d_boundary);
}

//---------------------------------------------------------------------------//
template <class EQ>
typename Sweeper2D<EQ>::SP_sweeper
Sweeper2D<EQ>::Create(SP_input       input,
                      SP_mesh        mesh,
                      SP_material    material,
                      SP_quadrature  quadrature,
                      SP_state       state,
                      SP_boundary    boundary,
                      SP_sweepsource sweepsource)
{
  SP_sweeper p(new Sweeper2D(input, mesh, material, quadrature,
                             state, boundary, sweepsource));
  return p;
}


//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//---------------------------------------------------------------------------//

TRANSPORT_INSTANTIATE_EXPORT(Sweeper2D<Equation_SD_2D>)
TRANSPORT_INSTANTIATE_EXPORT(Sweeper2D<Equation_DD_2D>)
TRANSPORT_INSTANTIATE_EXPORT(Sweeper2D<Equation_SC_2D>)
TRANSPORT_TEMPLATE_EXPORT(detran_utilities::SP<Sweeper2D<Equation_SD_2D> >)
TRANSPORT_TEMPLATE_EXPORT(detran_utilities::SP<Sweeper2D<Equation_DD_2D> >)
TRANSPORT_TEMPLATE_EXPORT(detran_utilities::SP<Sweeper2D<Equation_SC_2D> >)

} // end namespace detran


