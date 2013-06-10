//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Sweeper2DMOC.cc
 *  @author robertsj
 *  @date   Nov 7, 2012
 *  @brief  Sweeper2DMOC class definition.
 */
//---------------------------------------------------------------------------//

#include "transport/Sweeper2DMOC.hh"
#include "transport/Equation_SC_MOC.hh"

namespace detran
{

//---------------------------------------------------------------------------//
template <class EQ>
Sweeper2DMOC<EQ>::Sweeper2DMOC(SP_input input,
                               SP_mesh mesh,
                               SP_material material,
                               SP_quadrature quadrature,
                               SP_state state,
                               SP_boundary boundary,
                               SP_sweepsource sweepsource)
  : Base(input, mesh, material, quadrature, state, boundary, sweepsource)
  , d_boundary(boundary)
{
    //d_tracks = mesh->tracks();
}

//---------------------------------------------------------------------------//
template <class EQ>
typename Sweeper2DMOC<EQ>::SP_sweeper
Sweeper2DMOC<EQ>::Create(SP_input       input,
                         SP_mesh        mesh,
                         SP_material    material,
                         SP_quadrature  quadrature,
                         SP_state       state,
                         SP_boundary    boundary,
                         SP_sweepsource sweepsource)
{
  SP_sweeper p(new Sweeper2DMOC(input, mesh, material, quadrature,
                                state, boundary, sweepsource));
  return p;
}

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//---------------------------------------------------------------------------//

TRANSPORT_INSTANTIATE_EXPORT(Sweeper2DMOC<Equation_SC_MOC>)
TRANSPORT_TEMPLATE_EXPORT(detran_utilities::SP<Sweeper2DMOC<Equation_SC_MOC> >)

} // end namespace detran
