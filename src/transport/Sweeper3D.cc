//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Sweeper3D.cc
 *  @author robertsj
 *  @date   Nov 7, 2012
 *  @brief  Sweeper3D class definition.
 */
//---------------------------------------------------------------------------//

#include "Sweeper3D.hh"
#include "discretization/Equation_DD_3D.hh"
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
typename Sweeper3D<EQ>::SP_sweeper
Sweeper3D<EQ>::Create(SP_input       input,
                      SP_mesh        mesh,
                      SP_material    material,
                      SP_quadrature  quadrature,
                      SP_state       state,
                      SP_boundary    boundary,
                      SP_sweepsource sweepsource)
{
  SP_sweeper p(new Sweeper3D(input, mesh, material, quadrature,
                             state, boundary, sweepsource));
  return p;
}

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//---------------------------------------------------------------------------//

template class Sweeper3D<Equation_DD_3D>;
//template class Sweeper3D<Equation_SD_3D>;
//template class Sweeper3D<Equation_SC_3D>;

} // end namespace detran


