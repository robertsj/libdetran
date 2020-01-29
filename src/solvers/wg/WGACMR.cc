//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  WGACMR.cc
 *  @brief WGACMR member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "solvers/WGACMR.hh"

namespace detran
{

//----------------------------------------------------------------------------//
template <class D>
WGACMR<D>::WGACMR(SP_state                  state,
                  SP_material               material,
                  SP_boundary               boundary,
                  const vec_externalsource &q_e,
                  SP_fissionsource          q_f,
                  bool                      multiply)
  : Base(state, material, q_e, q_f, multiply)
{
  // Assign the angular order.  \todo Limited to zeroth for testing
  if (d_input->check("inner_acmr_order"))
    d_order = d_input->get<int>("inner_acmr_order");
  Assert(d_order == 0);

  // Create the acceleration mesh
  size_t level = 2;
  CoarseMesh<D> C(d_mesh, level);
  d_coarsemesh = C.get_coarse_mesh();

  // Create tally and give it to the sweeper


}

//----------------------------------------------------------------------------//
template <class D>
WGACMR<D>::~WGACMR()
{
  // TODO Auto-generated destructor stub
}

//----------------------------------------------------------------------------//
template <class D>
void WGACMR<D>::update(const size_t g, moments_type &phi)
{

}

//----------------------------------------------------------------------------//
// IMPLEMENTATION
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
template <class D>
void WGACMR<D>::build_operator()
{

  // Build the Nth order operator

}

//----------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//----------------------------------------------------------------------------//

SOLVERS_INSTANTIATE_EXPORT(WGACMR<_1D>)
SOLVERS_INSTANTIATE_EXPORT(WGACMR<_2D>)
SOLVERS_INSTANTIATE_EXPORT(WGACMR<_3D>)

} // end namespace detran
