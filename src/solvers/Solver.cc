//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Solver.cc
 *  @author robertsj
 *  @date   Oct 24, 2012
 *  @brief  Solver class definition.
 */
//---------------------------------------------------------------------------//

#include "Solver.hh"

namespace detran
{

//---------------------------------------------------------------------------//
template <class D>
Solver<D>::Solver(SP_state                  state,
                  SP_material               material,
                  SP_boundary               boundary,
                  const vec_externalsource &q_e,
                  SP_fissionsource          q_f)
  : d_state(state)
  , d_material(material)
  , d_boundary(boundary)
  , d_externalsources(q_e)
  , d_fissionsource(q_f)
  , d_maximum_iterations(100)
  , d_tolerance(1e-5)
  , d_print_level(2)
  , d_print_interval(10)
{
  // Preconditions
  Require(d_state);
  // Extract input and mesh from state
  d_input = d_state->get_input();
  d_mesh = d_state->get_mesh();
  Require(d_input);
  Require(d_mesh);
  Require(d_material);
  Require(d_boundary);

  Assert(d_input->check("number_groups"));
  d_number_groups = d_input->get<int>("number_groups");
}

//---------------------------------------------------------------------------//
template <class D>
Solver<D>::~Solver()
{
  /* ... */
}

template class Solver<_1D>;
template class Solver<_2D>;
template class Solver<_3D>;

}


