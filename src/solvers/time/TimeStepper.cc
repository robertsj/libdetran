//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   TimeStepper.cc
 *  @brief  TimeStepper
 *  @author Jeremy Roberts
 *  @date   Nov 16, 2012
 */
//---------------------------------------------------------------------------//

#ifndef detran_TIMESTEPPER_CC_
#define detran_TIMESTEPPER_CC_

#include "TimeStepper.hh"
#include "kinetics/SyntheticDiscreteSource.hh"
#include "kinetics/SyntheticMomentSource.hh"
#include <cmath>

namespace detran
{

//---------------------------------------------------------------------------//
template <class D>
TimeStepper<D>::TimeStepper(SP_input       input,
                            SP_material    material,
                            SP_mesh        mesh,
                            bool           multiply)
  : d_input(input)
  , d_material(material)
  , d_mesh(mesh)
  , d_multiply(multiply)
  , d_discrete(false)
  , d_number_groups(0)
  , d_dt(1.0)
  , d_final_time(10.0)
  , d_number_steps(10)
  , d_scheme(BDF1)
  , d_do_output(false)
{
  // Preconditions
  Require(d_input);
  Require(d_material);
  Require(d_mesh);

  // Create the fixed source solver and extract the quadrature.
  d_solver = new Fixed_T(d_input, d_material, d_mesh, d_multiply);
  d_quadrature = d_solver->quadrature();

  // Set the number of energy groups
  d_number_groups = d_material->number_groups();

  // Is this a moment or discrete problem?
  if (d_input->check("ts_discrete"))
  {
    d_discrete = d_input->template get<int>("time_discrete");
  }

  // Create synthetic source
  if (d_discrete)
  {
    Insist(d_quadrature,
           "Can't do time stepping using psi without a quadrature");
    d_syntheticsource = new SyntheticDiscreteSource(d_number_groups,
                                                    d_mesh,
                                                    d_quadrature,
                                                    d_material);
  }
  else
  {
    d_syntheticsource = new SyntheticMomentSource(d_number_groups,
                                                  d_mesh,
                                                  d_material);
  }

  // Add the synthetic source to the source vector.  This is the
  // *always* the first source.
  add_source(d_syntheticsource);

  // Get the maximum time and the time step to use.
  if (d_input->check("ts_final_time"))
    d_final_time = d_input->template get<double>("ts_final_time");
  Assert(d_final_time > 0.0);
  if (d_input->check("ts_step_size"))
    d_dt = d_input->template get<double>("ts_step_size");
  Assert(d_dt > 0.0);
  // Compute the number of steps.  This may lead to a slightly longer
  // final time.  There are tons of options we can put in later, but
  // a fixed time step will suffice for now.
  d_number_steps = std::ceil(d_final_time / d_dt);

  // Get the integration scheme.
  if (d_input->check("ts_scheme"))
    d_scheme = d_input->template get<int>("ts_scheme");
  Require(d_scheme < END_TIME_SCHEMES);



  // Check whether to include precursors.  If the material
  // has a nonzero beta, and fission is included, then
  // the precursors have to be included to balance



  // Size and initialize the state and precursor vectors for
  // previous time steps.
  size_t number_terms = d_scheme;
  if (d_scheme == 0) number_terms++;
  d_states.resize(number_terms);
  if (d_multiply) d_precursors.resize(number_terms);
  for (int i = 0; i < number_terms; ++i)
  {
    d_states[i] = new State(d_input, d_mesh, d_quadrature);
    if (d_multiply)
    {
      d_precursors[i] = new Precursors(d_mesh->number_cells(),
                                       d_material->number_precursor_groups());
    }
  }

  if (d_do_output) d_silooutput = new detran_ioutils::SiloOutput(d_mesh);
}

//---------------------------------------------------------------------------//
template <class D>
void TimeStepper<D>::add_source(SP_externalsource source)
{
  // Preconditions
  Require(source);

  d_sources.push_back(source);
}

//---------------------------------------------------------------------------//
template <class D>
void TimeStepper<D>::solve(SP_state initial_state)
{
  // Preconditions
  Require(initial_state);

  // Set the state and inititalize the precursors if necessary.
  d_state = initial_state;
  initialize_precursors();

  // Output the initial state
  if (d_do_output) d_silooutput->write_time_flux(0, d_state, d_discrete);

  // Fill previous states with backward Euler steps
  // if required.

  // Perform time steps
  double t = 0.0;
  for (int i = 0; i < d_number_steps; ++i)
  {
    // Increment the time
    t += d_dt;

    double factor = 1.0;
    if (d_scheme == IMP) factor = 0.5;

    // Update the material
    d_material->update(t, factor * d_dt);

    // Update the t
    update_external(t);

    // Update the synthetic source
    update_synthetic(t, factor * d_dt);

    // Update the multigroup solver
    d_solver->update();

    // Take the step and iterate if necessary.
    solve(factor * d_dt);
    // >>> Insert iteration <<<

    SP_state current_state = d_solver->state();

    // Output the initial state
    if (d_do_output) d_silooutput->write_time_flux(i+1, d_state, d_discrete);

    // Cycle the states and precursors.

  }


}

//---------------------------------------------------------------------------//
template <class D>
void TimeStepper<D>::solve(const double dt)
{
  // Solve the MG problem
  d_solver->solve();

  // Cycle the states and precursors.
  cycle_states_precursors();

  // Get the current state and copy to first element.
  d_state = d_solver->state();
  *d_states[0] = *d_state;

  // Compute the precursors

}

//---------------------------------------------------------------------------//
template <class D>
void TimeStepper<D>::initialize_precursors()
{
  // Preconditions
  Require(d_state);

  // Skip if we have no multiplication
  if (!d_multiply) return;
  Assert(d_precursors.size() > 0);

  /*
   *  Given the flux, phi_g, the precursors at steady state
   *  are defined by
   *    lambda_i * C_i = beta_i * sum_g^G nsf_g * phi_g
   *  The easiest way to accomplish this is with the
   *  fission source
   */
  d_fissionsource->update();
  const State::moments_type &fd = d_fissionsource->density();

  const vec_int &mt = d_mesh->mesh_map("MATERIAL");

  for (int i = 0; i < d_material->number_precursor_groups(); ++i)
  {
    double inv_lambda = 1.0 / d_material->lambda(i);
    for (int cell = 0; cell < d_mesh->number_cells(); ++cell)
    {
      d_precursors[0]->C(i)[cell] =
        inv_lambda * d_material->beta(mt[cell], i) * fd[cell];
    }
  }

}


//---------------------------------------------------------------------------//
template <class D>
void TimeStepper<D>::cycle_states_precursors()
{
  SP_state tmp_state;
  SP_precursors tmp_precursors;

  /*
   *  Consider the state vector
   *    [state1 state2 state3]
   *  where 1, 2, and 3 represent the previous iterates starting
   *  with the most recent.  We cycle as follows:
   *    states -> [state3 state1 state2]
   *  using pointers.  We do likewise for precursors.  Then,
   *  the current state and precursors are copied into the
   *  first element.
   *
   */

  // Save the first element.
  tmp_state = d_states[0];
  if (d_precursors.size())tmp_precursors = d_precursors[0];

  for (int i = 0; i < d_states.size() - 1; ++i)
  {
    d_states[i] = d_states[i+1];
    if (d_precursors.size()) d_precursors[i] = d_precursors[i+1];
  }
  d_states[d_number_previous - 1] = tmp_state;
  if (d_precursors.size()) d_precursors[d_number_previous-1] = tmp_precursors;

}

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//---------------------------------------------------------------------------//

template class TimeStepper<_1D>;
template class TimeStepper<_2D>;
template class TimeStepper<_3D>;

} // end namespace detran

#endif // TIMESTEPPER_CC_ 

//---------------------------------------------------------------------------//
//              end of file TimeStepper.cc
//---------------------------------------------------------------------------//
