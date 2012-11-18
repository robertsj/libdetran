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
  , d_step_factor(1.0)
  , d_final_time(10.0)
  , d_number_steps(10)
  , d_scheme(BDF1)
  , d_do_output(false)
  , d_fixup(false)
{
  // Preconditions
  Require(d_input);
  Require(d_material);
  Require(d_mesh);

  //-------------------------------------------------------------------------//
  // SETUP FIXED SOLVER
  //-------------------------------------------------------------------------//

  d_solver = new Fixed_T(d_input, d_material, d_mesh, d_multiply);
  d_solver->setup();

  // Extract the quadrature
  d_quadrature = d_solver->quadrature();
  Ensure(d_quadrature);

  // Set the material state.
  d_material->set_state(d_solver->state());

  // Set the number of energy groups
  d_number_groups = d_material->number_groups();

  //-------------------------------------------------------------------------//
  // SETUP SYNTHETIC SOURCE
  //-------------------------------------------------------------------------//

  if (d_input->check("ts_discrete"))
    d_discrete = d_input->template get<int>("ts_discrete");

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
  d_solver->set_source(d_syntheticsource);

  //-------------------------------------------------------------------------//
  // SETUP TIME DISCRETIZATION
  //-------------------------------------------------------------------------//

  // Get the maximum time and the time step to use.
  if (d_input->check("ts_final_time"))
    d_final_time = d_input->template get<double>("ts_final_time");
  Assert(d_final_time > 0.0);
  if (d_input->check("ts_step_size"))
    d_dt = d_input->template get<double>("ts_step_size");
  Assert(d_dt > 0.0);

  // Compute the number of steps.  May result in longer time than requested!
  d_number_steps = std::ceil(d_final_time / d_dt);

  if (d_input->check("ts_max_steps"))
  {
   int max_steps = d_input->template get<int>("ts_max_steps");
   Assert(max_steps > 1);
   if (d_number_steps > max_steps) d_number_steps = max_steps;
  }

  // Get the integration scheme and set the order.
  if (d_input->check("ts_scheme"))
    d_scheme = d_input->template get<int>("ts_scheme");
  Require(d_scheme < END_TIME_SCHEMES);
  d_order = d_scheme;
  if (d_scheme == 0) d_order++;

  // If implicit midpoint, half the step and set step factor
  if (d_scheme == IMP)
  {
    d_dt *= 0.5;
    d_step_factor = 2.0;
  }

  //-------------------------------------------------------------------------//
  // SETUP STATE AND PRECURSOR VECTORS
  //-------------------------------------------------------------------------//

  d_states.resize(d_order);
  if (d_multiply) d_precursors.resize(d_order);
  for (int i = 0; i < d_order; ++i)
  {
    d_states[i] = new State(d_input, d_mesh, d_quadrature);
    if (d_multiply)
    {
      d_precursors[i] = new Precursors(d_mesh->number_cells(),
                                       d_material->number_precursor_groups());
    }
  }
  if (d_multiply)
  {
    d_precursor = new Precursors(d_mesh->number_cells(),
                                 d_material->number_precursor_groups());
  }

  //-------------------------------------------------------------------------//
  // SETUP OUTPUT
  //-------------------------------------------------------------------------//

  if (d_input->check("ts_output"))
  {
    d_do_output = d_input->template get<int>("ts_output");
    if (d_do_output) d_silooutput = new detran_ioutils::SiloOutput(d_mesh);
  }

}

//---------------------------------------------------------------------------//
template <class D>
void TimeStepper<D>::add_source(SP_tdsource source)
{
  // Preconditions
  Require(source);

  // Add the source to my vector and add it to the fixed solver.  We keep
  // our own copy to update the time-dependence.
  d_sources.push_back(source);
  d_solver->set_source(source);
}

//---------------------------------------------------------------------------//
template <class D>
void TimeStepper<D>::solve(SP_state initial_state)
{
  // Preconditions
  Require(initial_state);

  // Set the state and initialize the precursors if necessary.
  d_state = initial_state;
  d_state->display();
  initialize_precursors();
  *d_states[0] = *d_state;

  // Output the initial state
  if (d_do_output) d_silooutput->write_time_flux(0, d_state, true);

  // Set the solver
  d_solver->set_solver();

  // Fill previous states with backward Euler steps if required.
  std::cout << " phi=" << d_state->phi(0)[0] << std::endl;
  // Perform time steps
  double t = 0.0;
  for (int i = 0; i < d_number_steps; ++i)
  {
    t += d_step_factor * d_dt;

    std::cout << " time =  " << t << std::endl;

    size_t iteration = 0;
    for (; iteration < 1; ++iteration)
    {
      // Update the material, sources, and solver

      d_material->update(t, d_dt);
      d_material->display();
      update_sources(t, d_dt);
      d_solver->update();

      // Solve the MG problem for the new state and update the precursors
      d_solver->solve();
      d_state = d_solver->state();
      //d_state->display();
      std::cout << " phi=" << d_state->phi(0)[0] << std::endl;
      update_precursors();
      if (d_scheme == IMP) extrapolate();

      // Cycle the previous iterates and copy the current solution
      cycle_states_precursors();
      *d_states[0] = *d_state;
      if (d_multiply) *d_precursors[0] = *d_precursor;

      // >>> check convergence <<<

    } // end iterations

    // Output the initial state
    Assert(d_do_output);
    if (d_do_output) d_silooutput->write_time_flux(i+1, d_state, true);

  } // end time steps

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
void TimeStepper<D>::update_precursors()
{
  // Skip if we have no multiplication
  if (!d_multiply) return;

  // Update the fission density.
  d_fissionsource->update();

  /*
   *  The precursors are defined via
   *    f(t, C_i(t)) =  dC/dt = -lambda_i*C_i + beta_i * sum_g X_ig F_g phi_g
   *  For the BDF and IMP schemes, we have
   *   (1/Delta) sum_{j=0}^m a_j * C_i(n+j+1) = -lambda_i*C_i(n+1) + beta_i * sum_g X_ig F_g phi_g
   *   (a_j/Delta + lambda)*C_i(n+1) = -(1/Delta) sum_{j=1}^{m} a_J * C_i(n+j+1) + beta_i * sum_g X_ig F_g phi_g
   */

  const detran_utilities::vec_int &mt = d_mesh->mesh_map("MATERIAL");

  for (size_t i = 0; i < d_material->number_precursor_groups(); ++i)
  {
    double lambda = d_material->lambda(i);
    double A = 1.0 / (bdf_coefs[d_order-1][0] / d_dt + lambda);

    for (int cell = 0; cell < d_mesh->number_cells(); ++cell)
    {
      // Add flux contribution
      double value = 0.0;
      for (size_t g = 0; g < d_material->number_groups(); ++g)
      {
        value += A * d_material->beta(mt[cell], i) *
                 d_material->chi_d(mt[cell], i, g);
      }
      // Add previous time step contributions
      for (size_t j = 1; j <= d_order; ++j)
      {
        value += (-A/d_dt) * bdf_coefs[d_order-1][j] * d_precursors[j-1]->C(i)[cell];
      }
      d_precursor->C(i)[cell] = value;
    }
  }
}

//---------------------------------------------------------------------------//
template <class D>
void TimeStepper<D>::cycle_states_precursors()
{
  // Preconditions
  Require(d_states.size() == d_order);

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

  for (int i = 0; i < d_order - 1; ++i)
  {
    d_states[i] = d_states[i + 1];
    if (d_precursors.size()) d_precursors[i] = d_precursors[i + 1];
  }
  d_states[d_order - 1] = tmp_state;
  if (d_precursors.size()) d_precursors[d_order - 1] = tmp_precursors;

}

//---------------------------------------------------------------------------//
template <class D>
void TimeStepper<D>::update_sources(const double t, const double dt)
{
  // Update the synthetic source.
  d_syntheticsource->build(dt, d_states, d_precursors, d_order);
  for (int cell = 0; cell < d_mesh->number_cells(); ++cell)
  {
    std::cout << " q[" << cell << "]=" << d_syntheticsource->source(cell, 0, 0)
              << std::endl;
  }
  // Update the external sources.
  for (int i = 0; i < d_sources.size(); ++i)
  {
    d_sources[i]->set_time(t);
//    std::cout << " EXT SOURCE " << i << std::endl;
//    for (int cell = 0; cell < d_mesh->number_cells(); ++cell)
//    {
//      std::cout << " q[" << cell << "]=" << d_sources[i]->source(cell, 0)
//                << std::endl;
//    }
  }
}

//---------------------------------------------------------------------------//
template <class D>
void TimeStepper<D>::extrapolate()
{
  // Extrapolate psi(n+1) = 2*psi(n+1/2) - psi(n)
  if (d_discrete)
  {
    for (size_t g = 0; g < d_number_groups; ++g)
    {
      // Clear the scalar flux.  We'll rebuild it with the extrapolated
      // angular flux.
      for (size_t i = 0; i < d_mesh->number_cells(); ++i)
        d_state->phi(g)[i] = 0.0;
      for (size_t o = 0; o < d_quadrature->number_octants(); ++o)
      {
        for (size_t a = 0; a < d_quadrature->number_angles_octant(); ++a)
        {
          int angle = d_quadrature->index(o, a);
          State::angular_flux_type &psi  = d_state->psi(g, o, a);
          const State::angular_flux_type &psi0 = d_states[0]->psi(g, o, a);
          for (size_t i = 0; i < d_mesh->number_cells(); ++i)
          {
            psi[i] = 2.0 * psi[i] - psi0[i];
            if (d_fixup and psi[i] < 0.0) psi[i] = 0.0;
            d_state->phi(g)[i] += d_quadrature->weight(a) * psi[i];
          } // end cell
        } // end angle
      } // end octant
    } // end group
  }
  else
  {
    for (size_t g = 0; g < d_number_groups; ++g)
    {
      State::moments_type &phi  = d_state->phi(g);
      const State::moments_type &phi0 = d_states[0]->phi(g);
      for (size_t i = 0; i < d_mesh->number_cells(); ++i)
        phi[i] = 2.0 * phi[i] - phi0[i];
    } // end group
  }

  if (d_multiply)
  {
    // update precursors
  }

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
