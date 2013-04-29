//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  TimeStepper.cc
 *  @brief TimeStepper class member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#ifndef detran_TIMESTEPPER_CC_
#define detran_TIMESTEPPER_CC_

#include "TimeStepper.hh"
#include "kinetics/SyntheticDiscreteSource.hh"
#include "kinetics/SyntheticMomentSource.hh"
#include "utilities/MathUtilities.hh"
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
  , d_discrete(false)
  , d_number_groups(0)
  , d_multiply(multiply)
  , d_dt(1.0)
  , d_step_factor(1.0)
  , d_final_time(10.0)
  , d_number_steps(10)
  , d_scheme(BDF1)
  , d_do_output(false)
  , d_fixup(false)
  , d_no_extrapolation(false)
  , d_residual_norm(0.0)
  , d_monitor(NULL)
  , d_monitor_level(1)
  , d_tolerance(1e-4)
  , d_maximum_iterations(1)
  , d_update_multiphysics_rhs(NULL)
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

  // Extract the quadrature, state, and fission.  This lets us fill the
  // state without necessarily building a separate solver.
  d_quadrature = d_solver->quadrature();
  d_state = d_solver->state();
  d_fissionsource = d_solver->fissionsource();
  Ensure(d_state);

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
   Assert(max_steps >= 0);
   if (d_number_steps > max_steps) d_number_steps = max_steps;
  }

  // Get the integration scheme and set the order.
  if (d_input->check("ts_scheme"))
    d_scheme = d_input->template get<int>("ts_scheme");
  Require(d_scheme < END_TIME_SCHEMES);
  d_order = d_scheme;
  if (d_scheme == 0) d_order++;

  if (d_input->check("ts_no_extrapolation") && d_scheme != IMP)
    d_no_extrapolation = d_input->template get<int>("ts_no_extrapolation");

  // Get the convergence criteria
  if (d_input->check("ts_max_iters"))
    d_maximum_iterations = d_input->template get<int>("ts_max_iters");
  if (d_input->check("ts_tolerance"))
    d_tolerance = d_input->template get<double>("ts_tolerance");

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
      d_precursors[i] = new Precursors(d_material->number_precursor_groups(),
                                       d_mesh->number_cells());
    }
  }
  d_state_0 = new State(d_input, d_mesh, d_quadrature);
  if (d_multiply)
  {
    d_precursor   = new Precursors(d_material->number_precursor_groups(),
                                   d_mesh->number_cells());
    d_precursor_0 = new Precursors(d_material->number_precursor_groups(),
                                   d_mesh->number_cells());
  }

  //-------------------------------------------------------------------------//
  // SETUP MONITOR AND OUTPUT
  //-------------------------------------------------------------------------//

  // Default monitor
  set_monitor(ts_default_monitor<D>, NULL);

  if (d_input->check("ts_monitor_level"))
  {
    d_monitor_level = d_input->template get<int>("ts_monitor_level");
  }

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

  // Set the state and initialize the precursors if necessary.  For
  // now, we assume steady state for the first order steps.
  // Update the material, sources, and solver
  d_material->update(0.0, 0, 1, false);
  d_state = initial_state;
  *d_solver->state() = *d_state;

  initialize_precursors();
  *d_states[0] = *d_state;
  if (d_precursors.size()) *d_precursors[0] = *d_precursor;
  if (d_multiphysics) *d_vec_multiphysics[0] = *d_multiphysics;


  // Output the initial state
  if (d_do_output) d_silooutput->write_time_flux(0, d_state, d_discrete);

  // Set the solver
  d_solver->set_solver();

  // Call the monitor, if present.  [data, this, step, time, dt, order, conv]
  if (d_monitor_level) d_monitor(d_monitor_data, this, 0, 0.0, d_dt, 1, true);

  // Perform time steps
  double  t = 0.0;
  double dt = 0.0;
  for (size_t i = 1; i <= d_number_steps; ++i)
  {
    t += d_dt;

    // Determine the order.
    size_t order = d_order;
    if (i < d_order) order = i;

    // Determine extrapolation flag.  By default, we extrapolate
    // if doing the first step of a higher order BDF method.  The
    // user can explicitly turn extrapolation off.
    bool flag = false;
    //if (d_scheme == IMP) flag = true;
    if (d_scheme == IMP || (order == 1 && d_order > 1)) flag = true;
    if (d_no_extrapolation && !(d_scheme == IMP)) flag = false;

    // Set the temporary time step
    dt = d_dt;
    if (flag) dt = 0.5 * d_dt;

    // Perform fixed-point iterations
    size_t iteration = 1;
    for (; iteration <= d_maximum_iterations; ++iteration)
    {
      // Perform the time step
      step(t, dt, order, flag);

      bool converged = check_convergence();
      if (iteration == d_maximum_iterations) converged = true;

      // Call the monitor, if present.
      if (d_monitor_level)
        d_monitor(d_monitor_data, this, i, t, d_dt, iteration, converged);

      if (converged) break;

    } // end iterations

    // Cycle the previous iterates and copy the current solution
    cycle_states_precursors(order);
    *d_states[0] = *d_state;
    if (d_multiply) *d_precursors[0] = *d_precursor;
    if (d_multiphysics) *d_vec_multiphysics[0] = *d_multiphysics;

    // Output the initial state
    if (d_do_output) d_silooutput->write_time_flux(i+1, d_state, true);

  } // end time steps

}

//---------------------------------------------------------------------------//
template <class D>
void TimeStepper<D>::step(const double t,
                          const double dt,
                          const size_t order,
                          const bool   flag)
{
  // Adjust the time for evaluating materials and sources if
  // we're in a half step to be extrapolated.  The dt passed
  // *is* the half step.
  double t_eval = t;
  if (flag) t_eval -= dt;

  // Update the material, sources, and solver
  d_material->update(t_eval, dt, order, true);
  update_sources(t_eval, dt, order);
  d_solver->update();

  // Save old state
  *d_state_0 = *d_state;
  if (d_multiply) *d_precursor_0 = *d_precursor;
  if (d_multiphysics) *d_multiphysics_0 = *d_multiphysics;

  // Solve the MG problem for the new state
  d_solver->solve();
  d_state = d_solver->state();

  // Update the precursors
  update_precursors(t_eval, dt, order);

  // Update the physics.  First, evaluate the right hand side
  // contribution, i.e. dY/dt = rhs(t, Y).  Given this
  // initial value, we multiply by delta_t and any BDF terms.
  if (d_multiphysics) update_multiphysics(t_eval, dt, order);

  if (flag) extrapolate();
}

//---------------------------------------------------------------------------//
template <class D>
void TimeStepper<D>::initialize_precursors()
{
  // Preconditions
  Require(d_state);

  // Skip if we have no multiplication or if we have no precursors
  if (!d_multiply || !d_material->number_precursor_groups()) return;
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
      d_precursor->C(i)[cell] =
        inv_lambda * d_material->beta(mt[cell], i) * fd[cell];
      //printf("%16.9f %16.9f %16.9f %16.9f \n", fd[cell], inv_lambda, d_material->beta(mt[cell], i), d_precursor->C(i)[cell]);
    }
  }


}

//---------------------------------------------------------------------------//
template <class D>
void TimeStepper<D>::update_precursors(const double t,
                                       const double dt,
                                       const size_t order)
{
  // Skip if we have no multiplication
  if (!d_multiply || !d_material->number_precursor_groups()) return;

  // Update the materials to eliminate the synthetic component.
  d_material->update(t, dt, order, false);

  // Update the solver
  d_solver->update();

  // Update the fission density.
  d_fissionsource->update();
  const State::moments_type &fd = d_fissionsource->density();
  const vec_int &mt = d_mesh->mesh_map("MATERIAL");

  /*
   *  The precursors are defined via
   *    f(t, C_i(t)) =  dC/dt = -lambda_i*C_i + beta_i * sum_g X_ig F_g phi_g
   *  For the BDF and IMP schemes, we have
   *   (1/Delta) sum_{j=0}^m a_j * C_i(n+j+1) = -lambda_i*C_i(n+1) + beta_i * sum_g X_ig F_g phi_g
   *   (a_j/Delta + lambda)*C_i(n+1) = -(1/Delta) sum_{j=1}^{m} a_J * C_i(n+j+1) + beta_i * sum_g X_ig F_g phi_g
   */

  for (size_t i = 0; i < d_material->number_precursor_groups(); ++i)
  {
    double lambda = d_material->lambda(i);
    double A = dt / (bdf_coefs[order-1][0]  + dt * lambda);

    for (int cell = 0; cell < d_mesh->number_cells(); ++cell)
    {
      // Add flux contribution
      double value = 0.0;
      value +=  A * d_material->beta(mt[cell], i) * fd[cell];

      // Add previous time step contributions
      for (size_t j = 1; j <= order; ++j)
      {
        //std::cout << " C cont " << d_precursors[j-1]->C(i)[cell] <<  " " << (A / dt) * bdf_coefs[order-1][j] * d_precursors[j-1]->C(i)[cell] << std::endl;
        value += (A / dt) * bdf_coefs[order-1][j] *
                 d_precursors[j-1]->C(i)[cell];
      }
      d_precursor->C(i)[cell] = value;
    }
  }
}

//---------------------------------------------------------------------------//
template <class D>
void TimeStepper<D>::update_multiphysics(const double t,
                                         const double dt,
                                         const size_t order)
{
  // Update the right hand side.  The result is placed into
  // the working vector d_multiphysics
  std::cout << " P before = " << d_multiphysics->variable(0)[0] - 300.0 << std::endl;
  d_update_multiphysics_rhs(d_multiphysics_data, this, t, dt);
  std::cout << " P after = " << d_multiphysics->variable(0)[0] << std::endl;

  // Loop through and compute
  //  y(n+1) = (1/a0) * ( dt*rhs + sum of bdf terms )
  for (size_t i = 0; i < d_multiphysics->number_variables(); ++i)
  {

    // Reference to P(n+1)
    MultiPhysics::vec_dbl &P   = d_multiphysics->variable(i);

    //std::cout << " Pold[0]=" << P[0] << std::endl;
    printf("delP = %18.12e \n", P[0]);
    printf("Pold[0] = %18.12e \n", d_vec_multiphysics[0]->variable(0)[0]);

    // Loop over all elements (usually spatial)
    for (int j = 0; j < P.size(); ++j)
    {

      double v = dt * P[j];
      for (size_t k = 1; k <= order; ++k)
        v += bdf_coefs[order-1][k] * d_vec_multiphysics[k-1]->variable(i)[j];
      P[j] = v / bdf_coefs[order-1][0];
    } // end element loop
    std::cout << " Pnew[0]=" << P[0] - 300.0 << std::endl;
  } // end variable loop
}

//---------------------------------------------------------------------------//
template <class D>
void TimeStepper<D>::cycle_states_precursors(const size_t order)
{
  // Preconditions
  Require(d_states.size() == d_order);
  Require(order <= d_order);

  SP_state        tmp_state;
  SP_precursors   tmp_precursors;
  SP_multiphysics tmp_multiphysics;

  // Save the first element.
  tmp_state = d_states[d_order - 1];
  if (d_precursors.size())
    tmp_precursors = d_precursors[d_order - 1];
  if (d_vec_multiphysics.size())
    tmp_multiphysics = d_vec_multiphysics[d_order - 1];

  for (size_t i = 0; i < d_order - 1; ++i)
  {
    size_t j = d_order - i - 1;

//    std::cout << " j=" << j << " order=" << order
//              << " dorder=" << d_order << " phi=" << d_states[j]->phi(0)[0] << std::endl;


    d_states[j] = d_states[j - 1];
    if (d_precursors.size())
      d_precursors[j] = d_precursors[j-1];
    if (d_vec_multiphysics.size())
      d_vec_multiphysics[j] = d_vec_multiphysics[j-1];

  }
  d_states[0] = tmp_state;
  if (d_precursors.size()) d_precursors[0] = tmp_precursors;
  if (d_vec_multiphysics.size()) d_vec_multiphysics[0] = tmp_multiphysics;

}

//---------------------------------------------------------------------------//
template <class D>
void TimeStepper<D>::update_sources(const double t,
                                    const double dt,
                                    const size_t order)
{
  // Update the synthetic source.
  d_syntheticsource->build(dt, d_states, d_precursors, order);
//  for (int cell = 0; cell < d_mesh->number_cells(); ++cell)
//  {
//    printf(" q[%3i]= %16.12f \n", cell, d_syntheticsource->source(cell, 0, 0));
//  }

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
  //THROW("bad!!!");
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
          // Updated flux
          State::angular_flux_type &psi  = d_state->psi(g, o, a);
          // Previous flux
          const State::angular_flux_type &psi0 = d_states[0]->psi(g, o, a);
          for (size_t i = 0; i < d_mesh->number_cells(); ++i)
          {
            psi[i] = 2.0 * psi[i] - psi0[i];
            if (d_fixup && psi[i] < 0.0) psi[i] = 0.0;
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

  // Extrapolate the precursors.
  if (d_multiply && d_precursor->number_precursor_groups())
  {
    for (size_t i = 0; i < d_precursor->number_precursor_groups(); ++i)
    {
      Precursors::vec_dbl &C = d_precursor->C(i);
      const Precursors::vec_dbl &C0 = d_precursors[0]->C(i);
      for (size_t j = 0; j < d_mesh->number_cells(); ++j)
        C[j] = 2.0 * C[j] - C0[j];
    } // end group
  }

  // Extrapolate multiphysics variables.
  if (d_multiphysics)
  {
    for (size_t i = 0; i < d_multiphysics->number_variables(); ++i)
    {
      MultiPhysics::vec_dbl &P = d_multiphysics->variable(i);
      const Precursors::vec_dbl &P0 = d_vec_multiphysics[0]->variable(i);
      for (size_t j = 0; j < P.size(); ++j)
        P[j] = 2.0 * P[j] - P0[j];
    } // end variable
  }

}

//---------------------------------------------------------------------------//
template <class D>
bool TimeStepper<D>::check_convergence()
{
  // Currently, this implementation only checks convergence
  // based on successive flux iterates.

  double d_residual_norm = 0.0;
  for (size_t g = 0; g < d_material->number_groups(); ++g)
  {
    double nr_g = detran_utilities::
      norm_relative_residual(d_state->phi(g), d_state_0->phi(g));
    d_residual_norm += nr_g * nr_g;
  }
  d_residual_norm = std::sqrt(d_residual_norm);

  if (d_residual_norm < d_tolerance)
    return true;
  return false;

}

//---------------------------------------------------------------------------//
template <class D>
void TimeStepper<D>::
set_multiphysics(SP_multiphysics ic,
                 multiphysics_pointer update_multiphysics_rhs,
                 void* multiphysics_data)
{
  // Preconditions
  Require(ic);
  Require(update_multiphysics_rhs);

  std::cout << " ic T=" << ic->variable(0)[0] << std::endl;

  d_multiphysics            = ic;
  d_update_multiphysics_rhs = update_multiphysics_rhs;
  d_multiphysics_data       = multiphysics_data;

  // Create previous physics states
  d_vec_multiphysics.resize(d_order);
  d_multiphysics_0 = new MultiPhysics(*ic);
  for (int i = 0; i < d_order; ++i)
  {
    d_vec_multiphysics[i] = new MultiPhysics(*ic);
  }
  std::cout << d_vec_multiphysics[0]->variable(0)[0] << std::endl;
  std::cout << " ic T=" << d_multiphysics->variable(0)[0] << std::endl;
}


//---------------------------------------------------------------------------//
template <class D>
void ts_default_monitor(void* data,
                        TimeStepper<D>* ts,
                        int step,
                        double t,
                        double dt,
                        int it,
                        bool converged)
{
  Require(ts);
  if (!ts->monitor_level()) return;
  if (step == 0 && it == 1)
  {
    printf(" step        t       dt   iter \n");
    printf("-------------------------------\n");
  }
  //            0   0.1000   0.0500      1
  printf(" %4i %8.4f %8.4f   %4i \n", step, t, dt, it);
}

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//---------------------------------------------------------------------------//

SOLVERS_INSTANTIATE_EXPORT(TimeStepper<_1D>)
SOLVERS_INSTANTIATE_EXPORT(TimeStepper<_2D>)
SOLVERS_INSTANTIATE_EXPORT(TimeStepper<_3D>)

} // end namespace detran

#endif // TIMESTEPPER_CC_ 

//---------------------------------------------------------------------------//
//              end of file TimeStepper.cc
//---------------------------------------------------------------------------//
