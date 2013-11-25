//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  MGDSA.cc
 *  @brief MGDSA member definitions
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "MGDSA.hh"
#include "callow/solver/LinearSolverCreator.hh"
#include "ioutils/StdOutUtils.hh"

namespace detran
{

//----------------------------------------------------------------------------//
MGDSA::MGDSA(SP_input         input,
             SP_material      material,
             SP_mesh          mesh,
             SP_scattersource ssource,
             SP_fissionsource fsource,
             size_t           cutoff,
             bool             include_fission,
             bool             adjoint)
  : Base(input, material, mesh, ssource, fsource, cutoff,
         include_fission, adjoint, "MG-DSA")
{
  if (d_input->check("mgdsa_disable_fission") && d_include_fission)
  {
    d_include_fission = (0 == d_input->get<int>("mgdsa_disable_fission"));
  }

  // Compute the diffusion coefficients.
  d_material->compute_diff_coef();

  // Check for outer preconditioner db
  SP_input db;
  if (d_input->check("outer_pc_db"))
    db = d_input->get<SP_input>("outer_pc_db");

  // Create the loss operator for this group
  d_operator = new Operator_T(d_input,
                              d_material,
                              d_mesh,
                              d_include_fission,
                              d_group_cutoff,
                              d_adjoint, // adjoint
                              1.0);      // keff

  // Create the linear solver for this group.
  d_solver = callow::LinearSolverCreator::Create(db);

  // Set the operators for this group.  The database is used
  // to set the preconditioner parameters for the diffusion solves.
  d_solver->set_operators(d_operator, db);

  // Set the size of the operator
  d_size = d_operator->number_columns();
}

//----------------------------------------------------------------------------//
void MGDSA::apply(Vector &V, Vector &V_out)
{
  // Currently, DSA is only used on the flux moments, not on the boundaries.
  size_t size_moments = d_mesh->number_cells();
  V_out.set(0.0);

  // Copy input vector to a multigroup flux; only the Krylov block is used.
  State::vec_moments_type
    phi(d_number_groups, State::moments_type(size_moments, 0.0));
  for (int g = d_group_cutoff; g < d_number_groups; ++g)
    for (int i = 0; i < size_moments; ++i)
      phi[g][i] = V[(g - d_group_cutoff) * size_moments + i];

  //--------------------------------------------------------------------------//
  // Create the scatter/fission source: S_V <-- S * V
  //--------------------------------------------------------------------------//

  Vector S_V(size_moments * d_number_active_groups, 0.0);
  for (int g = d_group_cutoff; g < d_number_groups; g++)
  {
    State::moments_type source(size_moments, 0.0);
    // Add scatter
    d_scattersource->build_total_group_source(g, d_group_cutoff, phi, source);
    // Add fission
    if (d_include_fission)
      d_fissionsource->build_total_group_source(g, phi, source);
    for (int i = 0; i < size_moments; i++)
      S_V[(g - d_group_cutoff) * size_moments + i] = source[i];
  }

  //--------------------------------------------------------------------------//
  // Solve the diffusion equation:  invC_S_V <-- inv(C) * S * V
  //--------------------------------------------------------------------------//

  Vector invC_S_V(size_moments * d_number_active_groups, 0.0);
  d_solver->solve(S_V, invC_S_V);

  //--------------------------------------------------------------------------//
  // Add the result to the output vector: V_out <-- (I + inv(C) * S) * V
  //--------------------------------------------------------------------------//

  for (int i = 0; i < size_moments * d_number_active_groups; ++i)
    V_out[i] = V[i] + invC_S_V[i];
}

//----------------------------------------------------------------------------//
void MGDSA::build(const double keff, SP_state state)
{
  // Print the operator for debugging
  if (d_input->check("mgpc_print_operators"))
  {
    d_operator->print_matlab("diffusion.out");
  }
}

} // end namespace detran

//----------------------------------------------------------------------------//
//              end of file MGDSA.cc
//----------------------------------------------------------------------------//
