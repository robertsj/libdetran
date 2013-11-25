//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  MGCMDSA.cc
 *  @brief MGCMDSA member definitions
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "MGCMDSA.hh"
#include "callow/solver/LinearSolverCreator.hh"
#include "utilities/MathUtilities.hh"

namespace detran
{

//----------------------------------------------------------------------------//
MGCMDSA::MGCMDSA(SP_input         input,
                 SP_material      material,
                 SP_mesh          mesh,
                 SP_scattersource ssource,
                 SP_fissionsource fsource,
                 size_t           cutoff,
                 bool             include_fission,
                 bool             adjoint)
  : Base(input, material, mesh, ssource, fsource, cutoff, adjoint, "MG-CMDSA")
  , d_include_fission(include_fission)
{
  if (d_input->check("mgdsa_disable_fission") and d_include_fission)
    d_include_fission = (0 == d_input->get<int>("mgdsa_disable_fission"));

}

//----------------------------------------------------------------------------//
void MGCMDSA::build(const double keff, SP_state state)
{
  // Create coarse mesh and material
  Base::build(keff, state);

  // Check for outer preconditioner db
  SP_input db;
  if (d_input->check("outer_pc_db"))
    db = d_input->get<SP_input>("outer_pc_db");

  // Create coarse mesh diffusion operator
  d_operator = new Operator(d_input,
                            d_c_material,
                            d_coarsemesh,
                            d_include_fission,
                            d_group_cutoff,
                            d_adjoint,
                            keff);

  // Create the linear solver for inverting diffusion operator
  d_solver = callow::LinearSolverCreator::Create(db);

  // Set the operators for this group.  The database is used
  // to set the preconditioner parameters for the diffusion solves.
  d_solver->set_operators(d_operator, db);
}

//----------------------------------------------------------------------------//
void MGCMDSA::apply(Vector &V_h, Vector &V_h_out)
{
  std::cout << " hello" << std::endl;
  /*
   *  fine mesh dsa:    (I-inv(C)*S)*phi
   *  coarse mesh dsa:  (I-P*inv(C)*R*S)*phi
   */

  // Currently, DSA is only used on the flux moments, not on the boundaries.
  size_t size_moments = d_mesh->number_cells();
  V_h_out.set(0.0);
  // Copy input vector to a multigroup flux; only the Krylov block is used.
  State::vec_moments_type
    phi(d_number_groups, State::moments_type(size_moments, 0.0));
  for (int g = d_group_cutoff; g < d_number_groups; ++g)
    for (int i = 0; i < size_moments; ++i)
      phi[g][i] = V_h[(g - d_group_cutoff) * size_moments + i];

  //--------------------------------------------------------------------------//
  // Construct scatter source: SV_h <-- S * V_h
  //--------------------------------------------------------------------------//

  // Create the total group source on the fine mesh
  Vector SV_h(size_moments * d_number_active_groups, 0.0);
  for (int g = d_group_cutoff; g < d_number_groups; g++)
  {
    State::moments_type source(size_moments, 0.0);
    d_scattersource->build_total_group_source(g, d_group_cutoff, phi, source);
    for (int i = 0; i < size_moments; i++)
      SV_h[(g - d_group_cutoff) * size_moments + i] = source[i];
  }

  //--------------------------------------------------------------------------//
  // Restrict scatter source: SV_H <-- R * SV_h
  //--------------------------------------------------------------------------//

  Vector SV_H(d_size_coarse, 0.0);
  d_restrict->multiply(SV_h, SV_H);


  //--------------------------------------------------------------------------//
  // Solve coarse mesh diffusion equation: invC_SV_H <-- inv(C_H) * R * SV_h
  //--------------------------------------------------------------------------//

  Vector invC_SV_H(size_moments * d_number_active_groups, 0.0);
  d_solver->solve(SV_H, invC_SV_H);

  //--------------------------------------------------------------------------//
  // Project the result: invC_SV_h <-- P * inv(C_H) * R * SV_h
  //--------------------------------------------------------------------------//

  Vector invC_SV_h(d_size_fine, 0.0);
  d_restrict->multiply(invC_SV_H, invC_SV_h);

  //--------------------------------------------------------------------------//
  // Add result to output: V_h_out <--- V_h + P * inv(C_H) * R * SV_h
  //--------------------------------------------------------------------------//

  V_h_out.copy(V_h);
  for (int i = 0; i < size_moments * d_number_active_groups; ++i)
    V_h_out[i] += invC_SV_h[i];

}

} // end namespace detran

//----------------------------------------------------------------------------//
//              end of file MGCMDSA.cc
//----------------------------------------------------------------------------//
