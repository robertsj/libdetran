//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  MGCMDSA.cc
 *  @brief MGCMDSA member definitions
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "MGCMDSA.hh"
#include "callow/solver/LinearSolverCreator.hh"
#include "callow/solver/GaussSeidel.hh"
#include "callow/solver/Richardson.hh"
#include "callow/solver/Jacobi.hh"

#include "utilities/MathUtilities.hh"

#define COUT(c) std::cout << c << std::endl;

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
  : Base(input, material, mesh, ssource, fsource,
         cutoff, include_fission, adjoint, "MG-CMDSA")
{
  if (d_input->check("mgpc_disable_fission") and d_include_fission)
    d_include_fission = (0 == d_input->get<int>("mgpc_disable_fission"));


  d_SF = new MGScatterFissionOperator(d_input,
                                      d_material,
                                      d_mesh,
                                      d_scattersource,
                                      d_fissionsource,
                                      d_group_cutoff,
                                      d_include_fission,
                                      d_adjoint);

  bool use_smoother = false;
  if (d_input->check("mgpc_cmdsa_use_smoothing"))
    use_smoother = (0 != d_input->get<int>("mgpc_cmdsa_use_smoothing"));

  if (use_smoother)
  {
    int n = 2;
    if (d_input->check("mgpc_cmdsa_smoothing_iters"))
      n = d_input->get<int>("mgpc_cmdsa_smoothing_iters");

    double w = 1.0;
    if (d_input->check("mgpc_cmdsa_smoothing_relax"))
      w = d_input->get<double>("mgpc_cmdsa_smoothing_relax");

    // A full diffusion operator for use in a weighted Richardson
    d_smoothing_operator = new Operator(d_input,
                                        d_material,
                                        d_mesh,
                                        d_include_fission,
                                        d_group_cutoff,
                                        d_adjoint,
                                        1.0);

    // Gauss-Seidel for the smoother
    //d_smoothing_solver = new callow::GaussSeidel(0.0, 0.0, n, w);
    //d_smoothing_solver = new callow::Richardson(0.0, 0.0, n, w);
    d_smoothing_solver = new callow::Jacobi(0.0, 0.0, n, w);
    d_smoothing_solver->set_operators(d_smoothing_operator);
    d_smoothing_solver->set_monitor_level(0);
  }
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

  // Print the operator for debugging
  if (d_input->check("mgpc_print_operators"))
  {
    d_operator->print_matlab("diffusion_cm.out");
    d_SF->compute_explicit("SF.out");
  }
}

//----------------------------------------------------------------------------//
void MGCMDSA::apply(Vector &V_h, Vector &V_h_out)
{
  /*
   *  fine mesh dsa:    (I-inv(C)*S)*phi
   *  coarse mesh dsa:  (I-P*inv(C)*R*S)*phi
   */

  // Currently, DSA is only used on the flux moments, not on the boundaries.
  size_t size_moments = d_mesh->number_cells();
  V_h_out.set(0.0);

  //--------------------------------------------------------------------------//
  // Construct scatter source: SV_h <-- S * V_h
  //--------------------------------------------------------------------------//

  // Create the total group source on the fine mesh
  Vector SV_h(size_moments * d_number_active_groups, 0.0);
  d_SF->multiply(V_h, SV_h);

  //--------------------------------------------------------------------------//
  // Restrict scatter source: SV_H <-- R * SV_h
  //--------------------------------------------------------------------------//

  Vector SV_H(d_size_coarse, 0.0);
  d_restrict->multiply(SV_h, SV_H);

  //--------------------------------------------------------------------------//
  // Solve coarse mesh diffusion equation: invC_SV_H <-- inv(C_H) * R * SV_h
  //--------------------------------------------------------------------------//

  Vector invC_SV_H(d_size_coarse, 0.0);
  d_solver->solve(SV_H, invC_SV_H);

  //--------------------------------------------------------------------------//
  // Prolong the result: invC_SV_h <-- P * inv(C_H) * R * SV_h
  //--------------------------------------------------------------------------//

  Vector invC_SV_h(d_size_fine, 0.0);
  d_prolong->multiply(invC_SV_H, invC_SV_h);

  //--------------------------------------------------------------------------//
  // Smooth
  //--------------------------------------------------------------------------//

  if (d_smoothing_solver)
  {
    d_smoothing_solver->solve(SV_h, invC_SV_h);
  }

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
