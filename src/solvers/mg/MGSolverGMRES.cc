//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  MGSolverGMRES.cc
 *  @brief MGSolverGMRES member definitions
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "MGSolverGMRES.hh"
#include "MGDSA.hh"
#include "MGCMDSA.hh"
#include "MGTCDSA.hh"
#include "callow/solver/LinearSolverCreator.hh"

namespace detran
{

//----------------------------------------------------------------------------//
template <class D>
MGSolverGMRES<D>::MGSolverGMRES(SP_state                  state,
                                SP_material               material,
                                SP_boundary               boundary,
                                const vec_externalsource &q_e,
                                SP_fissionsource          q_f,
                                bool                      multiply)
  : Base(state, material, boundary, q_e, q_f, multiply)
  , d_moments_size(0)
  , d_moments_size_group(0)
  , d_boundary_size(0)
  , d_boundary_size_group(0)
  , d_lower(0)
  , d_upper(d_number_groups)
  , d_reflective_solve_iterations(0)
  , d_update_angular_flux(false)
{

  //--------------------------------------------------------------------------//
  // DETERMINE ENERGY GROUP BOUNDS
  //--------------------------------------------------------------------------//

  // Set the bounds for the downscatter GS portion and upscatter Krylov
  // portion.  The default is to use GS on the downscatter block and Krylov
  // on the upscatter block.  Note, if this is a multiplying problem, the
  // cutoff is automatically set to zero, i.e. GS is not used at all.
  d_krylov_group_cutoff = material->upscatter_cutoff(d_adjoint);
  if (d_input->check("outer_krylov_group_cutoff"))
  {
    d_krylov_group_cutoff =
      d_input->template get<int>("outer_krylov_group_cutoff");
    if (!d_adjoint)
    {
      Insist((d_krylov_group_cutoff >= 0) &&
             (d_krylov_group_cutoff <= d_material->upscatter_cutoff(d_adjoint)),
             "Upscatter cutoff must be >= 0 and <= material cutoff");
    }
    else
    {
      Insist((d_krylov_group_cutoff <= d_number_groups) &&
             (d_krylov_group_cutoff >= d_material->upscatter_cutoff(d_adjoint)),
             "Downscatter cutoff must be < # groups and >= material cutoff");
    }
  }
  if (multiply)
  {
    d_krylov_group_cutoff = 0;
    if (d_adjoint) d_krylov_group_cutoff = d_number_groups - 1;
  }

  if (d_adjoint)
  {
    d_lower = d_number_groups - 1;
    d_upper = -1;
  }

  d_number_active_groups = std::abs(d_upper - d_krylov_group_cutoff);

  //--------------------------------------------------------------------------//
  // SETUP SWEEPER FOR MULTIGROUP OPERATOR
  //--------------------------------------------------------------------------//

  d_sweeper = d_wg_solver->get_sweeper();
  d_sweepsource = d_wg_solver->get_sweepsource();

  if (d_number_active_groups)
  {
    //------------------------------------------------------------------------//
    // SETUP SOLVER
    //------------------------------------------------------------------------//

    // Create operator
    d_operator = new Operator_T(d_state,
                                d_boundary,
                                d_sweeper,
                                d_sweepsource,
                                d_krylov_group_cutoff,
                                d_adjoint);
    //d_operator->compute_explicit("transport.out");

    // Create temporary unknown and right hand size vectors
    d_x = new callow::Vector(d_operator->number_rows(), 0.0);
    d_b = new callow::Vector(d_operator->number_rows(), 0.0);

    // Get callow solver parameter database.  If not present, create a new
    // one based on default outer solver parameters.  This way, we can
    // use the basic outer_max_iters and outer_tolerance parameters.
    SP_input db;
    if (d_input->check("outer_solver_db"))
    {
      db = d_input->template get<SP_input>("outer_solver_db");
    }
    else
    {
      db = new detran_utilities::InputDB("mgsolvergmres_db");
      db->template put<double>("linear_solver_rtol", d_tolerance);
      db->template put<double>("linear_solver_atol", d_tolerance);
      db->template put<int>("linear_solver_maxit", d_maximum_iterations);
      db->template put<int>("linear_solver_monitor_level", d_print_level);
      d_input->template put<SP_input>("outer_solver_db", db);
    }
    d_solver = callow::LinearSolverCreator::Create(db);
    Assert(d_solver);

    // Set the transport operator.  Note, no second db argument
    // is given, since that is for setting PC's.  We do that
    // explicitly below.
    d_solver->set_operators(d_operator);

    d_moments_size_group  = d_operator->moments_size();
    d_moments_size        = d_moments_size_group * d_number_groups;
    d_boundary_size_group = d_operator->boundary_size();
    d_boundary_size       = d_boundary_size_group * d_number_groups;

    //------------------------------------------------------------------------//
    // PRECONDITIONER
    //------------------------------------------------------------------------//

    std::string pc_type = "none";
    if (d_input->check("outer_pc_type"))
    {
      pc_type = d_input->template get<std::string>("outer_pc_type");
    }
    std::cout << "Using MG-GMRES with PC-" << pc_type << std::endl;

    size_t pc_side = callow::LinearSolver::LEFT;
    if (d_input->check("outer_pc_side"))
      pc_side = d_input->template get<int>("outer_pc_side");

    // Multigroup DSA
    if (pc_type == "mgdsa")
    {
      Assert(d_sweepsource->get_scatter_source());
      d_pc = new MGDSA(d_input,
                       d_material,
                       d_mesh,
                       d_sweepsource->get_scatter_source(),
                       d_fissionsource,
                       d_krylov_group_cutoff,
                       d_multiply,
                       d_adjoint);
    }
    else if (pc_type == "mgcmdsa")
    {
      d_pc = new MGCMDSA(d_input,
                         d_material,
                         d_mesh,
                         d_sweepsource->get_scatter_source(),
                         d_fissionsource,
                         d_krylov_group_cutoff,
                         d_multiply,
                         d_adjoint);
    }
    else if (pc_type == "mgtcdsa")
    {
      SP_pc P(new MGCMDSA(d_input,
                          d_material,
                          d_mesh,
                          d_sweepsource->get_scatter_source(),
                          d_fissionsource,
                          d_krylov_group_cutoff,
                          d_multiply,
                          d_adjoint));

      d_pc = new MGTCDSA<D>(d_input,
                            d_material,
                            d_mesh,
                            d_sweepsource->get_scatter_source(),
                            d_fissionsource,
                            P,
                            d_operator,
                            d_krylov_group_cutoff,
                            d_multiply,
                            d_adjoint);
    }

    if (d_pc)
      d_solver->set_preconditioner(d_pc, pc_side);

  }

  // Check for boundary update.  In GMRES, this is not the
  // default, since it implies extra sweeps that may not be needed.
  if (d_input->check("compute_boundary_flux"))
  {
    if (d_input->template get<int>("compute_boundary_flux"))
    {
      d_update_angular_flux =
        d_input->template get<int>("compute_boundary_flux");
    }
  }
  if (d_state->store_angular_flux())
    d_update_angular_flux = true;
}

//----------------------------------------------------------------------------//
template <class D>
int MGSolverGMRES<D>::number_sweeps() const
{
  Require(d_sweeper);
  return d_sweeper->number_sweeps();
}

//----------------------------------------------------------------------------//
template <class D>
typename MGSolverGMRES<D>::SP_operator MGSolverGMRES<D>::get_operator()
{
  return d_operator;
}

//----------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//----------------------------------------------------------------------------//

template class MGSolverGMRES<_1D>;
template class MGSolverGMRES<_2D>;
template class MGSolverGMRES<_3D>;

} // end namespace detran

//----------------------------------------------------------------------------//
//              end of MGSolverGMRES.cc
//----------------------------------------------------------------------------//
