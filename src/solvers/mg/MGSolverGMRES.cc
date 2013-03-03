//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   MGSolverGMRES.cc
 *  @author robertsj
 *  @date   Jun 19, 2012
 *  @brief  MGSolverGMRES member definitions.
 */
//---------------------------------------------------------------------------//

#include "MGSolverGMRES.hh"
#include "MGDSA.hh"
#include "callow/solver/LinearSolverCreator.hh"

namespace detran
{

//---------------------------------------------------------------------------//
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
  , d_reflective_solve_iterations(0)
  , d_update_boundary_flux(false)
{

  //-------------------------------------------------------------------------//
  // DETERMINE ENERGY GROUP BOUNDS
  //-------------------------------------------------------------------------//

  // Set the bounds for the downscatter GS portion and upscatter Krylov
  // portion.  The default is to use GS on the downscatter block and Krylov
  // on the upscatter block.
  d_krylov_group_cutoff = material->upscatter_cutoff();
  if (d_input->check("outer_krylov_group_cutoff"))
  {
    d_krylov_group_cutoff = d_input->template get<int>("outer_krylov_group_cutoff");
    Insist((d_krylov_group_cutoff >= 0) and
           (d_krylov_group_cutoff <= d_material->upscatter_cutoff()),
           "Upscatter cutoff must be >= 0 and <= material upscatter cutoff");
  }
  d_number_active_groups = d_number_groups - d_krylov_group_cutoff;

  //-------------------------------------------------------------------------//
  // SETUP SWEEPER FOR MULTIGROUP OPERATOR
  //-------------------------------------------------------------------------//

  d_sweeper = d_wg_solver->get_sweeper();
  d_sweepsource = d_wg_solver->get_sweepsource();

  if (d_number_active_groups)
  {
    //-----------------------------------------------------------------------//
    // SETUP SOLVER
    //-----------------------------------------------------------------------//

    // Create operator
    d_operator = new Operator_T(d_state,
                                d_boundary,
                                d_sweeper,
                                d_sweepsource,
                                d_krylov_group_cutoff);

    //d_operator->compute_explicit("MGTO.out");

    // Create temporary unknown and right hand size vectors
    d_x = new callow::Vector(d_operator->number_rows(), 0.0);
    d_b = new callow::Vector(d_operator->number_rows(), 0.0);

    // Get callow solver parameter database
    SP_input db;
    if (d_input->check("outer_solver_db"))
    {
      db = d_input->template get<SP_input>("outer_solver_db");
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

    if (pc_type == "mgdsa")
    {
      Assert(d_sweepsource->get_scatter_source());
      d_pc = new MGDSA(d_input,
                       d_material,
                       d_mesh,
                       d_sweepsource->get_scatter_source(),
                       d_krylov_group_cutoff,
                       d_multiply);
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
      d_update_boundary_flux =
        d_input->template get<int>("compute_boundary_flux");
    }
  }

}

//---------------------------------------------------------------------------//
template <class D>
int MGSolverGMRES<D>::number_sweeps() const
{
  Require(d_sweeper);
  return d_sweeper->number_sweeps();
}

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//---------------------------------------------------------------------------//

template class MGSolverGMRES<_1D>;
template class MGSolverGMRES<_2D>;
template class MGSolverGMRES<_3D>;

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of MGSolverGMRES.cc
//---------------------------------------------------------------------------//
