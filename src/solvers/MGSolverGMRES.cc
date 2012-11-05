//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   MGSolverGMRES.cc
 *  @author robertsj
 *  @date   Jun 19, 2012
 *  @brief  MGSolverGMRES member definitions.
 */
//---------------------------------------------------------------------------//

#include "MGSolverGMRES.hh"
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
{

  //-------------------------------------------------------------------------//
  // DETERMINE ENERGY GROUP BOUNDS
  //-------------------------------------------------------------------------//

  // Set the bounds for the downscatter GS portion and upscatter Krylov
  // portion.  The default is to use GS on the downscatter block and Krylov
  // on the upscatter block.
  d_upscatter_cutoff = material->upscatter_cutoff();
  if (d_input->check("outer_upscatter_cutoff"))
  {
    d_upscatter_cutoff = d_input->template get<int>("outer_upscatter_cutoff");
    Insist((d_upscatter_cutoff >= 0) and
           (d_upscatter_cutoff <= d_material->upscatter_cutoff()),
           "Upscatter cutoff must be >= 0 and <= material upscatter cutoff");
  }
  d_upscatter_size = d_number_groups - d_upscatter_cutoff;


  //-------------------------------------------------------------------------//
  // SETUP SWEEPER FOR MULTIGROUP OPERATOR
  //-------------------------------------------------------------------------//

  d_sweeper = d_wg_solver->get_sweeper();
  d_sweepsource = d_wg_solver->get_sweepsource();

  //-------------------------------------------------------------------------//
  // SETUP SOLVER
  //-------------------------------------------------------------------------//

  // Create operator
  d_operator = new MGTransportOperator<D>(d_state,
                                          d_boundary,
                                          d_sweeper,
                                          d_sweepsource,
                                          d_upscatter_cutoff);

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

  d_moments_size_group = d_operator->moments_size();
  d_moments_size = d_moments_size_group * d_number_groups;
  d_boundary_size_group = d_operator->boundary_size();
  d_boundary_size = d_boundary_size_group * d_number_groups;
  //--------------------------------------------------------------------------//
  // PRECONDITIONER
  //--------------------------------------------------------------------------//

  /*
   *  Eventual Options:
   *    1. fine mesh diffusion (= DSA)
   *    2. coarse mesh diffusion (= CMDSA)
   *    3. coarse mesh transport (= ~S2 on coarse grid)
   */

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
