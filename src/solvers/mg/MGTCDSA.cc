//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  MGTCDSA.cc
 *  @brief MGTCDSA member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "MGTCDSA.hh"
#include "solvers/FixedSourceManager.hh"
#include "solvers/mg/MGSolverGMRES.hh"

namespace detran
{

#define COUT(c) std::cout << c << std::endl;

//----------------------------------------------------------------------------//
template <class D>
MGTCDSA<D>::MGTCDSA(SP_input         input,
                    SP_material      material,
                    SP_mesh          mesh,
                    SP_scattersource ssource,
                    SP_fissionsource fsource,
                    SP_base          P,
                    SP_mgto          A,
                    size_t           cutoff,
                    bool             include_fission,
                    bool             adjoint)
  : Base(input, material, mesh, ssource, fsource, cutoff,
         include_fission, adjoint)
  , d_A(A)
  , d_P(P)
  , d_number_coarse_corrections(0)
  , d_use_fine_correction(true)
  , d_A_count(0)
  , d_A_tilde_count(0)
  , d_P_count(0)
{
  Require(d_A);
  Require(d_P);

  if (d_input->check("mgpc_tcdsa_number_coarse_corrections"))
  {
    d_number_coarse_corrections =
      d_input->template get<int>("mgpc_tcdsa_number_coarse_corrections");
  }
  if (d_input->check("mgpc_tcdsa_use_fine_correction"))
    {
      d_use_fine_correction =
        0 != d_input->template get<int>("mgpc_tcdsa_use_fine_correction");
    }

  // Build the approximate transport preconditioner
  {
    // Copy the real input but assign coarser quadrature
    SP_input db = detran_utilities::InputDB::Create();
    *db = *d_input; // copy the input
    int na = 1, np = 1;
    if (d_input->check("mgpc_number_polar_octant"))
      np = d_input->get<int>("mgpc_number_polar_octant");
    if (d_input->check("mgpc_number_azimuth_octant"))
      na = d_input->get<int>("mgpc_number_azimuth_octant");
    db->put<int>("quad_number_azimuth_octant", na);
    db->put<int>("quad_number_polar_octant",   np);
    // Turn off the preconditioner to avoid a nasty infinite recursion
    db->put<std::string>("outer_pc_type", "none");

    // Build a new fixed source solver for the coarse quadrature
    FixedSourceManager<D> manager(db, d_material, d_mesh, d_include_fission);
    manager.setup();
    manager.set_solver();

    // Get the manager's multigroup operator
    MGSolverGMRES<D>* mg_solver =
      dynamic_cast<MGSolverGMRES<D>*>(manager.solver().bp());
    Assert(mg_solver);

    d_A_tilde = mg_solver->get_operator();
    Assert(d_A_tilde);
  }

  d_size = d_P->size();
}

//----------------------------------------------------------------------------//
template <class D>
void MGTCDSA<D>::apply(Vector &b, Vector &x)
{
  apply(b, x, d_number_coarse_corrections);   // x <--  P*b
  if (d_use_fine_correction)
  {
    Vector y(b);
    d_A->multiply(x, y);                      // x <-- A*P*b
    y.add_a_times_x(-2.0, b);                 // x <-- A*P*b - 2b
    y.scale(-1);                              // x <-- (2I-AP)b
    apply(y, x, d_number_coarse_corrections); // x <-- P(2I-AP)b
    ++d_A_count;
  }

//  COUT(" A count  = "  << d_A_count)
//  COUT(" At count = "  << d_A_tilde_count)
//  COUT(" P count  = "  << d_P_count)
//  COUT(" At sweeps = " << d_A_tilde->sweeper()->number_sweeps())
}

//----------------------------------------------------------------------------//
template <class D>
void MGTCDSA<D>::build(const double k, SP_state state)
{
  d_P->build(k, state);
  // update the fission source's eigenvalue if present
  if (d_A_tilde->sweepsource()->get_fission_source())
    d_A_tilde->sweepsource()->get_fission_source()->set_scale(1.0 / k);
}

//----------------------------------------------------------------------------//
template <class D>
void MGTCDSA<D>::apply(Vector &x_in, Vector &x_out, const size_t k)
{
  Vector y(x_in);
  if (k > 0)
  {
    apply(x_in, x_out, k-1);        // x <-- P*x
    d_A_tilde->multiply(x_out, y);  // x <-- A*P*x
    y.add_a_times_x(-2.0, x_in);    // x <-- A*P*x - 2I
    y.scale(-1);                    // x <-- (2I-AP)x
    apply(y, x_out, k-1);           // x <-- P(2I-AP)x
    ++d_A_tilde_count;
  }
  else
  {
    d_P->apply(x_in, x_out);
    ++d_P_count;
  }
}

//----------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//----------------------------------------------------------------------------//

template class MGTCDSA<_1D>;
template class MGTCDSA<_2D>;
template class MGTCDSA<_3D>;

} // end namespace detran

//----------------------------------------------------------------------------//
//              end of file MGTCDSA.cc
//----------------------------------------------------------------------------//




