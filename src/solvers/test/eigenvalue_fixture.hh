//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   eigenproblem_fixture.hh
 * \author Jeremy Roberts
 * \date   Jul 28, 2012
 * \brief  Eigenproblem for testing.
 */
//---------------------------------------------------------------------------//

#ifndef EIGENPROBLEM_FIXTURE_HH_
#define EIGENPROBLEM_FIXTURE_HH_

// Detran
#include "PowerIteration.hh"
#include "PinCell.hh"

// Detran utilities
#include "DBC.hh"

// Fixtures
#include "angle/test/quadrature_fixture.hh"
#include "material/test/material_fixture.hh"
#include "geometry/test/mesh_fixture.hh"

using namespace detran_test;
using namespace detran_postprocess;
using namespace detran;
using namespace std;

namespace detran_test
{

//typedef detran::Eigensolver<_2D>::SP_solver SP_eigensolver;

//// Return a power iteration solver for the test pin cell
//static SP_eigensolver eigenproblem_fixture()
//{
//  typedef PowerIteration<_2D> solver_T;
//
//  // Get fixtures
//  solver_T::SP_material mat = material_fixture_2g();
//  PinCell::SP_pincell pin = pincell_fixture();
//  solver_T::SP_mesh mesh = pin->mesh();
//  solver_T::SP_quadrature quad = quadruplerange_fixture();
//
//  // Input
//  solver_T::SP_input inp(new InputDB());
//  inp->put<string>("bc_left",   "reflect");
//  inp->put<string>("bc_right",  "reflect");
//  inp->put<string>("bc_bottom", "reflect");
//  inp->put<string>("bc_top",    "reflect");
//
//  // State
//  solver_T::SP_state state(new State(inp, mesh, quad));
//
//  // Boundary
//  solver_T::SP_boundary
//    boundary(new Boundary<_2D>(inp, mesh, quad));
//
//  // Fission source
//  solver_T::SP_fissionsource
//    q_f(new FissionSource(state, mesh, mat));
//
//  // Solver
//  SP_eigensolver
//    solver(new solver_T(inp, state, mesh, mat,
//                        quad, boundary, q_f));
//
//  return solver;
//
//}  // material_fixture_1g



} // end namespace detran_test

#endif /* MATERIAL_FIXTURE_HH_ */

//---------------------------------------------------------------------------//
//              end of material_fixture.hh
//---------------------------------------------------------------------------//
