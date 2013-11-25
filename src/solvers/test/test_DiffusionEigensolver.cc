//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  test_DiffusionEigensolver.cc
 *  @brief Test of DiffusionEigensolver class
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                        \
        FUNC(test_DiffusionEigensolver)

#include "utilities/TestDriver.hh"
#include "solvers/DiffusionEigensolver.hh"
#include "geometry/Mesh1D.hh"
#include "geometry/Mesh2D.hh"
#include "geometry/Mesh3D.hh"
#include "callow/utils/Initialization.hh"
#include <iostream>

using namespace detran_test;
using namespace detran;
using namespace detran_material;
using namespace detran_external_source;
using namespace detran_geometry;
using namespace detran_utilities;
using namespace std;

int main(int argc, char *argv[])
{
  callow_initialize(argc, argv);
  RUN(argc, argv);
  callow_finalize();
}

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

int test_DiffusionEigensolver(int argc, char *argv[])
{

  // Material (kinf = 5/4 = 1.25)
  DiffusionEigensolver<_1D>::SP_material mat(new Material(1, 1));
  mat->set_sigma_t(0, 0,    1.0);
  mat->set_sigma_s(0, 0, 0, 0.6);
  mat->set_sigma_f(0, 0,    0.5);
  mat->set_chi(0, 0, 1.0);
  mat->set_diff_coef(0, 0,  1.0/3.0);
  mat->compute_sigma_a();
  mat->finalize();
  // Mesh
  vec_dbl cm(2, 0.0); cm[1] = 10.0;
  vec_int fm(1, 10);
  vec_int mat_map(1, 0);
  DiffusionEigensolver<_1D>::SP_mesh mesh(new Mesh1D(fm, cm, mat_map));

  typedef BoundaryDiffusion<_1D> B_T;
  typedef B_T::SP_boundary       SP_boundary;
  typedef BoundaryValue<_1D>     BV_T;

  // fixed volume
  if (1)
  {

    // Input
    DiffusionEigensolver<_1D>::SP_input input;
    input = new InputDB();
    input->put<int>(     "number_groups",               1);
    input->put<string>(  "bc_west",                     "reflect");
    input->put<string>(  "diffusion_eigensolver_type",  "power");
    // State
    DiffusionEigensolver<_1D>::SP_state state(new State(input, mesh));
    // solver
    DiffusionEigensolver<_1D> solver(input, mat, mesh, state);
    // Solve.
    solver.solve();
    // Get the boundary
    SP_boundary boundary;
    boundary = solver.boundary();
    // compute total source and absorptions
    double gain = 0;
    double absorption = 0;

    for (int i = 0; i < mesh->number_cells(); i++)
    {
      gain += state->phi(0)[i] * mat->sigma_a(0, 0) / state->eigenvalue();
      absorption += state->phi(0)[i] * mat->sigma_a(0, 0);
    }
    // compute the leakage (just the p-current out the vacuum side)
    double leakage = BV_T::value((*boundary)(Mesh::EAST, 0, B_T::OUT));
    // compute net balance (should be 0!!)
    double net = gain - (absorption+leakage);
    cout << "  eigenvalue = " << state->eigenvalue() << endl;
    cout << "        gain = " << gain << endl;
    cout << "  absorption = " << absorption << endl;
    cout << "     leakage = " << leakage << endl;
    cout << " NET BALANCE = " << gain-(absorption+leakage) << endl;
    TEST(soft_equiv(net, 0.0));

  }

  return 0;
}
//---------------------------------------------------------------------------//
//              end of test_DiffusionEigensolver.cc
//---------------------------------------------------------------------------//
