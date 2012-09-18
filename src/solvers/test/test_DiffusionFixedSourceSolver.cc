//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_DiffusionFixedSourceSolver.cc
 * \author robertsj
 * \date   Apr 4, 2012
 * \brief  test_DiffusionFixedSourceSolver class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                                \
        FUNC(test_DiffusionFixedSourceSolver_1D)

// Detran headers
#include "TestDriver.hh"
#include "DiffusionFixedSourceSolver.hh"
#include "Mesh1D.hh"
#include "Mesh2D.hh"
#include "Mesh3D.hh"
#include "external_source/ConstantSource.hh"

// Setup
#include "angle/test/quadrature_fixture.hh"
#include "geometry/test/mesh_fixture.hh"
#include "material/test/material_fixture.hh"
#include "external_source/test/external_source_fixture.hh"

using namespace detran_test;
using namespace detran;
using namespace detran_material;
using namespace detran_external_source;
using namespace detran_geometry;
using namespace detran_utilities;
using namespace std;

int main(int argc, char *argv[])
{
  PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
  RUN(argc, argv);
  PetscFinalize();
}

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

int test_DiffusionFixedSourceSolver_1D(int argc, char *argv[])
{

  // Material (kinf = 5/4 = 1.25)
  SP_material mat(new Material(1, 1));
  mat->set_sigma_t(0, 0,    1.0);
  mat->set_sigma_s(0, 0, 0, 0.6);
  mat->set_sigma_f(0, 0,    0.5);
  mat->set_chi(0, 0, 1.0);
  mat->set_diff_coef(0, 0,  0.33);
  mat->compute_sigma_a();
  mat->finalize();
  // Mesh
  vec_dbl cm(2, 0.0); cm[1] = 10.0;
  vec_int fm(1, 10);
  vec_int mat_map(1, 0);
  SP_mesh mesh(new Mesh1D(fm, cm, mat_map));

  typedef BoundaryDiffusion<_1D> B_T;
  typedef B_T::SP_boundary       SP_boundary;
  typedef BoundaryValue<_1D>     BV_T;

  // fixed volume
  if (0){

    // Constant unit source.
    ConstantSource::SP_externalsource q_e(new ConstantSource(1, mesh, 1.0));

    // Input
    DiffusionFixedSourceSolver<_1D>::SP_input input;
    input = new InputDB();
    input->put<int>(     "number_groups",         1);
    input->put<string>(  "bc_west",               "reflect");
    // State
    DiffusionFixedSourceSolver<_1D>::SP_state state(new State(input, mesh));
    // solver
    DiffusionFixedSourceSolver<_1D> solver(input, mat, mesh, state);
    // Build source
    solver.build_source(q_e);
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
      gain += q_e->source(i, 0);
      absorption += state->phi(0)[i] * mat->sigma_a(0, 0);
    }
    // compute the leakage (just the p-current out the vacuum side)
    double leakage = BV_T::value((*boundary)(Mesh::EAST, 0, B_T::OUT));
    // compute net balance (should be 0!!)
    double net = gain - (absorption+leakage);
    cout << "        gain = " << gain << endl;
    cout << "  absorption = " << absorption << endl;
    cout << "     leakage = " << leakage << endl;
    cout << " NET BALANCE = " << gain-(absorption+leakage) << endl;
    TEST(soft_equiv(net, 0.0));

  }

  // boundary source (RMM example
  {

    // Constant unit source.
    ConstantSource::SP_externalsource q_e(new ConstantSource(1, mesh, 1.0));

    // Input
    DiffusionFixedSourceSolver<_1D>::SP_input input;
    input = new InputDB();
    input->put<int>("number_groups",         1);
    input->put<int>("diffusion_fixed_type",  1);

    // State
    DiffusionFixedSourceSolver<_1D>::SP_state state(new State(input, mesh));

    // solver
    DiffusionFixedSourceSolver<_1D> solver(input, mat, mesh, state);

    // extract boundary
    SP_boundary boundary;
    boundary = solver.boundary();
    // set a unit source on the left
    BV_T::value((*boundary)(Mesh::WEST, 0, B_T::IN)) = 1.0;
    solver.build_source();

    // extract the loss operator
    DiffusionFixedSourceSolver<_1D>::SP_lossoperator M = solver.lossoperator();

    /*
     *   simple 1d slab response matrix
     *
     *    R =  | r t |     vac|slab 1|slab 2|slab 3|vac
     *         | t r |
     *    [JoutL, JoutR]' = R*[JinL, JinR]'
     *    JoutR = JinL * t + JinR * r
     *    M =
     *      B 0  0 0  0 0
     *      0 0  1 0  0 0
     *      0 1  0 0  0 0
     *      0 0  0 0  1 0
     *      0 0  0 1  0 0
     *      0 0  0 0  0 B
     *
     */
    double Ji[3][2] = { {0, 0}, {0, 0}, {0, 0} };
    double Jo[3][2] = { {0, 0}, {0, 0}, {0, 0} };

    for (int i = 0; i < 3; i++)
    {
      // set keff and solve
      M->construct(1 + i/3.0);
      M->display();

      solver.solve();
      // define responses
      double r = BV_T::value((*boundary)(Mesh::WEST, 0, B_T::OUT));
      double t = BV_T::value((*boundary)(Mesh::EAST, 0, B_T::OUT));
      double F = 0;
      double A = 0;
      for (int c = 0; c < mesh->number_cells(); c++)
      {
        F += state->phi(0)[i] * mesh->dx(i) * mat->nu_sigma_f(0, 0);
        A += state->phi(0)[i] * mesh->dx(i) * mat->sigma_a(0, 0);
      }

    }

  }

  return 0;
}
//---------------------------------------------------------------------------//
//              end of test_DiffusionFixedSourceSolver.cc
//---------------------------------------------------------------------------//
