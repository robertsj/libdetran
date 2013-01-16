//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   test_EigenvalueManager.cc
 *  @author robertsj
 *  @date   Apr 4, 2012
 *  @brief  test_EigenvalueManager class definition.
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                           \
        FUNC(test_EigenvalueManager_1D)     \
        FUNC(test_EigenvalueManager_2D)     \
        FUNC(test_EigenvalueManager_3D)

// Detran headers
#include "TestDriver.hh"
#include "EigenvalueManager.hh"
#include "Mesh1D.hh"
#include "Mesh2D.hh"
#include "Mesh3D.hh"
#include "callow/utils/Initialization.hh"
#include "boundary/BoundaryTraits.hh"
#include "boundary/BoundarySN.hh"
#include "utilities/Profiler.hh"

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
using std::cout;
using std::endl;

int main(int argc, char *argv[])
{
  START_PROFILER();
  callow_initialize(argc, argv);
  RUN(argc, argv);
  callow_finalize();
  STOP_PROFILER();
}

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

SP_material test_EigenvalueManager_material()
{
  // Material (kinf = 1.0)
//  SP_material mat(new Material(1, 1));
//  mat->set_sigma_t(0, 0, 1.0);
//  mat->set_sigma_s(0, 0, 0, 0.5);
//  mat->set_sigma_f(0, 0, 0.5);
//  mat->set_chi(0, 0, 1.0);
//  mat->set_diff_coef(0, 0, 0.33);
//  mat->compute_sigma_a();
//  mat->finalize();

  // Total
//  SP_material mat(new Material(1, 2));
//  mat->set_sigma_t(0, 0, 0.2263); // (matid, g, value);
//  mat->set_sigma_t(0, 1, 1.0119);
//  // Fission
//  mat->set_sigma_f(0, 0, 0.0067);
//  mat->set_sigma_f(0, 1, 0.1241);
//  mat->set_chi(0, 0, 1.0);
//  mat->set_chi(0, 1, 0.0);
//  // Scattering
//  mat->set_sigma_s(0, 0, 0, 0.2006); // 1 <- 1
//  mat->set_sigma_s(0, 0, 1, 0.0000); // 1 <- 2
//  mat->set_sigma_s(0, 1, 0, 0.0161); // 2 <- 1
//  mat->set_sigma_s(0, 1, 1, 0.9355); // 2 <- 2
//  mat->compute_diff_coef();
//  mat->finalize();
  // kinf =  0.738135984663063
  //SP_material mat = material_fixture_2g();
// mat->compute_diff_coef();
  // kinf =  0.738135984663063
  SP_material mat = material_fixture_2g();
  mat->compute_sigma_a();
  mat->compute_diff_coef();
  return mat;
}

SP_mesh test_EigenvalueManager_mesh(int d)
{
  vec_dbl cm(2, 0.0); cm[1] = 100.0;
  vec_dbl cc(2, 0.0); cc[1] = 0.1;
  vec_int fm(1, 10);
  vec_int ff(1, 1);
  vec_int mat_map(1, 1);
  SP_mesh mesh;
  if (d == 1) mesh = new Mesh1D(fm, cm, mat_map);
  if (d == 2) mesh = new Mesh2D(fm, fm, cm, cm, mat_map);
  if (d == 3) mesh = new Mesh3D(fm, ff, ff, cm, cc, cc, mat_map);
  Assert(mesh);
  return mesh;
}

InputDB::SP_input test_EigenvalueManager_input()
{
  InputDB::SP_input inp(new InputDB());
  // BASIC INPUT
  inp->put<string>("problem_type",            "eigenvalue");
  inp->put<string>("equation",                "diffusion");
  inp->put<string>("bc_west",                 "reflect");
  inp->put<string>("bc_east",                 "vacuum");
  inp->put<string>("bc_south",                "reflect");
  inp->put<string>("bc_north",                "vacuum");
  inp->put<string>("bc_bottom",               "reflect");
  inp->put<string>("bc_top",                  "reflect");
  inp->put<int>("sweeper_print_octants", 1);
  // QUADRATURE
  inp->put<int>("quad_number_polar_octant",   2);
  inp->put<int>("quad_number_azimuth_octant", 2);
  // INNER
  inp->put<string>("inner_solver",            "SI");
  inp->put<double>("inner_tolerance",         1e-9);
  inp->put<int>("inner_print_level",          0);
  inp->put<double>("inner_tolerance",         1e-9);
  inp->put<int>("inner_max_iters",            10000);

  // OUTER
  inp->put<string>("outer_solver",            "GS");
  inp->put<int>("outer_print_level",          0);
  inp->put<int>("outer_max_iters",            1000);
  inp->put<int>("outer_krylov_group_cutoff",  0);
  inp->put<double>("outer_tolerance",         1e-9);
  // EIGEN
  inp->put<int>("eigen_print_level",          2);
  inp->put<string>("eigen_solver",            "diffusion");
  inp->put<double>("eigen_tolerance",         1e-8);
  inp->put<int>("eigen_print_interval",       1);
  inp->put<int>("eigen_max_iters",            100);
  // CALLOW DATABASE (APPLIED TO ALL)
  {
    InputDB::SP_input callowdb(new InputDB("callowdb"));
    callowdb->put<std::string>("eigen_solver_type",               "slepc");
    callowdb->put<int>("eigen_solver_monitor_level",              1000);
    callowdb->put<int>("eigen_solver_maxit",                      4000);
    callowdb->put<double>("eigen_solver_tol",                     1e-10);
    callowdb->put<std::string>("linear_solver_type",              "petsc");
    callowdb->put<double>("linear_solver_atol",                   1e-12);
    callowdb->put<double>("linear_solver_rtol",                   1e-12);
    callowdb->put<int>("linear_solver_maxit",                     4000);
    callowdb->put<int>("linear_solver_gmres_restart",             80);
    callowdb->put<int>("linear_solver_monitor_level",             0);
    callowdb->put<std::string>("pc_type",                         "petsc_pc");
    callowdb->put<std::string>("petsc_pc_type",                   "lu");
    callowdb->put<int>("petsc_pc_factor_levels",                  4);
    callowdb->put<int>("eigen_solver_monitor_level",              2);
    inp->put<InputDB::SP_input>("eigen_solver_db",  callowdb);
    inp->put<InputDB::SP_input>("outer_solver_db",  callowdb);
    inp->put<InputDB::SP_input>("inner_solver_db",  callowdb);
  }
  return inp;
}

template <class D>
int test_EigenvalueManager_T()
{
  typedef EigenvalueManager<D>       Manager_T;

  // Input, materials, and mesh
  typename Manager_T::SP_input input = test_EigenvalueManager_input();
  SP_material mat = test_EigenvalueManager_material();
  SP_mesh mesh = test_EigenvalueManager_mesh(D::dimension);
  input->template put<int>("number_groups", mat->number_groups());
  input->template put<int>("dimension", D::dimension);

  // Manager
  Manager_T manager(input, mat, mesh);

  // Solve.
  bool flag = manager.solve();
  TEST(flag);

  typename Manager_T::SP_state state = manager.state();
  //state->display();

  return 0;
}

int test_EigenvalueManager_1D(int argc, char *argv[])
{
  return test_EigenvalueManager_T<_1D>();
}

int test_EigenvalueManager_2D(int argc, char *argv[])
{
  return test_EigenvalueManager_T<_2D>();
}

int test_EigenvalueManager_3D(int argc, char *argv[])
{
  return test_EigenvalueManager_T<_3D>();
}

//---------------------------------------------------------------------------//
//              end of test_EigenvalueManager.cc
//---------------------------------------------------------------------------//
