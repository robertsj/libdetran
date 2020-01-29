//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_MGDSA.cc
 *  @brief Test of MGSolverGMRES
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                            \
        FUNC(test_MGDSA_1g)                  \
        FUNC(test_MGDSA_7g_forward)          \
        FUNC(test_MGDSA_7g_forward_multiply) \
        FUNC(test_MGDSA_7g_adjoint)          \
        FUNC(test_MGDSA_7g_adjoint_multiply)

#include "TestDriver.hh"
#include "solvers/FixedSourceManager.hh"
#include "solvers/test/fixedsource_fixture.hh"

using namespace detran_test;
using namespace detran;
using namespace detran_utilities;
using namespace std;

int main(int argc, char *argv[])
{
  callow_initialize(argc, argv);
  RUN(argc, argv);
  callow_finalize();
}

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
void set_data(InputDB::SP_input input)
{
  input->put<std::string>("outer_solver", "GMRES");
  input->put<double>("inner_tolerance",   1e-14);
  input->put<double>("outer_tolerance",   1e-14);
  input->put<int>("inner_max_iters",      1000000);
  input->put<int>("outer_max_iters",      100);

  // set pc info
  input->put<string>("outer_pc_type",           "mgdsa");
  input->put<int>("outer_pc_side",              1);

  // outer preconditioner parameters
  InputDB::SP_input db(new InputDB("outer_pc_db"));
  db->put<double>("linear_solver_atol", 0.0);
  db->put<double>("linear_solver_rtol", 1e-14);
  db->put<string>("linear_solver_type", "gmres");
  db->put<int>("pc_side", 1);
  db->put<string>("pc_type", "ilu0");
  db->put<int>("linear_solver_maxit", 1000);
  db->put<int>("linear_solver_gmres_restart", 30);
  db->put<int>("linear_solver_monitor_level", 0);
  input->put<InputDB::SP_input>("outer_pc_db", db);
}

//----------------------------------------------------------------------------//
int test_MGDSA_1g(int argc, char *argv[])
{
  FixedSourceData data = get_fixedsource_data(1, 1);
  set_data(data.input);
  data.input->put<int>("outer_krylov_group_cutoff", 0);
  // solve
  FixedSourceManager<_1D> manager(data.input, data.material, data.mesh);
  manager.setup();
  manager.set_source(data.source);
  manager.set_solver();
  manager.solve();
  double ref = 3.6060798202396702e+00;
  printf(" %24.16e %24.16e \n", ref, manager.state()->phi(0)[0]);
  TEST(soft_equiv(manager.state()->phi(0)[0], 3.6060798202396702e+00));
  return 0;
}

//----------------------------------------------------------------------------//
int test_MGDSA_7g_forward(int argc, char *argv[])
{
  FixedSourceData data = get_fixedsource_data(1, 7);
  set_data(data.input);
  //
  FixedSourceManager<_1D> manager(data.input, data.material, data.mesh);
  manager.setup();
  manager.set_source(data.source);
  manager.set_solver();
  manager.solve();

  double ref[] = {3.8222482411752114e+00, 6.0724967034635284e+00,
      4.8091867341762669e+00, 3.4227342847400251e+00, 4.3296860657059604e+00,
      3.1727425353628869e+00, 2.0010257885859750e+00};

  for (int g = 0; g < 7; ++g)
  {
    printf("%4i %24.16e  %24.16e \n", g, ref[g], manager.state()->phi(g)[0]);
    TEST(soft_equiv(ref[g], manager.state()->phi(g)[0]));
  }

  return 0;
}

//----------------------------------------------------------------------------//
int test_MGDSA_7g_forward_multiply(int argc, char *argv[])
{
  FixedSourceData data = get_fixedsource_data(1, 7);
  set_data(data.input);
  //
  FixedSourceManager<_1D> manager(data.input, data.material, data.mesh, true);
  manager.setup();
  manager.set_source(data.source);
  manager.set_solver();
  manager.solve();
  double ref[] = {1.024469035161e+01, 1.373208084228e+01, 4.885822252151e+00,
      3.423610715484e+00, 4.329712880984e+00, 3.172743507104e+00,
      2.001025825635e+00};
  for (int g = 0; g < 7; ++g)
  {
    printf("%4i %20.12e  %20.12e \n", g, ref[g], manager.state()->phi(g)[0]);
    TEST(soft_equiv(ref[g], manager.state()->phi(g)[0]));
  }
  return 0;
}

//----------------------------------------------------------------------------//
int test_MGDSA_7g_adjoint(int argc, char *argv[])
{
  FixedSourceData data = get_fixedsource_data(1, 7);
  set_data(data.input);
  data.input->put<int>("adjoint", 1);
  //
  FixedSourceManager<_1D> manager(data.input, data.material, data.mesh);
  manager.setup();
  manager.set_source(data.source);
  manager.set_solver();
  manager.solve();
  double ref[] = {4.822821210169e+00, 5.195447583845e+00, 4.808868531433e+00,
      3.462846664756e+00, 4.370889488896e+00, 3.095008727518e+00,
      1.941071200909e+00};
  for (int g = 0; g < 7; ++g)
  {
    printf("%4i %20.12e  %20.12e \n", g, ref[g], manager.state()->phi(g)[0]);
    TEST(soft_equiv(ref[g], manager.state()->phi(g)[0]));
  }
  return 0;
}

//----------------------------------------------------------------------------//
int test_MGDSA_7g_adjoint_multiply(int argc, char *argv[])
{
  FixedSourceData data = get_fixedsource_data(1, 7);
  set_data(data.input);
  data.input->put<int>("adjoint", 1);
  //
  FixedSourceManager<_1D> manager(data.input, data.material, data.mesh, true);
  manager.setup();
  manager.set_source(data.source);
  manager.set_solver();
  manager.solve();
  double ref[] = {5.299599439525e+00, 5.263719813201e+00, 5.277382000990e+00,
      4.410382972775e+00, 5.686798568201e+00, 7.079652532382e+00,
      7.915816059562e+00};
  for (int g = 0; g < 7; ++g)
  {
    printf("%4i %20.12e  %20.12e \n", g, ref[g], manager.state()->phi(g)[0]);
    TEST(soft_equiv(ref[g], manager.state()->phi(g)[0]));
  }
  return 0;
}

//----------------------------------------------------------------------------//
//              end of test_MGDSA.cc
//----------------------------------------------------------------------------//
