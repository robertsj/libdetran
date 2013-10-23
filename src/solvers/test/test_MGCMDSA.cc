//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_MGCMDSA.cc
 *  @brief Test of MGSolverGMRES
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                            \
        FUNC(test_MGCMDSA_1g)                  \
        FUNC(test_MGCMDSA_7g_forward)          \
        FUNC(test_MGCMDSA_7g_forward_multiply) \
        FUNC(test_MGCMDSA_7g_adjoint)          \
        FUNC(test_MGCMDSA_7g_adjoint_multiply)

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
InputDB::SP_input get_pcdb()
{
  InputDB::SP_input db(new InputDB("outer_pc_db"));
  db->put<double>("linear_solver_atol", 1e-15);
  db->put<double>("linear_solver_rtol", 1e-15);
  db->put<string>("linear_solver_type", "gmres");
  db->put<int>("pc_side", 1);
  db->put<string>("pc_type", "ilu0");
  db->put<int>("linear_solver_maxit", 1000);
  db->put<int>("linear_solver_gmres_restart", 30);
  db->put<int>("linear_solver_monitor_level", 0);
  return db;
}

//----------------------------------------------------------------------------//
int test_MGCMDSA_1g(int argc, char *argv[])
{
  FixedSourceData data = get_fixedsource_data(1, 1);
  data.input->put<std::string>("outer_solver", "GMRES");
  data.input->put<double>("inner_tolerance", 1e-14);
  data.input->put<double>("outer_tolerance", 1e-14);
  data.input->put<int>("inner_max_iters", 1000000);
  data.input->put<int>("outer_max_iters", 100);
  // set pc info
  data.input->put<string>("outer_pc_type", "mgcmdsa");
  data.input->put<int>("outer_pc_side", 1);
  data.input->put<int>("outer_krylov_group_cutoff",   0);
  // outer preconditioner parameters
  InputDB::SP_input pcdb =  get_pcdb();
  data.input->put<InputDB::SP_input>("outer_pc_db", pcdb);
  // solve
  FixedSourceManager<_1D> manager(data.input, data.material, data.mesh);
  manager.setup();
  manager.set_source(data.source);
  manager.set_solver();
  manager.solve();
  cout << " " << manager.state()->phi(0)[0] << endl;
  //manager.state()->display();
  //TEST(soft_equiv(manager.state()->phi(0)[0], 3.6060798202396613));
  return 0;
}

int test_MGCMDSA_7g_forward(int argc, char *argv[])
{
  FixedSourceData data = get_fixedsource_data(1, 7);
  //
  data.input->put<std::string>("outer_solver", "GS");
  data.input->put<double>("inner_tolerance", 1e-14);
  data.input->put<double>("outer_tolerance", 1e-14);
  data.input->put<int>("inner_max_iters", 1000000);
  data.input->put<int>("outer_max_iters", 1000000);
  data.input->put<int>("outer_print_level", 2);
  // set pc info
  data.input->put<int>("outer_pc_side",                     1);
  data.input->put<int>("outer_krylov_group_cutoff",         0);
  data.input->put<int>("print_transport_operator",          0);
  data.input->put<int>("print_preconditioner_operator",     0);
  // outer preconditioner parameters
  InputDB::SP_input pcdb = get_pcdb();
  data.input->put<InputDB::SP_input>("outer_pc_db",   pcdb);

  data.input->put<int>("mgpc_cmdsa_use_smoothing",      1);
  data.input->put<int>("mgpc_cmdsa_smoothing_iters",    2);
  data.input->put<double>("mgpc_cmdsa_smoothing_relax", 1.0);
  //data.input->put<string>("outer_pc_type",              "mgcmdsa");
  data.input->put<int>("mgpc_coarse_mesh_level",        2);

  FixedSourceManager<_1D> manager(data.input, data.material, data.mesh);
  manager.setup();
  manager.set_source(data.source);
  manager.set_solver();
  manager.solve();

  double ref[] = {5.9656743969764925e+00, 1.3299132737782095e+01,
      7.7114982411329898e+00, 4.1562283068907915e+00, 6.8097447636815422e+00,
      3.9593934960362973e+00, 2.1230735301683157e+00};

  for (int g = 0; g < 7; ++g)
  {
    printf("%4i %24.16e  %24.16e \n", g, ref[g], manager.state()->phi(g)[0]);
    TEST(soft_equiv(ref[g], manager.state()->phi(g)[0]));
  }

  return 0;
}

int test_MGCMDSA_7g_forward_multiply(int argc, char *argv[])
{
  FixedSourceData data = get_fixedsource_data(1, 7);
  data.input->put<std::string>("outer_solver", "GMRES");
  data.input->put<std::string>("bc_west", "reflect");
  data.input->put<std::string>("bc_east", "reflect");
  data.input->put<double>("inner_tolerance", 1e-14);
  data.input->put<double>("outer_tolerance", 1e-14);
  data.input->put<int>("inner_max_iters", 1000000);
  data.input->put<int>("outer_max_iters", 1000000);
  data.input->put<int>("outer_print_level", 2);
  FixedSourceManager<_1D> manager(data.input, data.material, data.mesh, true);
  manager.setup();
  manager.set_source(data.source);
  manager.set_solver();
  manager.solve();
  double ref[] = {3.646729598901197e+02, 5.352648103971697e+03,
                  3.309487533450470e+02, 1.856704021668497e+01,
                  2.763765763929116e+01, 1.018645586459539e+01,
                  4.020322297305944e+00 };
  for (int g = 0; g < 7; ++g)
  {
    printf("%4i %20.12e  %20.12e \n", g, ref[g], manager.state()->phi(g)[0]);
    TEST(soft_equiv(ref[g], manager.state()->phi(g)[0]));
  }
  return 0;
}

int test_MGCMDSA_7g_adjoint(int argc, char *argv[])
{
  FixedSourceData data = get_fixedsource_data(1, 7);
  data.input->put<std::string>("outer_solver", "GMRES");
  data.input->put<std::string>("bc_west", "reflect");
  data.input->put<std::string>("bc_east", "reflect");
  data.input->put<double>("inner_tolerance", 1e-14);
  data.input->put<double>("outer_tolerance", 1e-14);
  data.input->put<int>("inner_max_iters", 1000000);
  data.input->put<int>("outer_max_iters", 1000000);
  data.input->put<int>("adjoint", 1);
  data.input->put<int>("outer_krylov_group_cutoff", 6);
  FixedSourceManager<_1D> manager(data.input, data.material, data.mesh);
  manager.setup();
  manager.set_source(data.source);
  manager.set_solver();
  manager.solve();
  double ref[] = {1.859700683043188e+02, 1.976212385028667e+02,
                  3.498588283183322e+01, 1.129601285153168e+01,
                  2.693961991801324e+01, 8.478379540507643e+00,
                  3.681286723043155e+00};
  for (int g = 0; g < 7; ++g)
  {
    printf("%4i %20.12e  %20.12e \n", g, ref[g], manager.state()->phi(g)[0]);
    TEST(soft_equiv(ref[g], manager.state()->phi(g)[0]));
  }
  return 0;
}

int test_MGCMDSA_7g_adjoint_multiply(int argc, char *argv[])
{
  FixedSourceData data = get_fixedsource_data(1, 7);
  data.input->put<std::string>("outer_solver", "GMRES");
  data.input->put<std::string>("bc_west", "reflect");
  data.input->put<std::string>("bc_east", "reflect");
  data.input->put<double>("inner_tolerance", 1e-14);
  data.input->put<double>("outer_tolerance", 1e-14);
  data.input->put<int>("inner_max_iters", 100000);
  data.input->put<int>("outer_max_iters", 1000000);
  data.input->put<int>("adjoint", 1);
  FixedSourceManager<_1D> manager(data.input, data.material, data.mesh, true);
  manager.setup();
  manager.set_source(data.source);
  manager.set_solver();
  manager.solve();
  double ref[] = {8.164042429244962e+02, 6.027576162725934e+02,
                  4.583647621311337e+02, 3.956844300911188e+02,
                  1.145654520054639e+03, 1.333127219999030e+03,
                  1.356688501751730e+03 };
  for (int g = 0; g < 7; ++g)
  {
    printf("%4i %20.12e  %20.12e \n", g, ref[g], manager.state()->phi(g)[0]);
    TEST(soft_equiv(ref[g], manager.state()->phi(g)[0]));
  }
  return 0;
}

//----------------------------------------------------------------------------//
//              end of test_MGCMDSA.cc
//----------------------------------------------------------------------------//
