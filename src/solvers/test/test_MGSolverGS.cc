//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_MGSolverGS.cc
 *  @brief Test of MGSolverGS
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                                 \
        FUNC(test_MGSolverGS_1g)                  \
        FUNC(test_MGSolverGS_7g_forward)          \
        FUNC(test_MGSolverGS_7g_forward_multiply) \
        FUNC(test_MGSolverGS_7g_adjoint)          \
        FUNC(test_MGSolverGS_7g_adjoint_multiply)

#include "TestDriver.hh"
#include "solvers/FixedSourceManager.hh"
#include "solvers/test/fixedsource_fixture.hh"

using namespace detran_test;
using namespace detran;
using namespace detran_utilities;
using namespace std;
using std::cout;
using std::endl;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

typedef FixedSourceManager<_1D> Manager;
typedef Manager::SP_manager SP_manager;

void set_data(InputDB::SP_input db)
{
  db->put<std::string>("outer_solver", "GS");
  db->put<double>("inner_tolerance", 1e-14);
  db->put<double>("outer_tolerance", 1e-14);
  db->put<int>("inner_max_iters", 1000000);
  db->put<int>("outer_max_iters", 1000000);
  db->put<std::string>("bc_west", "reflect");
  db->put<std::string>("bc_east", "reflect");
}

SP_manager get_manager(FixedSourceData &data, bool fiss)
{
  SP_manager manager(new Manager(data.input, data.material, data.mesh, fiss));
  manager->setup();
  manager->set_source(data.source);
  manager->set_solver();
  manager->solve();
  return manager;
}

int test_MGSolverGS_1g(int argc, char *argv[])
{
  callow_initialize(argc, argv);
  {
    FixedSourceData data = get_fixedsource_data(1, 1);
    set_data(data.input);
    data.input->put<std::string>("bc_west", "vacuum");
    data.input->put<std::string>("bc_east", "vacuum");
    SP_manager manager = get_manager(data, false);
    TEST(soft_equiv(manager->state()->phi(0)[0], 3.6060798202396613));
  }
  callow_finalize();
  return 0;
}

int test_MGSolverGS_7g_forward(int argc, char *argv[])
{
  callow_initialize(argc, argv);
  {
    FixedSourceData data = get_fixedsource_data(1, 7);
    set_data(data.input);
    SP_manager manager = get_manager(data, false);
    double ref[] = {1.983654685392368e+01, 3.441079047626809e+02,
         5.302787426165165e+01, 1.125133608569081e+01, 2.662710276585539e+01,
         1.010604145062320e+01, 4.015682491688769e+00};
    for (int g = 0; g < 7; ++g)
    {
      printf(" %10.12e  %10.12e \n", ref[g], manager->state()->phi(g)[0]);
      TEST(soft_equiv(ref[g], manager->state()->phi(g)[0]));
    }
  }
  callow_finalize();
  return 0;
}

int test_MGSolverGS_7g_forward_multiply(int argc, char *argv[])
{
  callow_initialize(argc, argv);
  {
    FixedSourceData data = get_fixedsource_data(1, 7);
    set_data(data.input);
    SP_manager manager = get_manager(data, true);
    double ref[] =
    { 3.646729598901197e+02, 5.352648103971697e+03, 3.309487533450470e+02,
        1.856704021668497e+01, 2.763765763929116e+01, 1.018645586459539e+01,
        4.020322297305944e+00 };
    for (int g = 0; g < 7; ++g)
    {
      TEST(soft_equiv(ref[g], manager->state()->phi(g)[0]));
    }
  }
  callow_finalize();
  return 0;
}

int test_MGSolverGS_7g_adjoint(int argc, char *argv[])
{
  callow_initialize(argc, argv);
  {
    FixedSourceData data = get_fixedsource_data(1, 7);
    set_data(data.input);
    data.input->put<int>("adjoint", 1);
    SP_manager manager = get_manager(data, false);
    double ref[] =
    { 1.859700683043188e+02, 1.976212385028667e+02, 3.498588283183322e+01,
        1.129601285153168e+01, 2.693961991801324e+01, 8.478379540507643e+00,
        3.681286723043155e+00 };
    for (int g = 0; g < 7; ++g)
    {
      printf("%4i %20.12e  %20.12e \n", g, ref[g], manager->state()->phi(g)[0]);
      TEST(soft_equiv(ref[g], manager->state()->phi(g)[0]));
    }
  }
  callow_finalize();
  return 0;
}

int test_MGSolverGS_7g_adjoint_multiply(int argc, char *argv[])
{
  callow_initialize(argc, argv);
  {
    FixedSourceData data = get_fixedsource_data(1, 7);
    set_data(data.input);
    data.input->put<int>("adjoint", 1);
    SP_manager manager = get_manager(data, true);
    double ref[] = {8.164042429244962e+02, 6.027576162725934e+02,
                    4.583647621311337e+02, 3.956844300911188e+02,
                    1.145654520054639e+03, 1.333127219999030e+03,
                    1.356688501751730e+03 };
    for (int g = 0; g < 7; ++g)
    {
      TEST(soft_equiv(ref[g], manager->state()->phi(g)[0]));
    }
  }
  callow_finalize();
  return 0;
}

//----------------------------------------------------------------------------//
//              end of test_MGSolverGS.cc
//----------------------------------------------------------------------------//
