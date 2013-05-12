//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_MGSolverGS.cc
 *  @brief Test of MGSolverGS
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST             \
        FUNC(test_MGSolverGS)

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
  callow_initialize(argc, argv);
  RUN(argc, argv);
  callow_finalize();
}

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

int test_MGSolverGS(int argc, char *argv[])
{
  int flag = 0;

  // 1-D, 1 group
  {
    FixedSourceData data = get_fixedsource_data(1, 1);
    data.input->put<std::string>("outer_solver", "GS");
    data.input->put<double>("inner_tolerance",   1e-14);
    data.input->put<double>("outer_tolerance",   1e-14);
    data.input->put<int>("inner_max_iters",      1000000);
    data.input->put<int>("outer_max_iters",      1000000);
    FixedSourceManager<_1D> manager(data.input, data.material, data.mesh);
    manager.setup();
    manager.set_source(data.source);
    manager.set_solver();
    manager.solve();
    TEST(soft_equiv(manager.state()->phi(0)[0], 3.6060798202396613));
  }

  // 1-D, 7 group, pure reflecting, forward
  {
    FixedSourceData data = get_fixedsource_data(1, 7);
    data.input->put<std::string>("outer_solver", "GS");
    data.input->put<std::string>("bc_west",      "reflect");
    data.input->put<std::string>("bc_east",      "reflect");
    data.input->put<double>("inner_tolerance",   1e-14);
    data.input->put<double>("outer_tolerance",   1e-14);
    data.input->put<int>("inner_max_iters",      1000000);
    data.input->put<int>("outer_max_iters",      1000000);
    FixedSourceManager<_1D> manager(data.input, data.material, data.mesh);
    manager.setup();
    manager.set_source(data.source);
    manager.set_solver();
    manager.solve();
    double ref[] = {1.983654685392368e+01, 3.441079047626809e+02,
                    5.302787426165165e+01, 1.125133608569081e+01,
                    2.662710276585539e+01, 1.010604145062320e+01,
                    4.015682491688769e+00};
    for (int g = 0; g < 7; ++g)
    {
      TEST(soft_equiv(ref[g], manager.state()->phi(g)[0]));
    }
  }

  // 1-D, 7 group, pure reflecting, forward with multiplication
  {
    FixedSourceData data = get_fixedsource_data(1, 7);
    data.input->put<std::string>("outer_solver", "GS");
    data.input->put<std::string>("bc_west",      "reflect");
    data.input->put<std::string>("bc_east",      "reflect");
    data.input->put<double>("inner_tolerance",   1e-14);
    data.input->put<double>("outer_tolerance",   1e-14);
    data.input->put<int>("inner_max_iters",      1000000);
    data.input->put<int>("outer_max_iters",      1000000);
    FixedSourceManager<_1D> manager(data.input, data.material, data.mesh, true);
    manager.setup();
    manager.set_source(data.source);
    manager.set_solver();
    manager.solve();
    //manager.state()->display();
    double ref[] = {3.646729598901197e+02, 5.352648103971697e+03,
                    3.309487533450470e+02, 1.856704021668497e+01,
                    2.763765763929116e+01, 1.018645586459539e+01,
                    4.020322297305944e+00};
    for (int g = 0; g < 7; ++g)
    {
      TEST(soft_equiv(ref[g], manager.state()->phi(g)[0]));
    }
  }

  // 1-D, 7 group, pure reflecting, adjoint
  {
    FixedSourceData data = get_fixedsource_data(1, 7);
    data.input->put<std::string>("outer_solver", "GS");
    data.input->put<std::string>("bc_west",      "reflect");
    data.input->put<std::string>("bc_east",      "reflect");
    data.input->put<double>("inner_tolerance",   1e-14);
    data.input->put<double>("outer_tolerance",   1e-14);
    data.input->put<int>("inner_max_iters",      1000000);
    data.input->put<int>("outer_max_iters",      1000000);
    data.input->put<int>("adjoint",              1);
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
      TEST(soft_equiv(ref[g], manager.state()->phi(g)[0]));
    }
  }

  // 1-D, 7 group, pure reflecting, adjoint with multiplication
  {
    FixedSourceData data = get_fixedsource_data(1, 7);
    data.input->put<std::string>("outer_solver", "GS");
    data.input->put<std::string>("bc_west",      "reflect");
    data.input->put<std::string>("bc_east",      "reflect");
    data.input->put<double>("inner_tolerance",   1e-14);
    data.input->put<double>("outer_tolerance",   1e-14);
    data.input->put<int>("inner_max_iters",      100000);
    data.input->put<int>("outer_max_iters",      1000000);
    data.input->put<int>("adjoint",              1);
    FixedSourceManager<_1D> manager(data.input, data.material, data.mesh, true);
    manager.setup();
    manager.set_source(data.source);
    manager.set_solver();
    manager.solve();
    double ref[] = {8.164042429244962e+02, 6.027576162725934e+02,
                    4.583647621311337e+02, 3.956844300911188e+02,
                    1.145654520054639e+03, 1.333127219999030e+03,
                    1.356688501751730e+03};
    for (int g = 0; g < 7; ++g)
    {
      TEST(soft_equiv(ref[g], manager.state()->phi(g)[0]));
    }
  }

  return flag;
}


//----------------------------------------------------------------------------//
//              end of test_MGSolverGS.cc
//----------------------------------------------------------------------------//
