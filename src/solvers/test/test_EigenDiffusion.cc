//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_EigenPI.cc
 *  @brief Test of EigenPI
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                              \
        FUNC(test_EigenPI_1g)                  \
        FUNC(test_EigenPI_7g_forward)          \
        FUNC(test_EigenPI_7g_adjoint)

#include "TestDriver.hh"
#include "solvers/EigenvalueManager.hh"
#include "solvers/test/eigenvalue_fixture.hh"
#include "utilities/MathUtilities.hh"

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

int test_EigenPI_1g(int argc, char *argv[])
{
  EigenvalueData data = get_eigenvalue_data(1, 1);
  data.input->put<std::string>("eigen_solver", "PI");
  EigenvalueManager<_1D> manager(data.input, data.material, data.mesh);
  manager.solve();
  TEST(soft_equiv(manager.state()->eigenvalue(), 1.0));
  return 0;
}

int test_EigenPI_7g_forward(int argc, char *argv[])
{
  EigenvalueData data = get_eigenvalue_data(1, 7);
  data.input->put<std::string>("eigen_solver", "PI");
  EigenvalueManager<_1D> manager(data.input, data.material, data.mesh);
  manager.solve();
  double ref[] =
  { 7.641447918995387e-02, 9.959998182719878e-01, 4.630704325129503e-02,
      9.200150683392864e-04, 2.626285787433469e-05, 3.927767793415337e-07,
      8.670823991668230e-09 };
  vec_dbl phi(7, 0.0);
  for (int g = 0; g < 7; ++g) phi[g] = manager.state()->phi(g)[0];
  detran_utilities::vec_scale(phi, 1.0 / detran_utilities::norm(phi, "L2"));
  for (int g = 0; g < 7; ++g)
  {
    printf("%4i %20.12e  %20.12e \n", g, ref[g], phi[g]);
    TEST(soft_equiv(ref[g], phi[g]));
  }
  TEST(soft_equiv(1.038797451683334, manager.state()->eigenvalue()));
  return 0;
}

int test_EigenPI_7g_adjoint(int argc, char *argv[])
{
  EigenvalueData data = get_eigenvalue_data(1, 7);
  data.input->put<std::string>("eigen_solver", "PI");
  data.input->put<int>("adjoint", 1);
  EigenvalueManager<_1D> manager(data.input, data.material, data.mesh);
  manager.solve();
  vec_dbl phi(7, 0.0);
  for (int g = 0; g < 7; ++g) phi[g] = manager.state()->phi(g)[0];
  detran_utilities::vec_scale(phi, 1.0 / detran_utilities::norm(phi, "L2"));
  double ref[] =
  { 3.869207620414996e-01, 2.823112768002036e-01, 2.433579222846369e-01,
      2.493696176828925e-01, 1.037168852716319e-01, 5.397326649609824e-01,
      5.891653761163530e-01 };
  for (int g = 0; g < 7; ++g)
  {
    printf("%4i %20.12e  %20.12e \n", g, ref[g], phi[g]);
    //TEST(soft_equiv(ref[g], phi[g]));
  }
  //TEST(soft_equiv(1.038797451683334, manager.state()->eigenvalue()));
  return 0;
}

//----------------------------------------------------------------------------//
//              end of test_EigenPI.cc
//----------------------------------------------------------------------------//
