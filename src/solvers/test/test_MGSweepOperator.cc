//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_MGSweepOperator.cc
 *  @brief Test of MGSolverGMRES
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                     \
        FUNC(test_MGSweepOperator_1g)

#include "TestDriver.hh"
#include "solvers/mg/MGSweepOperator.hh"
#include "solvers/mg/MGTransportSolver.hh"
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
int test_MGSweepOperator_1g(int argc, char *argv[])
{
  typedef MGTransportSolver<_1D> MGS;
  FixedSourceData data = get_fixedsource_data(1, 1, 100);
  data.input->put<std::string>("outer_solver", "GMRES");
  data.input->put<int>("inner_max_iters", 1000000);
  data.input->put<int>("outer_max_iters", 100);
  data.input->put<int>("outer_krylov_group_cutoff",   0);
  FixedSourceManager<_1D> M(data.input, data.material, data.mesh);
  M.setup();
  M.set_solver();
  MGS* solver = dynamic_cast<MGS*>(M.solver().bp());
  TEST(solver);
  MGSweepOperator<_1D> S(M.state(),
                               M.boundary(),
                               solver->sweeper(),
                               solver->sweepsource(),
                               0, false);
  S.compute_explicit("sweeper.out");
  return 0;
}


//----------------------------------------------------------------------------//
//              end of test_MGSweepOperator.cc
//----------------------------------------------------------------------------//
