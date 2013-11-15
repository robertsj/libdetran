//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_MGSolverGS.cc
 *  @brief Test of MGSolverGS
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                \
        FUNC(test_MGCMFD)

#include "TestDriver.hh"
#include "solvers/FixedSourceManager.hh"
#include "solvers/test/fixedsource_fixture.hh"
#include "solvers/mg/MGSolverGS.hh"
#include "solvers/wg/WGSolverSI.hh"
#include "transport/CurrentTally.hh"
#include "transport/CoarseMesh.hh"
#include "transport/Homogenize.hh"

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
  db->put<int>("inner_max_iters", 1);
  db->put<int>("outer_max_iters", 0);
  db->put<std::string>("bc_west", "vacuum");
  db->put<std::string>("bc_east", "vacuum");
}

int test_MGCMFD(int argc, char *argv[])
{
  callow_initialize(argc, argv);
  {
    int ng = 1;
    FixedSourceData data = get_fixedsource_data(1, ng, 100);
    set_data(data.input);
    data.input->put<std::string>("bc_west", "vacuum");
    data.input->put<std::string>("bc_east", "vacuum");
    data.material->compute_diff_coef();
    data.material->compute_sigma_a();

    // Build manager and extract the solvers
    SP_manager manager(new Manager(data.input, data.material, data.mesh, true));
    manager->setup();
    manager->set_source(data.source);
    manager->set_solver();

    typedef MGSolverGS<_1D> GS;
    typedef WGSolverSI<_1D> SI;
    GS &mgs = *(dynamic_cast<GS* >(&(*manager->solver())));
    SI &wgs = *(dynamic_cast<SI* >(&(*mgs.wg_solver())));

    // Build the coarse mesh and the current tally and set the sweeper
    CoarseMesh::SP_coarsemesh coarse(new CoarseMesh(data.mesh, 10));
    typedef CurrentTally<_1D> Tally;
    Tally::SP_tally tally(new Tally(coarse, manager->quadrature(), ng));
    mgs.sweeper()->set_tally(tally);
    tally->display();

    // Perform one multigroup sweep
    mgs.solve();

    // Homogenize
    Homogenize H(data.material, Homogenize::PHI_D);
    SP_material cmat = H.homogenize(mgs.state(), mgs.mesh(), "COARSEMESH");
    cmat->display();

    //

  }
  callow_finalize();
  return 0;
}

//----------------------------------------------------------------------------//
//              end of test_MGSolverGS.cc
//----------------------------------------------------------------------------//
