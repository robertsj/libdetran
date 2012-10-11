//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_FixedSourceManager.cc
 * \author robertsj
 * \date   Apr 4, 2012
 * \brief  test_FixedSourceManager class definition.
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                        \
        FUNC(test_FixedSourceManager_1D)

// Detran headers
#include "TestDriver.hh"
#include "FixedSourceManager.hh"
#include "Mesh1D.hh"
#include "Mesh2D.hh"
#include "Mesh3D.hh"
#include "external_source/ConstantSource.hh"
#include "callow/utils/Initialization.hh"

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
  callow_initialize(argc, argv);
  RUN(argc, argv);
  callow_finalize();
}

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

int test_FixedSourceManager_1D(int argc, char *argv[])
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

  typedef BoundarySN<_1D> B_T;
  typedef B_T::SP_boundary       SP_boundary;
  typedef BoundaryValue<_1D>     BV_T;

  // fixed volume
  if (1)
  {

    // Constant unit source.
    ConstantSource::SP_externalsource q_e(new ConstantSource(1, mesh, 1.0));

    // Input
    FixedSourceManager<_1D>::SP_input input;
    input = new InputDB();
    input->put<int>(     "number_groups",         1);
    input->put<string>(  "bc_west",               "reflect");
    input->put<string>(  "problem_type",          "multiply");

    // Manager
    FixedSourceManager<_1D> manager(input, mat, mesh);

    // Build source
    manager.set_source(q_e);

    // Solve.
    manager.solve();

    // Get the state and boundary
    SP_boundary boundary;
    boundary = manager.boundary();
    FixedSourceManager<_1D>::SP_state state;
    state = manager.state();


    // compute total source and absorptions
    double gain = 0;
    double absorption = 0;
    for (int i = 0; i < mesh->number_cells(); i++)
    {
      cout << state->phi(0)[i] << " " << endl;
      gain += q_e->source(i, 0);
      absorption += state->phi(0)[i] * mat->sigma_a(0, 0);
    }
    // compute the leakage (just the p-current out the vacuum side)
    double leakage = 0;//BV_T::value((*boundary)(Mesh::EAST, 0, B_T::OUT));
    // compute net balance (should be 0!!)
    double net = gain - (absorption+leakage);
    cout << "        gain = " << gain << endl;
    cout << "  absorption = " << absorption << endl;
    cout << "     leakage = " << leakage << endl;
    cout << " NET BALANCE = " << gain-(absorption+leakage) << endl;
    TEST(soft_equiv(net, 0.0));

  }

  return 0;
}
//---------------------------------------------------------------------------//
//              end of test_FixedSourceManager.cc
//---------------------------------------------------------------------------//
