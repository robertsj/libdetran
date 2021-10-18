/*
 * test_SlabReactorTransientTransient.cc
 *
 *  Created on: Oct 17, 2021
 *  Author: rabab
 */

#define TEST_LIST  \
    FUNC(test_SlabReactor)\

#include "TestDriver.hh"
#include "Mesh1D.hh"
#include "callow/utils/Initialization.hh"
#include "callow/vector/Vector.hh"
#include "solvers/EigenvalueManager.hh"
#include "solvers/mg/DiffusionLossOperator.hh"
#include "solvers/mg/DiffusionGainOperator.hh"
#include "kinetics/TimeDependentMaterial.hh"
#include "SlabReactor_fixture.hh"

using namespace detran_test;
using namespace detran;
using namespace detran_material;
using namespace detran_geometry;
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


int test_SlabReactor(int argc, char *argv[])
{
  typedef TimeStepper<_1D> TS_1D;
  InputDB::SP_input inp = get_input();
  Mesh1D::SP_mesh mesh = get_mesh(3);

  TimeDependentMaterial::SP_material mat(new SlabMaterial('transport'));
  // steady state
  EigenvalueManager<_1D> manager(inp, mat, mesh);
  manager.solve();
  State::SP_state ic = manager.state();
  mat->set_eigenvalue(ic->eigenvalue());
  mat->update(0, 0, 1, false);

  vec_int matmap = mesh->mesh_map("MATERIAL");

  // Normalize state.
  double F = 0;
  for (int i = 0; i < mesh->number_cells(); ++i)
  {
    int m = matmap[i];
    F += (ic->phi(0)[i]) * mat->sigma_f(m, 0) +
          (ic->phi(1)[i]) * mat->sigma_f(m, 1);
  }

  ic->scale(1.0/F);

  printf(" %1.16f", ic->eigenvalue());

  // transient
  TS_1D stepper(inp, mat, mesh, true);
  stepper.solve(ic);


return 0;
}
