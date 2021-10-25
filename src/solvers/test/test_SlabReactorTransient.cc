/*
 * test_SlabReactorTransientTransient.cc
 *
 *  Created on: Oct 17, 2021
 *  Author: rabab
 */

#define TEST_LIST  \
    FUNC(test_SlabReactor)\

#include "TestDriver.hh"
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

//---------------------------------------------------------------------------//

void test_monitor(void* data, TimeStepper<_1D>* ts, int step, double t,
                  double dt, int it, bool conv)
{
  static double maxp = 0;
  double F = 0;
  TimeStepper<_1D>::vec_int matmap = ts->mesh()->mesh_map("MATERIAL");
  TimeStepper<_1D>::vec_int cmm = ts->mesh()->mesh_map("COARSEMESH");

  vec_dbl &T = ts->multiphysics()->variable(0);

  double T_avg = 0;
  double T_max = 0;
  double V = 0;
  for (int i = 0; i < ts->mesh()->number_cells(); ++i)
  {
    int m = matmap[i];
    F += ts->mesh()->volume(0) *
         (ts->state()->phi(0)[i] * ts->material()->sigma_f(m, 0) +
          ts->state()->phi(1)[i] * ts->material()->sigma_f(m, 1) );
    if (cmm[i] != 4)
    {
      V += ts->mesh()->volume(0);
      T_avg += ts->mesh()->volume(0) * T[i];
      if (T[i] > T_max) T_max = T[i];
    }
  }
  F *= KAPPA / V;
  if (F > maxp && conv) maxp = F;
  T_avg /= 17550.0;

  printf("** %5i  %16.13f  %18.12e  %18.12e  %18.12e  %5i \n", step, t, F, T_avg, T_max, it);
}

//---------------------------------------------------------------------------//


int test_SlabReactor(int argc, char *argv[])
{
  typedef TimeStepper<_1D> TS_1D;
  InputDB::SP_input inp = get_input();
  Mesh1D::SP_mesh mesh = get_mesh(3);

  //-------------------------------------------------------------------------//
  // MATERIAL
  //-------------------------------------------------------------------------//

  bool transport = false;
  if (inp->get<std::string>("equation") != "diffusion") transport = true;
  TS_1D::SP_material mat(new SlabMaterial(mesh, transport, 1));


  //-------------------------------------------------------------------------//
  // STEADY STATE
  //-------------------------------------------------------------------------//

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
  stepper.set_monitor(test_monitor);

  detran_utilities::SP<SlabMaterial> mat_slab;
  mat_slab = mat;

  stepper.set_multiphysics(mat_slab->physics(),
		                    update_T_rhs,
                             (void *) mat_slab.bp());

  stepper.solve(ic);


return 0;
}
