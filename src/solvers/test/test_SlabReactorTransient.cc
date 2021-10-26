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



InputDB::SP_input get_input()
{
  InputDB::SP_input inp(new InputDB("1D slab"));
  inp->put<int>("dimension",                      1);
  inp->put<int>("number_groups",                  2);
  inp->put<std::string>("equation",               "dd");
  inp->put<std::string>("bc_west",                "vacuum");
  inp->put<std::string>("bc_east",                "vacuum");
  inp->put<int>("bc_zero_flux",                   0);
  inp->put<double>("ts_final_time",              60);
  inp->put<double>("ts_step_size",                0.01);
  inp->put<int>("ts_max_steps",                   10000);
  inp->put<int>("ts_output",                      0);
  inp->put<int>("ts_monitor_level",               1);
  inp->put<int>("ts_no_extrapolation",            0);
  inp->put<int>("ts_max_iters",             1000);
  inp->put<double>("ts_tolerance",          1.0e-8);
  inp->put<string>("eigen_solver",                "arnoldi");
  inp->put<string>("outer_solver",                 "GMRES");
  inp->put<int>("outer_krylov_group_cutoff",      0);
  inp->put<double>("eigen_solver_tol", 1e-16);
  inp->put<int>("eigen_solver_maxit",  1000);
  inp->put<string>("inner_solver",          "SI");
  inp->put<double>("outer_tolerance",       1e-14);
  inp->put<double>("inner_tolerance",       1e-14);
  inp->put<int>("inner_max_iters",          1e5);

  inp->put<int>("inner_print_level",        0);
  inp->put<int>("outer_print_level",        1);
  // inner gmres parameters
  InputDB::SP_input db(new InputDB("inner_solver_db"));
  db->put<double>("linear_solver_atol",                 1e-16);
  db->put<double>("linear_solver_rtol",                 1e-15);
  db->put<string>("linear_solver_type",                 "gmres");
  db->put<string>("pc_type",                            "petsc_pc");
  db->put<string>("petsc_pc_type",                      "lu");
  db->put<int>("linear_solver_maxit",                   2000);
  db->put<int>("linear_solver_gmres_restart",           30);
  db->put<int>("linear_solver_monitor_level",           0);
  db->put<string>("eigen_solver_type",                  "slepc");
  db->put<double>("eigen_solver_tol",                   1e-15);
  inp->put<InputDB::SP_input>("inner_solver_db",        db);
  inp->put<InputDB::SP_input>("outer_solver_db",        db);
  inp->put<InputDB::SP_input>("eigen_solver_db",        db);
  inp->put<InputDB::SP_input>("rom_solver_db",          db);
  inp->put<int>("compute_boundary_flux",                1);
  inp->put<int>("linear_solver_maxit",   1000);
  inp->put<int>("linear_solver_monitor_level", 0);
  inp->put<int>("linear_solver_monitor_diverge", 0);
  if (inp->get<std::string>("equation") != "diffusion")
  {
    inp->put<int>("ts_discrete",              1);
    inp->put<int>("store_angular_flux",       1);
  }

return inp;
}


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
