#define TEST_LIST  \
    FUNC(test_TransientSolver_fom)\
    FUNC(test_TransientSolver_rom)\

#include "TestDriver.hh"
#include "Mesh1D.hh"
#include "callow/utils/Initialization.hh"
#include "callow/vector/Vector.hh"
#include "solvers/EigenvalueManager.hh"
#include "solvers/mg/DiffusionLossOperator.hh"
#include "solvers/mg/DiffusionGainOperator.hh"
#include "kinetics/TimeDependentMaterial.hh"
#include "TransientSolver.hh"
#include "1D_transient_fixture.hh"
#include "ROMBasis.hh"

using namespace detran_test;
using namespace detran;
using namespace detran_material;
using namespace detran_geometry;
using namespace detran_utilities;
using namespace std;
using std::cout;
using std::endl;

typedef callow::MatrixDense::SP_matrix SP_matrix;


int main(int argc, char *argv[])
{
  callow_initialize(argc, argv);
  RUN(argc, argv);
  callow_finalize();
}

int test_TransientSolver_fom(int argc, char *argv[])
{
  typedef TimeStepper<_1D> TS_1D;
  InputDB::SP_input inp = get_input();
  Mesh1D::SP_mesh mesh = get_mesh(3);

  TimeDependentMaterial::SP_material mat(new SlabMaterial('transport'));
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
  time_t begin, end;
  time(&begin);
  TS_1D stepper(inp, mat, mesh, true);
  stepper.solve(ic);

  time(&end);
  time_t elapsed = end - begin;
  printf("time elapsed %1.6f", elapsed);

return 0;
}

int test_TransientSolver_rom(int argc, char *argv[])
{
  typedef TimeStepper<_1D> TS_1D;
  InputDB::SP_input inp = get_input();
  Mesh1D::SP_mesh mesh = get_mesh(3);

  TimeDependentMaterial::SP_material mat(new SlabMaterial('transport'));
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
  printf(" %1.16f \n", ic->eigenvalue());

  time_t begin, end;
  time(&begin);

  const char* flux_basis = "./../../../source/src/solvers/test/rom_basis/1d_transient_flux_basis_r=5";
  const char* precursors_basis = "./../../../source/src/solvers/test/rom_basis/1d_transient_precursors_basis_r=5";

  int r = 5;

  SP_matrix basis_f;
  basis_f = new callow::MatrixDense(42, 2*r);
  ROMBasis::GetBasis(flux_basis, basis_f);

  SP_matrix basis_p;
  basis_p = new callow::MatrixDense(168, r);
  ROMBasis::GetBasis(precursors_basis, basis_p);

  TransientSolver R(inp, mesh, mat, basis_f, basis_p);

  R.Solve(ic);
  time(&end);
  time_t elapsed = end - begin;
  printf("time elapsed %1.6f\n", elapsed);

return 0;
}
