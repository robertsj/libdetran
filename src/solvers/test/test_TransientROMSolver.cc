#define TEST_LIST                       FUNC(test_TransientSolver)

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

int test_TransientSolver(int argc, char *argv[])
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
  /*
  std::cout << "//-----------------  fom --------------------------------//";
  time_t begin, end;
  time(&begin);
  TS_1D stepper(inp, mat, mesh, true);
  stepper.solve(ic);

  time(&end);
  time_t elapsed = end - begin;

  State::SP_state final = stepper.state();

  SP_matrix phi_mat;
  SP_matrix precursors_mat;
 */
  //phi_mat = stepper.flux_mat;
  //precursors_mat = stepper.precursors_mat;
  //phi_mat->print_matlab("1d_flux.txt");
  //precursors_mat->print_matlab("1d_precursors.txt");

  std::cout << "//-----------------  rom --------------------------------//";

  const char* flux_basis = "/home/rabab/opt/detran/source/src/solvers/test/1D_slab_transient_flux_basis_fm=3";
  const char* precursors_basis = "/home/rabab/opt/detran/source/src/solvers/test/1D_slab_transient_precurses_basis_fm=3";


  int r = 10;
  time_t begin2, end2;
  time(&begin2);
  SP_matrix basis_f;
  basis_f = new callow::MatrixDense(42, 2*r);
  ROMBasis::GetBasis(flux_basis, basis_f);

  SP_matrix basis_p;
  basis_p = new callow::MatrixDense(168, r);
  ROMBasis::GetBasis(precursors_basis, basis_p);

  //
  TransientSolver R(inp, mesh, mat, basis_f, basis_p);

  R.Solve(ic);
  time(&end2);
  time_t elapsed2 = end2 - begin2;

return 0;
}
