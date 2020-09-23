/*
 * test_TransientROMSolver.cc
 *
 *  Created on: Sep 10, 2020
 *      Author: rabab
 */

#define TEST_LIST               \
        FUNC(test_TransientSolver)

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

  //mat->display();

  const char* flux_basis = "/home/rabab/opt/detran/source/src/solvers/test/1D_slab_transient_flux_basis_fm=3";
  const char* precursors_basis = "/home/rabab/opt/detran/source/src/solvers/test/1D_slab_transient_precurses_basis_fm=3";

  SP_matrix basis_f;
  basis_f = new callow::MatrixDense(42, 20);
  ROMBasis::GetBasis(flux_basis, basis_f);

  SP_matrix basis_p;
  basis_p = new callow::MatrixDense(168, 10);
  ROMBasis::GetBasis(precursors_basis, basis_p);

  //
  TransientSolver R(inp, mesh, mat, basis_f, basis_p);

  //
  //get initial state

  EigenvalueManager<_1D> manager(inp, mat, mesh);
  manager.solve();
  State::SP_state ic = manager.state();
  mat->set_eigenvalue(ic->eigenvalue());
  std::cout << " *** updating materials after SS solve ******" << "\n";
  mat->update(0, 0, 1, false);

  // solve
  std::cout << " *** Normalize state ******" << "\n";
 // normalize the vector
  int n = 21;

  callow::Vector phi(42, 0.0);
  for (int g=0; g<2; g++)
  {
   for (int i =0; i< n; i++)
   {
    phi[i + g*n] = ic->phi(g)[i];
   }
  }

  for (int g=0; g<2; g++)
    {
     for (int i =0; i< n; i++)
     {
      phi[i + g*n] = ic->phi(g)[i];
     }
    }

  double norm  = phi.norm();
  for (int g=0; g<2; g++)
    {
     for (int i =0; i< n; i++)
     {
       ic->phi(g)[i] = ic->phi(g)[i]/norm;
     }
    }

  std::cout <<"not normed  = " << phi.norm() << "\n";

  std::cout << norm << "\n";

  // Normalize state.
  double F = 0;
  phi.scale(1/norm);

  std::cout <<"normed  = " << phi.norm() << "\n";

  vec_int matmap = mesh->mesh_map("MATERIAL");
  for (int i = 0; i < mesh->number_cells(); ++i)
  {
    int m = matmap[i];
    F += (ic->phi(0)[i]) * mat->sigma_f(m, 0)+
         (ic->phi(1)[i]) * mat->sigma_f(m, 1);
  }

  ic->scale(1.0/F);
  std::cout << ic->eigenvalue() << "\n";
  printf("eigen value %1.16f " ,ic->eigenvalue()) ;

  std::cout << " ***SOLVE ROM ******" << "\n";

  R.Solve(ic);

  std::cout << " ***reconstruct ******" << "\n";
  std::cout <<  "  ************* END TEST **************************" << "\n";

  R.reconstruct();

  std::cout <<  "  ************* END TEST **************************" << "\n";

return 0;
}





