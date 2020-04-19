//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  test_1D_transient_diffusion_matrix.cc
 *  @brief Test transient diffusion matrix
 *  @note  Copyright(C) 2020 Jeremy Roberts
 */
//---------------------------------------------------------------------------//


#define TEST_LIST               \
        FUNC(test_1D_transient_diffusion_matrix)

#include "TestDriver.hh"
#include "callow/utils/Initialization.hh"
#include "solvers/EigenvalueManager.hh"
#include "solvers/mg/DiffusionLossOperator.hh"
#include "solvers/mg/DiffusionGainOperator.hh"
#include "kinetics/KineticsMaterial.hh"
#include "callow/matrix/Matrix.hh"
#include "callow/vector/Vector.hh"
#include "1D_transient_matrix_fixture.hh"



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

int test_1D_transient_diffusion_matrix(int argc, char *argv[])
{
  //-------------------------------------------------------------------------//
  // INPUT
  //-------------------------------------------------------------------------//
  InputDB::SP_input inp = get_inp();

  KineticsMaterial::SP_material mat = get_mat();
  // mat->display();

  Mesh1D::SP_mesh mesh = get_mesh(2);

 // solve the steady state system to compute keff
  EigenvalueManager<_1D> manager(inp, mat, mesh);
  manager.solve();
  State::SP_state ic = manager.state();
  double keff = ic->eigenvalue();

  DiffusionLossOperator* d_D = new DiffusionLossOperator(inp, mat, mesh,
							 false, 0.0, false, 1.0);
  //d_D->print_matlab("loss_mat");

  DiffusionGainOperator* d_G = new DiffusionGainOperator(inp, mat, mesh, false);
  //d_G->print_matlab("gain_mat");

  int * rows_L = d_D->rows();
  int* cols_L = d_D->columns();
  double* v_L = d_D->values();

  int * rows_F = d_G->rows();
  int* cols_F = d_G->columns();
  double* v_F = d_G->values();

  int num_cells = mesh->number_cells();

  // the big matrix
  int  n = (mat->number_groups() + mat->number_precursor_groups()) * num_cells;
  callow::Matrix A(n, n);
  A.preallocate(668);

  int r, c; //row, column
  double value;

  const vec_int &mat_map = mesh->mesh_map("MATERIAL");

  //loop over energy group
  for (int g=0; g < 2; g++)
  {
   // loop over cells
    for (int cell=0; cell < num_cells; cell++)
    {
      int m = mat_map[cell];
      r = num_cells*g + cell; //row number

      // fill upper left: L - (1-beta)/k *F
      for (int p = rows_L[r]; p < rows_L[r + 1]; p++)
      {
        c = cols_L[p];
        value = -v_L[p]*mat->velocity(g);
        A.insert(r, c, value);
      }

      for (int p = rows_F[r]; p < rows_F[r + 1]; p++)
      {
        c = cols_F[p];
        value = v_F[p]*(1 - mat->beta_total(m))*mat->velocity(g)/keff;
        // add the value
        A.insert(r, c, value, 1);
      }

      //loop over precursor groups
      for (int p=0 ; p < 8; p++)
      {
        // upper right: delayed production
        c = num_cells*2 + cell + p*num_cells;
        value = mat->lambda(p)*mat->chi_d(m, p, g)*mat->velocity(g);
        A.insert(r, c, value);

        // lower left: precursors production
        int i = 2*num_cells + cell + num_cells*p;
        c = cell  + num_cells*g;
        value = mat->beta(m, p)* mat->nu_sigma_f(m, g);
        A.insert(i, c, value);

       // lower right: precursors decay
        if (g == 1)
        {
        c = 2*num_cells + cell + num_cells*p;
        A.insert(i, c, -mat->lambda(p));
        }
      }
    }

  }
  A.assemble();

  // test eigenvalue
  TEST(soft_equiv(ic->eigenvalue(), 1.00080817033683, 1.0e-14))
  //printf("%1.16f", ic->eigenvalue());


  callow::Vector X(n, 1.0);
  callow::Vector Y(n, 0.0);

  A.multiply(X, Y);

  // test 10 point in the resultant vector
  double ref[] =
  {125558366.70909092, -1445.7142857142844, -0.012467, -0.02829299,
  -0.042524, -0.133042, -0.292467, -0.6664888, -1.634781, -3.554601};

  for (int i = 0; i < 10; i++)
  {
    TEST(soft_equiv(Y[14*i], ref[i], 1.0e-12));
  }
  //A.print_matlab("big_matrix");
  return 0;

 }


