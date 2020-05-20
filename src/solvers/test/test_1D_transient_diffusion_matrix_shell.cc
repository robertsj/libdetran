#define TEST_LIST               \
        FUNC(test_1D_transient_diffusion_matrix_shell)

#include "TestDriver.hh"
#include "callow/utils/Initialization.hh"
#include "callow/vector/Vector.hh"
#include "solvers/EigenvalueManager.hh"
#include "solvers/mg/DiffusionLossOperator.hh"
#include "solvers/mg/DiffusionGainOperator.hh"
#include "kinetics/KineticsMaterial.hh"
#include "callow/matrix/MatrixShell.hh"
#include "callow/vector/Vector.hh"
#include "callow/matrix/Matrix.hh"
#include "1D_transient_matrix_fixture.hh"

using namespace detran_test;
using namespace detran;
using namespace detran_material;
using namespace detran_geometry;
using namespace detran_utilities;
using namespace callow;

using namespace std;
using std::cout;
using std::endl;

int main(int argc, char *argv[])
{
  callow_initialize(argc, argv);
  RUN(argc, argv);
  callow_finalize();
}


class Test_Transient_Matrix_shell_1:public MatrixShell
{
  public:
  typedef DiffusionGainOperator::SP_gainoperator    SP_gainoperator;
  typedef DiffusionLossOperator::SP_lossoperator    SP_lossoperator;
  Test_Transient_Matrix_shell_1(int fmm)
  :MatrixShell(this),
  d_mat(get_mat()),
  d_mesh(get_mesh(fmm)),
  d_input(get_inp())
  {
    k_eff();
  }

  void multiply(const Vector &x,  Vector &y)
  {

   SP_lossoperator d_D (new DiffusionLossOperator(d_input, d_mat, d_mesh, false, 0.0, false, 1.0));

   SP_gainoperator d_G (new DiffusionGainOperator(d_input, d_mat, d_mesh, false));

   const vec_int &mat_map = d_mesh->mesh_map("MATERIAL");

   int r, c;

   int num_cells = d_mesh->number_cells();

   callow::Vector y1(2*num_cells, 0.0);
   callow::Vector y2(2*num_cells, 0.0);

   d_D->multiply(x, y1);
   d_G->multiply(x, y2);

   for (int g=0; g < 2; g++)
   {
     for (int cell=0; cell < num_cells; cell++)
     {
       int m = mat_map[cell];
       r = num_cells*g + cell; //row
       y[r] += (-y1[r] + y2[r]*(1 - d_mat->beta_total(m))/d_keff)*d_mat->velocity(g);
     }  // end of cells loop
   }   // end of energy group loop
  }

 void multiply_transpose(const Vector &x,  Vector &y)
 {
   THROW("matrix transpose operator not implemented");
 }

  private:
  KineticsMaterial::SP_material d_mat;
  Mesh1D::SP_mesh d_mesh;
  InputDB::SP_input d_input;
  double d_keff;
  void k_eff()
  {
    // solve the steady state system
    EigenvalueManager<_1D> manager(d_input, d_mat, d_mesh);
    manager.solve();
    State::SP_state ic = manager.state();
    d_keff = ic->eigenvalue();
  }
};


class Test_Transient_Matrix_shell_2:public MatrixShell
{
  public:
  Test_Transient_Matrix_shell_2(int fmm)
  :MatrixShell(this),
  d_mat(get_mat()),
  d_mesh(get_mesh(fmm)),
  d_input(get_inp())
  {

  }

  void multiply(const Vector &x,  Vector &y)
  {

   const vec_int &mat_map = d_mesh->mesh_map("MATERIAL");

   int r, c;

   int num_cells = d_mesh->number_cells();


   for (int g=0; g < 2; g++)
   {
     for (int cell=0; cell < num_cells; cell++)
     {
       int m = mat_map[cell];
       r = num_cells*g + cell; //row
       for (int p=0 ; p < 8; p++)
       {
         // delayed neutron contribution
         y[r] += d_mat->lambda(p)* d_mat->chi_d(m, p, g)*d_mat->velocity(g)
                *x[p*num_cells + cell];
       }
     } // end of precurors loop
    }  // end of cells loop
   }   // end of energy group loop

 void multiply_transpose(const Vector &x,  Vector &y)
  {
    THROW("matrix transpose operator not implemented");
  }

  private:
  KineticsMaterial::SP_material d_mat;
  Mesh1D::SP_mesh d_mesh;
  InputDB::SP_input d_input;

};

class Test_Transient_Matrix_shell_3:public MatrixShell
{
  public:
  Test_Transient_Matrix_shell_3(int fmm)
  :MatrixShell(this),
  d_mat(get_mat()),
  d_mesh(get_mesh(fmm)),
  d_input(get_inp())
  {

  }

  void multiply(const Vector &x,  Vector &y)
  {

   const vec_int &mat_map = d_mesh->mesh_map("MATERIAL");

   int r, c;

   int num_cells = d_mesh->number_cells();

   for (int g=0; g < 2; g++)
   {
     for (int cell=0; cell < num_cells; cell++)
     {
       int m = mat_map[cell];
       for (int p=0 ; p < 8; p++)
       {
         // add precursors production
         y[p*num_cells + cell] += d_mat->chi(m, g)*d_mat->nu_sigma_f(m, g)*d_mat->beta(m,p)
                                                *y[(g+1)*cell];
       } // end of precurors loop
     }  // end of cells loop
   }   // end of energy group loop

  }
 void multiply_transpose(const Vector &x,  Vector &y)
  {
    THROW("matrix transpose operator not implemented");
  }

  private:
  KineticsMaterial::SP_material d_mat;
  Mesh1D::SP_mesh d_mesh;
  InputDB::SP_input d_input;

};

class Test_Transient_Matrix_shell_4:public MatrixShell
{
  public:
  Test_Transient_Matrix_shell_4(int fmm)
  :MatrixShell(this),
  d_mat(get_mat()),
  d_mesh(get_mesh(fmm)),
  d_input(get_inp())
  {

  }

  void multiply(const Vector &x,  Vector &y)
  {
   const vec_int &mat_map = d_mesh->mesh_map("MATERIAL");

   int r, c;

   int num_cells = d_mesh->number_cells();

   for (int cell=0; cell < num_cells; cell++)
   {
    int m = mat_map[cell];
    for (int p=0 ; p < 8; p++)
    {
     // add precursors decay
     y[p*num_cells + cell] += -d_mat->lambda(p)
                                         *x[p*num_cells + cell];
    }
   } // end of precurors loop
  }  // end of cells loop

 void multiply_transpose(const Vector &x,  Vector &y)
  {
    THROW("matrix transpose operator not implemented");
  }

  private:
  KineticsMaterial::SP_material d_mat;
  Mesh1D::SP_mesh d_mesh;
  InputDB::SP_input d_input;
};


int test_1D_transient_diffusion_matrix_shell(int argc, char *argv[])
{
  int fmm = 2;
  int num_cells = fmm*7;

  Test_Transient_Matrix_shell_1 A_1(fmm);
  Test_Transient_Matrix_shell_2 A_2(fmm);
  Test_Transient_Matrix_shell_3 A_3(fmm);
  Test_Transient_Matrix_shell_4 A_4(fmm);

  callow::Vector X(num_cells*10, 0.0);
  callow::Vector Y(num_cells*10, 0.0);


  callow::Vector Y_1(fmm*7*2, 0.0);
  callow::Vector Y_2(fmm*7*2, 0.0);
  callow::Vector Y_3(fmm*7*8, 0.0);
  callow::Vector Y_4(fmm*7*8, 0.0);

  callow::Vector phi(2*num_cells, 0.0);
  callow::Vector C(8*num_cells, 0.0);


  // arbitrary input vector
  for (int i=0; i<140; i++)
  {
    X[i] = i*0.5;
  }

  // the flux part of the input
  for (int i=0; i < 2*num_cells; i++)
  {
    phi[i] =  X[i];
  }

  for (int i=0; i < 8*num_cells; i++)
  {
   C[i] = X[2*num_cells + i];
  }

  A_1.multiply(phi, Y_1);
  A_2.multiply(C,   Y_2);
  A_3.multiply(phi, Y_3);
  A_4.multiply(C,   Y_4);

  for (int i=0; i < 2*7*fmm   ;i++)
  {
    Y[i] = Y_1[i] + Y_2[i];
  }
  
  for (int i=0; i < 8*7*fmm; i++)
  {
    Y[2*7*fmm +i] = Y_3[i] + Y_4[i];
  }

  double ref[] =  { 7341499939.8,
                    -57199.999999999985,
                    -0.174538,
                    -0.59415279,
                    -1.190672,
                    -4.65647,
                    -12.283613999999998,
                    -32.6579512,
                    -91.547736,
                    -223.939863};

  for (int i = 0; i < 10; i++)
  {
    TEST(soft_equiv(Y[14*i], ref[i], 1.0e-12));
  }

  return 0;
 }






