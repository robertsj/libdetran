#define TEST_LIST               \
        FUNC(test_1D_transient_diffusion_matrix)

#include "TestDriver.hh"
#include "Mesh1D.hh"
#include "callow/utils/Initialization.hh"
#include "callow/vector/Vector.hh"
#include "solvers/EigenvalueManager.hh"
#include "solvers/mg/DiffusionLossOperator.hh"
#include "solvers/mg/DiffusionGainOperator.hh"
#include "kinetics/KineticsMaterial.hh"
#include "callow/matrix/Matrix.hh"



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



KineticsMaterial::SP_material get_mat()
  {
    KineticsMaterial::SP_material mat = KineticsMaterial::Create(4, 2, 8);
    // materials are reflector, fuel, control-1, control-2
    //double xs_t[4][2] = {{0.222222,  0.666667}, {0.25641,  0.66667}, {0.25641,  0.1515 }, {0.25641,  0.1515}};
    double xs_t[4][2] = {{0.222222,  0.666667}, {0.25641,  0.66667}, {0.25641,  0.66667}, {0.25641,  0.66667}};
    double xs_f[4][2] = {{0.0, 0.0}, {0.0012, 0.076}, {0.0012, 0.076}, {0.0012, 0.076}};
    double xs_s[4][2][2] = {{{0.19002222, 0.0},
			    {0.032, 0.6566667}},
			   {{0.223910, 0.0},
			    {0.022, 0.55267}},
			    {{0.223910, 0.0},
			    {0.022, 0.515167}},
			    {{0.223910, 0.0},
			    {0.022, 0.515167}}};

    double diff_coeff[4][2] = {{1.5, 0.5}, {1.3, 0.5}, {1.3, 0.5}, {1.3, 0.5}};

    for (int m = 0; m < 4; ++m)
    {
      for (int g = 0; g < 2; ++g)
      {
        if (g == 0 & m > 0) mat->set_chi(m, g, 1.0);
	mat->set_sigma_t(m, g, xs_t[m][g]);
	mat->set_nu(m, g, 2.5);
	mat->set_sigma_f(m, g, xs_f[m][g]);
	mat->set_diff_coef(m, g, diff_coeff[m][g]);

	for (int gp = 0; gp < 2; ++gp)
	{
	  mat->set_sigma_s(m, gp, g, xs_s[m][gp][g]);
	}
      }
    };

   // kinetics
   mat->set_velocity(0,      2.0e7);
   mat->set_velocity(1,      2.2e5);
   double beta[8] = {2.18e-4, 1.02e-03, 6.05e-4, 1.31e-03, 2.20e-03, 6.00e-4, 5.40e-04, 1.52e-04};
   double lambda[8] = {1.246700e-02, 2.829299e-02, 4.2524000e-02, 1.330420e-01, 2.924670e-01, 6.664888e-01, 1.634781e+00, 3.554601e+00};

   for (int i=0; i < 8; i++)
   {
     mat->set_lambda(i, lambda[i]);
     mat->set_beta(i, beta[i]);
     for (int m = 0; m < 4; m++)
     {
       mat->set_chi_d(m, i, 0, 1.0);
     }
   }

    mat->compute_sigma_a();

    mat->finalize();

    return mat;

}

// -------------------------------------------------
Mesh1D::SP_mesh get_mesh(Mesh1D::size_t fmm = 2)
{
    Mesh1D::vec_dbl cm(8);
    for (int i= 0; i < 8; i++)
    {
      cm[i] = 10*i;
    }

    Mesh1D::vec_int fm(7);

    for (int i= 0; i < 7; i++)
    {
      fm[i] = fmm;
    }

    Mesh1D::vec_int mt(7);
    mt[0] = 0; mt[1] = 1; mt[2] = 2;
    mt[3] = 1; mt[4] = 3; mt[5] = 1;
    mt[6] = 0;
    Mesh1D::SP_mesh mesh = Mesh1D::Create(fm, cm, mt);

    return mesh;
}

int test_1D_transient_diffusion_matrix(int argc, char *argv[])
{
  //-------------------------------------------------------------------------//
  // INPUT
  //-------------------------------------------------------------------------//

  InputDB::SP_input inp(new InputDB("1D slab"));
  inp->put<int>("dimension",                      1);
  inp->put<int>("number_groups",                  2);
  inp->put<std::string>("equation",                "diffusion");
  inp->put<std::string>("bc_west",                 "vacuum");
  inp->put<std::string>("bc_east",                 "vacuum");
  inp->put<int>("bc_zero_flux",                   0);

  KineticsMaterial::SP_material mat = get_mat();

  mat->display();
  std::cout << mat->beta_total(0) << "\n";
  Mesh1D::SP_mesh mesh = get_mesh(2);

 // solve the steady state system
  EigenvalueManager<_1D> manager(inp, mat, mesh);
  manager.solve();
  State::SP_state ic = manager.state();
  double keff = ic->eigenvalue();

  DiffusionLossOperator* d_D = new DiffusionLossOperator(inp, mat, mesh,
							 false, 0.0, false, 1.0);
  d_D->print_matlab("loss_mat");

  DiffusionGainOperator* d_G = new DiffusionGainOperator(inp, mat, mesh, false);
  d_G->print_matlab("gain_mat");

  int * rows_L = d_D->rows();
  int* cols_L = d_D->columns();
  double* v_L = d_D->values();

  int * rows_F = d_G->rows();
  int* cols_F = d_G->columns();
  double* v_F = d_G->values();

  //d_D->display(true);
  int num_cells = mesh->number_cells();

  // the big matrix
  int  n = (mat->number_groups() + mat->number_precursor_groups()) * num_cells;
  callow::Matrix A(n, n);
  A.preallocate(668);

  int j;
  double value;
  int r;
  int c;

  const vec_int &mat_map = mesh->mesh_map("MATERIAL");

 //loop over each row
 for (int g=0; g < 2; g++)
 {
   // loop over each energy group
   for (int cell=0; cell < num_cells; cell++)
   {
     int m = mat_map[cell];

     // fill the upper half
     r = num_cells*g + cell; //row number
     for (int p = rows_L[r]; p < rows_L[r + 1]; p++)
     {
       int c = cols_L[p];
       value = -v_L[p]*mat->velocity(g);
       A.insert(r, c, value);
     }

     for (int p = rows_F[r]; p < rows_F[r + 1]; p++)
     {
       c = cols_F[p];
       value = v_F[p]*(1 - mat->beta_total(m))*mat->velocity(g)/keff;
       A.insert(r, c, value, 1);
     }

     // upper left
     //loop over precursor groups
     for (int p=0 ; p < 8; p++)
     {
       c = num_cells*2 + cell + p*num_cells;
       value = mat->lambda(p)*mat->chi_d(m, p, g)*mat->velocity(g);
       A.insert(r, c, value);

       // lower left
       int i = 2*num_cells + cell + num_cells*p;
       c = cell  + num_cells*g;
       double value = mat->beta(m, p)* mat->nu_sigma_f(m, g);
       A.insert(i, c, value);

       // lower right
       if (g == 1)
      {
       c = 2*num_cells + cell + num_cells*p;
       A.insert(i, c, -mat->lambda(p));
      }
    }
  }

 }
 //A.display(true);
 A.assemble();
 A.print_matlab("big_matrix");
 return 0;

 }


