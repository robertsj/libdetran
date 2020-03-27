#define TEST_LIST               \
        FUNC(test_1D_matrix)

#include "TestDriver.hh"
#include "Mesh1D.hh"
#include "callow/utils/Initialization.hh"
#include "callow/vector/Vector.hh"
#include "solvers/EigenvalueManager.hh"
#include "solvers/mg/DiffusionLossOperator.hh"
#include "solvers/mg/DiffusionGainOperator.hh"
#include "kinetics/KineticsMaterial.hh"


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



KineticsMaterial::SP_material get_mat(){

	    KineticsMaterial::SP_material mat = KineticsMaterial::Create(4, 2, 8);
		// materials are reflector, fuel, control-1, control-2
		double xs_t[4][2] = {{0.222222,  0.666667}, {0.25641,  0.66667}, {0.25641,  0.1515 }, {0.25641,  0.1515}};
		// double xs_a[4][2] = {{0.0002, 0.01}, {0.0105, 0.1140}, {0.0105, 0.1640}, {0.0105, 0.1640}};
		double xs_f[4][2] = {{0.0, 0.0}, {0.0012, 0.076}, {0.0012, 0.076}, {0.0012, 0.076}};
		double xs_s[4][2][2] = {{{0.190022, 0.0},
								 {0.032, 0.6566667}},
								{{0.223910, 0.0},
								 {0.022, 0.55267}},
								{{0.223910, 0.0},
								 {0.022, 0.55267}},
								 {{0.223910, 0.0},
								 {0.022, 0.55267}}};

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
		mat->set_velocity(0,      1.0e7);
		mat->set_velocity(1,      2.0e5);
		double beta[8] = {2.18e-4, 1.02e-03, 6.05e-4, 1.31e-03, 2.20e-03, 6.00e-4, 5.40e-04, 1.52e-04};
		double lambda[8] = {1.246700e-02, 2.829299e-02, 4.2524000e-02, 1.330420e-01, 2.924670e-01, 6.664888e-01, 1.634781e+00, 3.554601e+00};

		for (int i=0; i < 8; i++)
			{
			  mat->set_lambda(i, lambda[i]);
			  mat->set_beta(i, beta[i]);
			  for (int m = 1; m < 4; m++)
				  {
						mat->set_chi_d(m, i, 0, 1);
				  }

			 }

		 mat->compute_sigma_a();

		 mat->finalize();

		return mat;

	  }

// -------------------------------------------------
Mesh1D::SP_mesh get_mesh(Mesh1D::size_t fmm = 5)
{
    Mesh1D::vec_dbl cm(8);
    cm[0] =  0.0;
    cm[1] = 10.0;
    cm[2] = 20.0;
    cm[3] = 30.0;
    cm[4] = 40.0;
    cm[5] = 50.0;
    cm[6] = 60.0;
    cm[7] = 70.0;
    Mesh1D::vec_int fm(7);
    fm[0] = fmm;
    fm[1] = fmm;
    fm[2] = fmm;
    fm[3] = fmm;
    fm[4] = fmm;
    fm[5] = fmm;
    fm[6] = fmm;
    Mesh1D::vec_int mt(7);
    mt[0] = 0; mt[1] = 1; mt[2] = 2;
    mt[3] = 1; mt[4] = 3; mt[6] = 0;
    Mesh1D::SP_mesh mesh = Mesh1D::Create(fm, cm, mt);

    return mesh;
}

int test_1D_matrix(int argc, char *argv[])
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

  Mesh1D::SP_mesh mesh = get_mesh(5);

  DiffusionLossOperator* d_D = new DiffusionLossOperator(inp, mat, mesh,
									                       false, 0.0, false, 1.0);

  DiffusionGainOperator* d_G = new DiffusionGainOperator(inp, mat, mesh, false);

  d_D->print_matlab();
  d_G->display(true);

   return 0;

   }
