/*
 * 1D_transient_matrix_fixture.hh
 *
 *  Created on: Apr 17, 2020
 *  Author: Rabab Elzohery
 */

#ifndef SOLVERS_TEST_1D_TRANSIENT_MATRIX_FIXTURE_HH_
#define SOLVERS_TEST_1D_TRANSIENT_MATRIX_FIXTURE_HH_


#include "callow/utils/Initialization.hh"
#include "solvers/EigenvalueManager.hh"
#include "kinetics/KineticsMaterial.hh"
#include "Mesh1D.hh"


using namespace detran;
using namespace detran_material;
using namespace detran_geometry;
using namespace detran_utilities;

using namespace std;
using std::cout;
using std::endl;



KineticsMaterial::SP_material get_mat()
{
  KineticsMaterial::SP_material mat = KineticsMaterial::Create(4, 2, 8);
  // The materials are reflector, fuel, control 1, control 2.
  // the data given are the diffusion coefficients, out-scattering and absorption cross section,
  // so I'll use them to compute the data needed by Detran.
  // Data for the fuel and the control are the same except for the thermal absorption
  double f = 1.0/3;
  double diff_coeff[2][2] = {{1.5, 0.5}, {1.3, 0.5}};
  double xs_t[2][2] = {{f/diff_coeff[0][0], f/diff_coeff[0][1]},
	                 {f/diff_coeff[1][0], f/diff_coeff[1][1]}};
  double xs_f[2][2] = {{0.0, 0.0}, {0.0012, 0.076}};

  // Inscattering xs = total - outscattering - absorption
  double xs_s[4][2][2] = {{{xs_t[0][0]-0.032-0.0002, 0.0},
    			  {0.032, xs_t[0][1]-0.01}},
    			  {{xs_t[1][0]-0.0220-0.0105, 0.0},
    			   {0.022, xs_t[1][1]-0.1140}},
    			  {{xs_t[1][0]-0.0220-0.0105, 0.0},
    			   {0.022, xs_t[1][1]- 0.1515}},
			   {{xs_t[1][0]-0.0220-0.0105, 0.0},
			    {0.022, xs_t[1][1]-0.1515}}};

  for (int m = 0; m < 4; ++m)
  {
    for (int g = 0; g < 2; ++g)
    {
      if (g == 0 & m > 0) mat->set_chi(m, g, 1.0);
      if (m > 1)
      {
        mat->set_sigma_t(m, g, xs_t[1][g]);
        mat->set_diff_coef(m, g, diff_coeff[1][g]);
	mat->set_sigma_f(m, g, xs_f[1][g]);
      }
      else
      {
        mat->set_sigma_t(m, g, xs_t[m][g]);
        mat->set_diff_coef(m, g, diff_coeff[m][g]);
        mat->set_sigma_f(m, g, xs_f[m][g]);
      }
      mat->set_nu(m, g, 2.5);
      for (int gp = 0; gp < 2; ++gp)
      {
        mat->set_sigma_s(m, gp, g, xs_s[m][gp][g]);
      }
    }
 };

   // kinetics
  mat->set_velocity(0, 2.0e7);
  mat->set_velocity(1, 2.2e5);
  double beta[8] = {2.18e-4, 1.02e-03, 6.05e-4, 1.31e-03,
                    2.20e-03, 6.00e-4, 5.40e-04, 1.52e-04};

  double lambda[8] = {1.246700e-02, 2.829299e-02, 4.2524000e-02, 1.330420e-01,
                      2.924670e-01, 6.664888e-01, 1.634781e+00, 3.554601e+00};

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

 InputDB::SP_input get_inp()
 {
   InputDB::SP_input inp(new InputDB("1D slab"));
   inp->put<int>("dimension",                      1);
   inp->put<int>("number_groups",                  2);
   inp->put<std::string>("equation",                "diffusion");
   inp->put<std::string>("bc_west",                 "vacuum");
   inp->put<std::string>("bc_east",                 "vacuum");
   inp->put<int>("bc_zero_flux",                   0);
   inp->put<double>(  "inner_tolerance",   1e-16);
   inp->put<int>(     "inner_print_out",   0);
   inp->put<double>(  "outer_tolerance",   1e-16);
   inp->put<int>(     "outer_print_out",   0);
   inp->put<int>(     "eigen_max_iters",   100);
   inp->put<double>(  "eigen_tolerance",   1e-16);
   return inp;
 }


#endif /* SOLVERS_TEST_1D_TRANSIENT_MATRIX_FIXTURE_HH_ */
