//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   test_LinearMaterial.cc
 *  @author Jeremy Roberts
 *  @date   Apr 1, 2012
 *  @brief  Test of LinearMaterial class.
 */
//---------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "kinetics/LinearMaterial.hh"
#include "geometry/Mesh1D.hh"
#include <iostream>
#include <fstream>

using namespace detran;
using namespace detran_utilities;

// Test of basic public interface
TEST(LinearMaterial, Basic)
{
  // Kinetics data
  double lambda[]   = {0.1, 1.0};
  double beta[]     = {0.003, 0.004};
  //double beta_total =  0.007;
  double vel[2]     = {1.0e6, 1.0e4};

  // Cross sections
  double sigma_t[2][2]    = { {1.0, 1.1},  {1.01, 1.11} };
  double sigma_f[2][2]    = { {0.1, 0.6},  {0.11, 0.61} };
  double sigma_s[2][2]    = { {0.5, 0.01}, {0.1, 0.6} };
  double chi[2]           = { 0.1, 0.9};
  // Define the times
  vec_dbl times(2, 2.0);
  times[1] = 10.0;

  // Define the materials
  LinearMaterial::vec_material materials(2);
  materials[0]  = std::make_shared<KineticsMaterial>(1, 2, 2, "time 0");
  materials[1]  = std::make_shared<KineticsMaterial>(1, 2, 2, "time 1");

  // Loop over discrete materials
  for (int m = 0; m < 2; ++m)
  {
    for (int i = 0; i < 2; ++i)
    {
      materials[m]->set_beta(0, i, beta[i]);
      materials[m]->set_lambda(i, lambda[i]);
    }
    for (int g = 0; g < 2; ++g)
    {
      materials[m]->set_velocity(g, vel[g]);
      materials[m]->set_sigma_t(0, g, sigma_t[m][g]);
      materials[m]->set_sigma_f(0, g, sigma_f[m][g]);
      materials[m]->set_sigma_s(0, 0, g, sigma_s[0][g]);
      materials[m]->set_sigma_s(0, 1, g, sigma_s[1][g]);
      materials[m]->set_chi(0, g, chi[g]);
      for (int i = 0; i < 2; ++i)
        materials[m]->set_chi_d(0, g, i, chi[g]);
    }
    materials[m]->finalize();
  }

  // Get the 1g KineticsParameters.
  LinearMaterial poo(times, materials, "combo");

  materials[0]->display();
  materials[1]->display();

  poo.update(5.0, 0.1, 1, true);
  poo.display();
  LinearMaterial::SP_material
    mat(new LinearMaterial(times, materials));
}

//---------------------------------------------------------------------------//
//              end of test_KineticsMaterial.cc
//---------------------------------------------------------------------------//
