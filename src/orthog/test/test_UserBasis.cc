//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_UserBasis.cc
 *  @brief Test of UserBasis class
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                   \
        FUNC(test_UserBasis)

#include "utilities/TestDriver.hh"
#include "utilities/Definitions.hh"
#include "orthog/UserBasis.hh"
#include "callow/utils/Initialization.hh"
#include "utilities/DBC.hh"
#include <cmath>
#include <iostream>

using namespace detran_orthog;
using namespace detran_utilities;
using namespace detran_test;
using namespace std;

int main(int argc, char *argv[])
{
  callow_initialize(argc, argv);
  RUN(argc, argv);
  callow_finalize();
}

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

int test_UserBasis(int argc, char *argv[])
{
  UserBasis::Parameters p;
  p.size = 8;
  p.order = 7;
  p.orthonormal = true;

  double ref[] = {0.3535534, 0.5400617, 0.5400617, 0.4308202};

  double data [8][8] =
   {
  		 {0.3535534, 0.3535534, 0.3535534, 0.3535534, 0.3535534, 0.3535534, 0.3535534, 0.3535534},
  		 {0.5400617, 0.3857584, 0.2314550, 0.0771517, -.0771517, -.2314550, -.3857584, -.5400617},
  		 {0.5400617, 0.0771517, -.2314550, -.3857584, -.3857584, -.2314550, 0.0771517, 0.5400617},
  		 {0.4308202, -.3077284, -.4308202, -.1846372, 0.1846372, 0.4308202, 0.3077287, -.4308202},
  		 {0.2820380, -.5237849, -.1208734, 0.3626203, 0.3626203, -.1208734, -.5237849, 0.2820380},
  		 {0.1497862, -.4921546, 0.3637664, 0.3209704, -.3209704, -.3637664, 0.4921546, -.1497862},
  		 {0.0615457, -.3077287, 0.5539117, -.3077287, -.3077287, 0.5539117, -.3077287, 0.0615457},
  		 {0.0170697, -.1194880, 0.3584641, -.5974401, 0.5974401, -.3584641, 0.1194880, -.0170697},
   };

  vec_dbl w(8, 0.0);
  for (int i = 0; i < 8; ++i) w[i] = 1.0;

  p.w = w;
  for (int i = 0; i < 8; ++i)
  {
    vec_dbl temp_vec(8, 0.0);
	for (int j = 0; j < 8; ++j)
	{
	  temp_vec[j] = data[i][j];
	}
	p.db->put<vec_dbl>("vec"+AsString(i), temp_vec);
  }

  OrthogonalBasis::Factory_T::ShowRegistered();
  OrthogonalBasis::SP_basis P = OrthogonalBasis::Create("UserBasis", p);
  P->basis()->display();

  double ref_coef[] = {4.150579e-8, 2.098191, -2.805804e-8, -3.122626e-1,
		  1.702241e-9, 9.195001e-3, -2.498209e-11, -5.897288e-5};

  double ref_appx[] = {1.0, 0.9009688, 0.6234897, 0.2225209,
		  -0.2225209, -0.6234897, -0.9009688, -1.0};

  callow::Vector f(8, 0.0);
  for (int i = 0; i < 8; ++i)
    f[i] = std::cos(i*3.1415926/(8-1));
  callow::Vector ft(8, 0.0);
  callow::Vector fa(8, 0.0);

  P->transform(f, ft);

  f.display("F");
  ft.display("FT");
  for (int i = 0; i < 8; ++i)
  {
    TEST(soft_equiv(ft[i], ref_coef[i], 1.0e-6));
  }

  P->inverse(ft, fa);

  fa.display("FA");
  for (int i = 0; i < 8; ++i)
  {
    TEST(soft_equiv(fa[i], ref_appx[i], 1.0e-6));
  }

  return 0;
}


//----------------------------------------------------------------------------//
//              end of test_UserBasis.cc
//----------------------------------------------------------------------------//
