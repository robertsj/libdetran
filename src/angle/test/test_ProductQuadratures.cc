//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   test_ProductQuadratures.cc
 *  @author Jeremy Roberts
 *  @date   Apr 1, 2012
 *  @brief  Test of GaussLegendre class
 *  @note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                    \
        FUNC(test_ProductQuadratures)

// Detran headers
#include "TestDriver.hh"
#include "QuadratureFactory.hh"
#include "ProductQuadrature.hh"

using namespace detran_angle;
using namespace detran_utilities;
using namespace detran_test;
using namespace std;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

int test_ProductQuadratures(int argc, char *argv[])
{
  typedef SP<ProductQuadrature> SP_pq;

  SP_pq q;

  QuadratureFactory qf;
  InputDB::SP_input db = InputDB::Create();
  db->put<int>("quad_number_polar_octant", 3);
  db->put<int>("quad_number_azimuth_octant", 4);

  // LegendreDTN
  {
    db->put<string>("quad_type", "chebyshevlegendre");
    q = qf.build(db, 2);
    TEST(q);
    q->display();
    for (int a = 0; a < 3; ++a)
    {
      cout << " a=" << a << " phi=" << q->phi(a) << "  "
                         << "cphi=" << q->cos_phi(a)  << "  "
                         << "sphi=" << q->sin_phi(a) << endl;
    }
  }

  return 0;
}


//---------------------------------------------------------------------------//
//              end of test_ProductQuadratures.cc
//---------------------------------------------------------------------------//
