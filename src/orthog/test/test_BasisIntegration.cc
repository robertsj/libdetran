//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_BasisIntegration.cc
 *  @brief Test of test_BasisIntegration class
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//
// LIST OF TEST FUNCTIONS
#define TEST_LIST                   \
        FUNC(test_BasisIntegration)

#include "angle/QuadratureFactory.hh"
#include "angle/ProductQuadrature.hh"
#include "utilities/TestDriver.hh"
#include "utilities/Definitions.hh"
#include "orthog/OrthogonalBasis.hh"
#include "callow/utils/Initialization.hh"
#include <cmath>

using namespace detran_orthog;
using namespace detran_angle;
using namespace detran_utilities;
using namespace detran_test;
using namespace std;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
  callow_finalize();
}

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

int test_BasisIntegration(int argc, char *argv[])
{
  callow_initialize(argc, argv);

  // 1-D integration of Jacobi moments.
  //
  // We expand psi(mu) in the Jacobi polynomials for mu in [0, 1].  Internally,
  // this is scaled to [-1, 1].  Hence, we expect a DGL quadrature to work
  // very well.
  if (1)
  {
    cout << "---  1-D  ---" << endl;

    vector<string> qtypes;
    qtypes.push_back("gl");
    qtypes.push_back("dgl");
    qtypes.push_back("gc");
    qtypes.push_back("dgc");
    qtypes.push_back("asdr");
    qtypes.push_back("tg");

    for (int qt = 0; qt < qtypes.size(); ++qt)
    {
      cout << qtypes[qt] << endl;
      for (int nn = 3; nn < 7; ++nn)
      {
        int n = nn;

        // Define the 1-D quadrature to use
        InputDB::SP_input db(new InputDB());
        db->put<string>("quad_type", qtypes[qt]);
        db->put<int>("quad_number_polar_octant", n);
        QuadratureFactory::SP_quadrature q = QuadratureFactory::build(db, 1);

        // Define the basis set
        OrthogonalBasis::Parameters p;
        p.size = n;
        p.order = 2;
        p.x.resize(p.size, 0.0);
        p.qw.resize(p.size, 0.0);
        p.lower_bound = 0.0;
        p.upper_bound = 1.0;
        for (int i = 0; i < n; ++i)
        {
          p.x[i] = q->mu(0, i);
          p.qw[i] = q->weight(i);
        }
        OrthogonalBasis::SP_basis P = OrthogonalBasis::Create("jacobi", p);

        // Isotropic flux, expansion coefficient, and approximation
        callow::Vector f(p.size, 1.0);
        callow::Vector ft(p.order + 1, 0.0);
        callow::Vector fa(p.size, 0.0);
        P->transform(f, ft);

        callow::Vector ref(p.order + 1, 0.0);
        callow::Vector err(ft);
        ref[0] = 0.5;
        err.subtract(ref);
        printf(" %4i & %6.3e & %6.3e & %6.3e \\\\\n ", n, err[0], err[1], err[2]);
      }
      printf("\\\\ \n");
    }
  }
//
//  // Integration of Legendre zeroth order moment
//  if (0) {
//    vector<string> qtypes;
//    qtypes.push_back("gl");
//    qtypes.push_back("dgl");
//    qtypes.push_back("uniform");
//    qtypes.push_back("asdr");
//
//    for (int qt = 0; qt < qtypes.size(); ++qt)
//    {
//      for (int nn = 1; nn < 5; ++nn)
//      {
//        int n = nn;
//        printf(" %4i ", 2 * n);
//        // Define the 1-D quadrature to use
//        InputDB::SP_input db(new InputDB());
//        db->put<string>("quad_type", qtypes[qt]);
//        db->put<int>("quad_number_polar_octant", n);
//        QuadratureFactory::SP_quadrature q = QuadratureFactory::build(db, 1);
//
//        // Define the basis set
//        OrthogonalBasis::Parameters p;
//        p.size  = n; p.order = 2;
//        p.x.resize(p.size,  0.0); p.qw.resize(p.size, 0.0);
//        p.lower_bound = -1.0; p.upper_bound = 1.0;
//        int i = 0;
//        for (int o = 0; o < 2; ++o)
//        {
//          for (int a = 0; a < n; ++a, ++i)
//          {
//            int aa = a;
//            if (o == 0) aa = n-a-1;
//            p.x[i] = -q->mu(o, aa); p.qw[i] = q->weight(aa);
//            std::cout << i << " " << p.x[i] << " "<< p.qw[i] << std::endl;
//          }
//        }
//
////        for (int p = 0; p < p.size; ++p)
////        {
////          p.x[p]  = q->cos_theta(p);
////          p.qw[p] = q->polar_weight(p);
////        }
//
//        OrthogonalBasis::SP_basis P = OrthogonalBasis::Create("clp", p);
//        P->weights()->display("W");
//        // Isotropic flux, expansion coefficient, and approximation
//        callow::Vector f(p.size, 1.0);
//        callow::Vector ft(p.order + 1, 0.0);
//        callow::Vector fa(p.size, 0.0);
//
//        P->transform(f, ft);
//        callow::Vector ref(p.order + 1, 0.0);
//        callow::Vector err(ft);
//        ref[0] = 2.0;
//        err.subtract(ref);
//        printf(" & %9.2e & %9.2e  & %9.2e \\\\ \n", err[0], err[1], err[2]);
//      }
//
//    }
//  }

  // 2-D integration of Chebyshev moments
  if (1)
  {
    cout << "---  2-D  ---" << endl;

    vector<string> qtypes;
    qtypes.push_back("u-gl");
    qtypes.push_back("u-dgl");
    qtypes.push_back("u-asdr");
    qtypes.push_back("u-gc");
    qtypes.push_back("u-dgc");
    qtypes.push_back("u-tg");

    for (int qt = 0; qt < qtypes.size(); ++qt)
    {
      cout << qtypes[qt] << " " << endl;
      for (int nn = 3; nn < 7; ++nn)
      {
        int n = nn;
        printf(" %4i ", n);
        // Define the 2-D quadrature to use
        InputDB::SP_input db(new InputDB());
        db->put<string>("quad_type", qtypes[qt]);
        db->put<int>("quad_number_polar_octant", n);
        db->put<int>("quad_number_azimuth_octant", 2);
        ProductQuadrature::SP_quadrature q = QuadratureFactory::build(db, 2);
        //q->display();
        // Define the basis set
        OrthogonalBasis::Parameters p;
        p.size = n;
        p.order = 2;
        p.x.resize(p.size, 0.0);
        p.qw.resize(p.size, 0.0);
        p.lower_bound = -1.0;
        p.upper_bound = 1.0;
        p.even_only = true;
        p.orthonormal = false;
        int i = 0;
        for (int pol = 0; pol < q->number_polar(0); ++pol, ++i)
        {
          q->incident_index(0, 0, pol);
          p.x[i]  = q->cos_theta(pol);
          p.qw[i] = 2 * q->polar_weight(pol);
        }

        OrthogonalBasis::SP_basis P = OrthogonalBasis::Create("cheby", p);

        // Isotropic flux, expansion coefficient, and approximation
        callow::Vector f(p.size, 1.0);
        callow::Vector ft(p.order + 1, 0.0);
        callow::Vector fa(p.size, 0.0);

        P->transform(f, ft);
        //ft.display("F");
        callow::Vector ref(p.order + 1, 0.0);
        callow::Vector err(ft);
        ref[0] = 1.570796326794897;
        err.subtract(ref);
        printf(" & %9.2e & %9.2e & %9.2e  \\\\ \n ", err[0], err[1], err[2]);
      }
      //printf("\\\\ \n");
    }
  }
  std::cout << std::endl << std::endl;

  // 3-D integration of Chebyshev moments
  if (1)
  {
    cout << "---  3-D  ---" << endl;
    vector<string> qtypes;
//    qtypes.push_back("u-gl");
//    qtypes.push_back("u-dgl");
//    qtypes.push_back("u-asdr");
    qtypes.push_back("u-gc");
//    qtypes.push_back("u-dgc");
//    qtypes.push_back("u-tg");

    for (int qt = 0; qt < qtypes.size(); ++qt)
    {
      cout << qtypes[qt] << endl;
      for (int nn = 2; nn < 7; ++nn)
      {
        int n = 2 * nn;
        printf("%4i ", nn);
        // Define the 2-D quadrature to use
        InputDB::SP_input db(new InputDB());
        db->put<string>("quad_type", qtypes[qt]);
        db->put<int>("quad_number_polar_octant", nn);
        db->put<int>("quad_number_azimuth_octant", 2);
        //THROW("lala");
        ProductQuadrature::SP_quadrature q = QuadratureFactory::build(db, 3);
        //q->display();
        // Define the basis set
        OrthogonalBasis::Parameters p;
        p.size = n;
        p.order = 2;
        p.x.resize(p.size, 0.0);
        p.qw.resize(p.size, 0.0);
        p.lower_bound = -1.0;
        p.upper_bound = 1.0;
        //p.even_only = true;
        p.orthonormal = false;
        int i = 0;
//         q->display();
//         q->display_indices();
        int np = q->number_polar(0) / 2;
        for (int pol = 0; pol < nn; ++pol, ++i)
        {
          q->incident_index(0, 0, pol);
          p.x[nn - pol - 1]  = -q->cos_theta(pol);
          p.x[pol + nn]      =  q->cos_theta(pol);
          p.qw[nn - pol - 1] =  q->polar_weight(pol);
          p.qw[pol + nn]     =  q->polar_weight(pol);
        }
//         for (int i = 0; i < n; ++i)
//         {
//           std::cout << i << " " << p.x[i] << " " <<  p.qw[i] << std::endl;
//         }
//         std::cout << " <--" << std::endl;

        OrthogonalBasis::SP_basis P = OrthogonalBasis::Create("cheby", p);

        // Isotropic flux, expansion coefficient, and approximation
        callow::Vector f(p.size, 1.0);
        callow::Vector ft(p.order + 1, 0.0);
        callow::Vector fa(p.size, 0.0);

        P->transform(f, ft);
        //ft.display("F");
        callow::Vector ref(p.order + 1, 0.0);
        callow::Vector err(ft);
        ref[0] = 1.570796326794897;
        err.subtract(ref);
        printf(" & %9.2e & %9.2e & %9.2e  \\\\\n ", err[0], err[1], err[2]);
      }
      //printf("\\\\ \n");
    }
  }

  callow_finalize();
  return 0;
}

//----------------------------------------------------------------------------//
//              end of test_Jacobi01.cc
//----------------------------------------------------------------------------//
