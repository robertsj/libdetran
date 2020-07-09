/*
 * test_ProjectedOperator.cc
 *
 *  Created on: Jul 3, 2020
 *      Author: rabab
 */



#define TEST_LIST                      \
        FUNC(test_ProjectedOperator_general)

#include "TestDriver.hh"
#include "callow/utils/Initialization.hh"
#include "utilities/MathUtilities.hh"
#include "callow/solver/EigenSolverCreator.hh"
#include "solvers/EigenvalueManager.hh"
#include "FixedSourceManager.hh"
#include "projection_fixture.hh"
#include "ProjectedOperator.hh"
#include "solvers/eigen/EnergyDependentEigenLHS.hh"
#include "solvers/eigen/EnergyIndependentEigenOperator.hh"


using namespace detran_test;
using namespace detran;
using namespace detran_utilities;
#include <iostream>
#include <fstream>
using namespace std;


int main(int argc, char *argv[])
{
  callow_initialize(argc, argv);
  RUN(argc, argv);
  callow_finalize();
}


int test_ProjectedOperator_general(int argc, char *argv[])
{
 typedef detran_utilities::SP<callow::MatrixDense>  SP_matrix;

 SP_matrix A;
 A = new callow::MatrixDense(5, 5);

 double row0[] = {-2.0, 1.1, 0.0, 0.0, 0.0};
 double row1[] = {1.0, -2.0, 1.1, 0.0, 0.0};
 double row2[] = {0.0, 1.0, -2.0, 1.1, 0.0};
 double row3[] = {0.0, 0.0, 1.0, -2.0, 1.1};
 double row4[] = {0.0, 0.0, 0.0, 1.0, -2.0};

 A->insert_row(0, row0, A->INSERT);
 A->insert_row(1, row1, A->INSERT);
 A->insert_row(2, row2, A->INSERT);
 A->insert_row(3, row3, A->INSERT);
 A->insert_row(4, row4, A->INSERT);
 A->assemble();

 SP_matrix U;
 U = new callow::MatrixDense(5, 2);

 double col0[] = {1, 1, 1, 1, 1};
 double col1[] = {2, 2, 2, 2, 2};

 U->insert_col(0, col0, U->INSERT);
 U->insert_col(1, col1, U->INSERT);

 detran::ProjectedOperator projector(1);

 projector.SetOperators(A, U);

 SP_matrix Ar;
 Ar = new callow::MatrixDense(2, 2);
 projector.Project(Ar);

 TEST(soft_equiv((*Ar)(0, 0), -1.6 ));
 TEST(soft_equiv((*Ar)(1, 0), -3.2 ));
 TEST(soft_equiv((*Ar)(0, 1), -3.2));
 TEST(soft_equiv((*Ar)(1, 1), -6.4));

 Ar->print_matlab("Ar.txt");


 //ProjectedOperator<callow::MatrixDense> Ar(A, U);

return 0;
}

