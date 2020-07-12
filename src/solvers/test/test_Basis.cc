/*
 * test_Basis.cc
 *
 *  Created on: Jun 15, 2020
 *      Author: rabab
 */



// LIST OF TEST FUNCTIONS
#define TEST_LIST                              \
        FUNC(test_Basis)                  \


#include "TestDriver.hh"
#include "solvers/rom/ROMBasis.hh"
#include "utilities/MathUtilities.hh"
#include "solvers/EigenvalueManager.hh"

using namespace detran_test;
using namespace detran;
using namespace detran_utilities;
using namespace std;
using std::cout;
using std::endl;

typedef callow::MatrixDense::SP_matrix             SP_matrix;

int main(int argc, char *argv[])
{
  callow_initialize(argc, argv);
  RUN(argc, argv);
  callow_finalize();
}


int test_Basis(int argc, char *argv[])
{

SP_matrix U;
U = new callow::MatrixDense(20, 5);

ROMBasis Basis(1);

Basis.GetBasis("/home/rabab/Desktop/test_basis", U);
//std::cout << Basis.rank();

std::cout << (*U)(1, 1) << "\n";

std::cout << (*U)(2, 2) << "\n";

return 0;
}


//----------------------------------------------------------------------------//
//              end of test_EigenArnoldi.cc
//----------------------------------------------------------------------------//
