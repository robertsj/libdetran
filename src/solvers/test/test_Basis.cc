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

int main(int argc, char *argv[])
{
  callow_initialize(argc, argv);
  RUN(argc, argv);
  callow_finalize();
}


int test_Basis(int argc, char *argv[])
{
	std::string a = "n";
  ROMBasis U(20, 4, 99);

return 0;
}


//----------------------------------------------------------------------------//
//              end of test_EigenArnoldi.cc
//----------------------------------------------------------------------------//
