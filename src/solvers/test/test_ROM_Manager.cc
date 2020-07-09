/*
 * test_ROM_Manager.cc
 *
 *  Created on: Jul 9, 2020
 *      Author: rabab
 */


#define TEST_LIST                      \
        FUNC(test_ROM_Manager)


#include "TestDriver.hh"
#include "callow/utils/Initialization.hh"
#include "utilities/MathUtilities.hh"
#include "ROM_Manager.hh"
#include "projection_fixture.hh"
#include "callow/matrix/MatrixBase.hh"
#include "callow/matrix/MatrixDense.hh"

typedef callow::MatrixDense::SP_matrix SP_matrix;
using namespace detran_test;
using namespace detran;
using namespace detran_utilities;
using namespace std;


int main(int argc, char *argv[])
{
  callow_initialize(argc, argv);
  RUN(argc, argv);
  callow_finalize();
}

int test_ROM_Manager(int argc, char *argv[])
{
 Mesh1D::SP_mesh mesh = get_mesh();
 Material::SP_material mat = get_mat();
 InputDB::SP_input input = get_input();

 ROM_Manager ROM(input, mesh, mat, "EnergyDependent");

 SP_matrix U;
 U = new callow::MatrixDense(40, 5);
 double row0[] = {1, 1, 1, 1, 1};
 for (int i=0; i<20; i++)
{
 U->insert_row(i, row0, U->INSERT);
}

ROM.Solve(U);
 return 0;
}




