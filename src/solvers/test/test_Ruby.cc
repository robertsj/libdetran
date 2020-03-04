// Ruby is a nice name :)

// 0. Make db
// 1. Make materials
// 2. Make mesh
// 3. Make solver.

//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_EigenArnoldi.cc
 *  @brief Test of EigenPI
 *  @note  Copyright(C) 2020 me
 */
//----------------------------------------------------------------------------//
// LIST OF TEST FUNCTIONS
#define TEST_LIST                          \
        FUNC(test_Ruby_0)

#include "TestDriver.hh"
#include "solvers/EigenvalueManager.hh"
#include "solvers/test/eigenvalue_fixture.hh"
#include "utilities/MathUtilities.hh"

using namespace detran_test;
using namespace detran;
using namespace detran_utilities;
using namespace std;
using std::cout;
using std::endl;

int main(int argc, char *argv[]) {
	callow_initialize(argc, argv);
	RUN(argc, argv);
	callow_finalize();
}

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

int test_Ruby_0(int argc, char *argv[])
{
	int ng = 1;
    InputDB::SP_input db = InputDB::Create();
    db->put<std::string>("equation", "diffusion");
	db->put<int>("number_groups", ng);
	db->put<int>("outer_print_level", 1);
	db->put<int>("inner_print_level", 0);
	db->put<std::string>("bc_west", "reflect");
	db->put<std::string>("bc_east", "reflect");
	db->put<double>("inner_tolerance", 1e-18);
	db->put<double>("outer_tolerance", 1e-16);
	db->put<double>("eigen_tolerance", 1e-16);
	db->put<int>("inner_max_iters", 1000000);
	db->put<int>("outer_max_iters", 1000000);
	db->put<int>("eigen_max_iters", 1000000);

	Material::SP_material mat = Material::Create(1, ng);

	mat->set_sigma_t(0, 0, 1.0);
	mat ->set_sigma_f(0, 0, 0.3);
	mat->set_sigma_s(0, 0, 0, 0.5);
    mat->compute_sigma_a();
    mat->compute_diff_coef();
    mat->set_chi(0, 0, 1);
    mat->finalize();
    mat -> display();
    Mesh::SP_mesh mesh;
    vec_int fm(1, 5);
    vec_dbl cm(2, 0.0);
    cm[1] = 5.0;
    vec_int mt(1, 0);
    mesh = Mesh1D::Create(fm, cm, mt);

    EigenvalueManager<_1D> manager(db, mat, mesh);
    manager.solve();




	return 0;
}

