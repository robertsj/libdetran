//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_MGCoarseMeshPreconditioner.cc
 *  @brief Test of MGCoarseMeshPreconditioner
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                                          \
        FUNC(test_MGCoarseMeshPreconditioner_no_condense)  \
        FUNC(test_MGCoarseMeshPreconditioner_space)        \
        FUNC(test_MGCoarseMeshPreconditioner_energy)       \
        FUNC(test_MGCoarseMeshPreconditioner_space_energy)

#include "TestDriver.hh"
#include "solvers/mg/MGCoarseMeshPreconditioner.hh"
#include "callow/utils/Initialization.hh"
#include "solvers/test/fixedsource_fixture.hh"

using namespace detran_test;
using namespace detran;
using namespace detran_utilities;
using namespace callow;
using namespace std;
using std::cout;
using std::endl;
#define COUT(c) cout << c << endl;

int main(int argc, char *argv[])
{
  callow_initialize(argc, argv);
  RUN(argc, argv);
  callow_finalize();
}

/// Dummy preconditioner
class TestPC: public MGCoarseMeshPreconditioner
{
public:
  TestPC(SP_input         input,
         SP_material      material,
         SP_mesh          mesh)
  : MGCoarseMeshPreconditioner(input, material, mesh,
      SP_scattersource(0), SP_fissionsource(0), 0, false, false, "testpc")
  , d_phi(35, 1.0)
  {
    for (int g = 0; g < 7; ++g)
      for (int i = 0; i < 5; ++i)
        d_phi[5*g+i] = g+1;
  }
  void apply(Vector &b, Vector &x){}
  SP_matrix R(){return d_restrict;}
  SP_matrix P(){return d_prolong;}
  Vector &phi(){return d_phi;}
  void build(const double keff = 1.0, SP_state state = SP_state(0))
  {
    MGCoarseMeshPreconditioner::build(1.0);
  }
private:
  Vector d_phi;
};

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

/*
 *  This tests the restriction and prolongation matrices directly to ensure
 *  they each represent identity matrices of full dimension.
 */
int test_MGCoarseMeshPreconditioner_no_condense(int argc, char *argv[])
{
  FixedSourceData data = get_fixedsource_data(1, 7);

  // No condensation, so P = R = I.
  vec_int f_per_c(7, 1);
  data.input->put<int>("mgpc_coarse_mesh_level", 1);
  data.input->put<vec_int>("mgpc_fine_per_coarse_group", f_per_c);
  TestPC pc(data.input, data.material, data.mesh);
  pc.build();
  TEST(pc.R()->number_columns() == 35);
  TEST(pc.R()->number_rows() == 35);
  TEST(pc.P()->number_columns() == 35);
  TEST(pc.P()->number_rows() == 35);
  TEST(pc.R()->number_nonzeros() == 35);
  TEST(pc.P()->number_nonzeros() == 35);
  for (int i = 0; i < 35; ++i)
  {
    for (int j = 0; j < 35; ++j)
    {
      TEST(soft_equiv(i == j ? 1.0 : 0.0, (*pc.R())(i, j)));
      TEST(soft_equiv(i == j ? 1.0 : 0.0, (*pc.P())(i, j)));
    }
  }

  return 0;
}

/*
 *  This tests condensation in space (1-D).  The fine mesh has 5 cells,
 *  and we  use a level 2.  This implies the following mapping:
 *   fine:   [0 1 2 3 4]
 *   coarse: [0 0 0 1 1]
 *  The flux starts as [1 1 1 1 1 2 2 ... 7 7 7 7 7]
 *  and condenses to   [1 1 2 2 ... 7 7]
 */
int test_MGCoarseMeshPreconditioner_space(int argc, char *argv[])
{
  FixedSourceData data = get_fixedsource_data(1, 7);

  // mesh condensation, level 2 = [3, 2]
  vec_int f_per_c(7, 1);
  data.input->put<int>("mgpc_coarse_mesh_level", 2);
  data.input->put<vec_int>("mgpc_fine_per_coarse_group", f_per_c);
  TestPC pc(data.input, data.material, data.mesh);
  //pc.R()->display(true);

  pc.R()->print_matlab("R.out");
  COUT("lala");
  pc.build();
  pc.R()->print_matlab("R.out");

  COUT("lala");

  //pc.P()->display(true);
  pc.P()->print_matlab("P.out");
  TEST(pc.R()->number_columns() == 35);
  TEST(pc.R()->number_rows()    == 14);
  TEST(pc.P()->number_columns() == 14);
  TEST(pc.P()->number_rows()    == 35);
  Vector phi_c(14, 0.0);
  for (int g = 0; g < 7; ++g)
    for (int i = 0; i < 2; ++i)
      phi_c[2 * g + i] = g + 1;
  Vector tmp_c(14, 0.0);
  pc.R()->multiply(pc.phi(), tmp_c);

  for (int i = 0; i < 14; ++i)
  {
    //DISP(phi_c[i] << " " << tmp_c[i]);
    TEST(soft_equiv(phi_c[i], tmp_c[i]));
  }
  Vector tmp_f(35, 0.0);
  pc.P()->multiply(tmp_c, tmp_f);
  for (int i = 0; i < 35; ++i)
  {
    //DISP(pc.phi()[i] << " " << tmp_f[i]);
    TEST(soft_equiv(pc.phi()[i], tmp_f[i]));
  }

  return 0;
}


/*
 *  This tests condensation in space (1-D).  The fine mesh has 5 cells,
 *  and we  use a level 2.  This implies the following mapping:
 *   fine:   [0 1 2 3 4]
 *   coarse: [0 0 0 1 1]
 *  The flux starts as [1 1 1 1 1 2 2 ... 7 7 7 7 7]
 *  and condenses to   [1 1 2 2 ... 7 7]
 */
int test_MGCoarseMeshPreconditioner_energy(int argc, char *argv[])
{
  FixedSourceData data = get_fixedsource_data(1, 7);

  // Linear energy spectrum and key
  vec_dbl spectrum(7, 0.0);
  for (int i = 0; i < 7; ++i)
    spectrum[i] = (double)(i+1) * 0.084515425472852;
  std::string spectrum_key = "ALL";
  vec_int spectrum_map(data.mesh->number_cells(), 0);
  data.mesh->add_mesh_map(spectrum_key, spectrum_map);

  // Mapping with [0 1 2 | 3 4 5 6] -> [0 | 1]
  //vec_int f_per_c(2, 3); f_per_c[1] = 4;
  vec_int f_per_c(1, 7);
  data.input->put<vec_int>("mgpc_fine_per_coarse_group", f_per_c);

  // Using our own spectrum
  data.input->put<int>("mgpc_condensation_option",  5);
  data.input->put<vec_dbl>("mgpc_spectrum",         spectrum);
  data.input->put<string>("mgpc_spectrum_key",      spectrum_key);

  // No coarse meshing in space
  data.input->put<int>("mgpc_coarse_mesh_level", 1);

  TestPC pc(data.input, data.material, data.mesh);

  pc.R()->print_matlab("R.out");
  COUT("lala");
  pc.build();
  pc.R()->print_matlab("R.out");

  COUT("lala");

  //pc.P()->display(true);
  pc.P()->print_matlab("P.out");
  return 0;
  TEST(pc.R()->number_columns() == 35);
  TEST(pc.R()->number_rows()    == 14);
  TEST(pc.P()->number_columns() == 14);
  TEST(pc.P()->number_rows()    == 35);
  Vector phi_c(14, 0.0);
  for (int g = 0; g < 7; ++g)
    for (int i = 0; i < 2; ++i)
      phi_c[2 * g + i] = g + 1;
  Vector tmp_c(14, 0.0);
  pc.R()->multiply(pc.phi(), tmp_c);
//
//  for (int i = 0; i < 14; ++i)
//  {
//    //DISP(phi_c[i] << " " << tmp_c[i]);
//    TEST(soft_equiv(phi_c[i], tmp_c[i]));
//  }
//  Vector tmp_f(35, 0.0);
//  pc.P()->multiply(tmp_c, tmp_f);
//  for (int i = 0; i < 35; ++i)
//  {
//    //DISP(pc.phi()[i] << " " << tmp_f[i]);
//    TEST(soft_equiv(pc.phi()[i], tmp_f[i]));
//  }

  return 0;
}

/*
 *  This tests condensation in space (1-D) and energy.  The fine mesh has
 *  5 cells, and we  use a level 2.  This implies the following mapping:
 *   fine:   [0 1 2 3 4]
 *   coarse: [0 0 0 1 1]
 *  In energy, we go from [0 1 2 3 4 5 6] to [0 0 0 1 1 1 1].
 *  The flux starts as [1 1 1 1 1 2 2 ... 7 7 7 7 7]
 *  and condenses to   [6 6 22 22]
 *  where 6 = 1+2+3 and 22 = 4+5+6+7
 */
int test_MGCoarseMeshPreconditioner_space_energy(int argc, char *argv[])
{
  FixedSourceData data = get_fixedsource_data(1, 7);

  // energy and mesh condensation
  vec_int f_per_c(2, 3);
  f_per_c[1] = 4;
  vec_int reg_map(5, 0);
  data.mesh->add_mesh_map("REGION", reg_map);
  data.input->put<int>("mgpc_coarse_mesh_level", 2);
  data.input->put<vec_int>("mgpc_fine_per_coarse_group", f_per_c);
  vec_dbl spectrum(7, 0.0);
  for (int g = 0; g < 7; ++g)
    spectrum[g] = g + 1;
  data.input->put<vec_dbl>("mgpc_spectrum", spectrum);
  data.input->put<std::string>("mgpc_spectrum_key", "REGION");
  data.input->put<int>("mgpc_condensation_option",
      TestPC::CONDENSE_WITH_USER_SPECTRUM);
  TestPC pc(data.input, data.material, data.mesh);
  pc.build();
  TEST(pc.R()->number_columns() == 35);
  TEST(pc.R()->number_rows() == 4);
  TEST(pc.P()->number_columns() == 4);
  TEST(pc.P()->number_rows() == 35);
  double phi_c[] =
  { 6, 6, 22, 22 };
  Vector tmp_c(4, 0.0);
  pc.R()->multiply(pc.phi(), tmp_c);
  for (int i = 0; i < 4; ++i)
  {
    TEST(soft_equiv(phi_c[i], tmp_c[i]));
  }
  Vector tmp_f(35, 0.0);
  pc.P()->multiply(tmp_c, tmp_f);
  for (int i = 0; i < 35; ++i)
  {
    TEST(soft_equiv(pc.phi()[i], tmp_f[i]));
  }

  return 0;
}

//----------------------------------------------------------------------------//
//              end of test_MGCoarseMeshPreconditioner.cc
//----------------------------------------------------------------------------//
