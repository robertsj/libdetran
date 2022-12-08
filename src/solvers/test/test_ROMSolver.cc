/*
 * test_ROMSolver.cc
 *
 *  Created on: Jul 9, 2020
 *      Author: rabab
 */


#define TEST_LIST \
        FUNC(test_ROM_diffusion)\
		FUNC(test_ROM_EnergyIndependent)\
		FUNC(test_ROM_EnergyDependent)


#include "TestDriver.hh"
#include "callow/utils/Initialization.hh"
#include "utilities/MathUtilities.hh"
#include "projection_fixture.hh"
#include "callow/matrix/MatrixBase.hh"
#include "callow/matrix/MatrixDense.hh"
#include "solvers/EigenvalueManager.hh"
#include "callow/vector/Vector.hh"
#include "solvers/rom/ROMSolver.hh"
#include "solvers/rom/ROMBasis.hh"

using namespace detran_test;
using namespace detran;
using namespace detran_utilities;
using namespace std;

typedef callow::MatrixDense::SP_matrix SP_matrix;
typedef callow::Vector::SP_vector      SP_vector;



int main(int argc, char *argv[])
{
  callow_initialize(argc, argv);
  RUN(argc, argv);
  callow_finalize();
}

int test_ROM_diffusion(int argc, char *argv[])
{
 Mesh1D::SP_mesh mesh = get_mesh(1, "core");
 Material::SP_material mat = get_mat();
 InputDB::SP_input input = get_input();
 input->put<std::string>("operator", "diffusion");
 input->put<std::string>("equation", "diffusion");
 input->display();

 int n = mesh->number_cells();
 int r = 5;

// get the basis
 SP_matrix U;
 U  = std::make_shared<callow::MatrixDense>(2*n, 2*r);
 ROMBasis::GetBasis("../../../source/src/solvers/test/flux_basis_core0_diff", U);

 // ROM
 ROMSolver<_1D> ROM(input, mesh, mat);
 SP_vector  ROM_flux;
 ROM_flux  = std::make_shared<callow::Vector>(2*n, 0.0);
 ROM.Solve(U, ROM_flux);
 double keff_rom = ROM.keff();

 // FOM
 EigenvalueManager<_1D> manager(input, mat, mesh);
 manager.solve();
 double keff_fom = manager.state()->eigenvalue();

// error and testing
 callow::Vector phi0_fom(n, 0.0);
 callow::Vector phi1_fom(n, 0.0);

 callow::Vector phi0_rom(n, 0.0);
 callow::Vector phi1_rom(n, 0.0);

 for (int i = 0; i < n; ++i)
{
  phi0_fom[i] = manager.state()->phi(0)[i];
  phi1_fom[i] = manager.state()->phi(1)[i];
  phi0_rom[i] = (*ROM_flux)[i];
  phi1_rom[i] = (*ROM_flux)[i + n];
}
 vec_dbl error1 (n, 0);
 vec_dbl error2 (n, 0);

 for (int i = 0; i < n; ++i)
 {
   error1[i] = phi0_fom[i]/phi0_fom.norm() - phi0_rom[i]/phi0_rom.norm();

   error2[i] = phi1_fom[i]/phi1_fom.norm() - phi1_rom[i]/phi1_rom.norm();
 }

 printf(" keff_ref %10.5e  rom %10.5e\n", keff_fom, keff_rom);
 std::cout << "The error in group 1 is " << detran_utilities::norm(error1, "L2") << "\n";
 std::cout << "The error in group 2 is " << detran_utilities::norm(error2, "L2") << "\n";
 std::cout << "The absolute error in the eigenvalue  " << abs(keff_rom - keff_fom) << "\n";

 TEST(soft_equiv(manager.state()->eigenvalue(), ROM.keff(), 1E-6));

 return 0;
}


int test_ROM_EnergyIndependent(int argc, char *argv[])
{
  Mesh1D::SP_mesh mesh = get_mesh(1, "core");
  Material::SP_material mat = get_mat();
  InputDB::SP_input input = get_input();
  input->put<std::string>("equation", "dd");
  input->put<std::string>("operator", "energy-independent");

  ROMSolver<_1D> ROM(input, mesh, mat);
  int n = mesh->number_cells();
  int r = 7;

  // get the basis
  SP_matrix U;
  U  = std::make_shared<callow::MatrixDense>(n, r);
  ROMBasis::GetBasis("../../../source/src/solvers/test/fission_density_core0_transport_r=7", U);
  SP_vector  fd_rom;

  // ROM
  fd_rom  = std::make_shared<callow::Vector>(n, 0.0);
  ROM.Solve(U, fd_rom);
  double keff_rom = ROM.keff();

 //FOM
  EigenvalueManager<_1D> manager(input, mat, mesh);
  manager.solve();
  double keff_fom = manager.state()->eigenvalue();

  // error and testing
  callow::Vector fd_fom(n, 0.0);
  for (int i = 0; i < n; ++i)
  {
	fd_fom[i] = manager.fissionsource()->density()[i];
  }

  callow::Vector error(n, 0.0);

  for (int i = 0; i < n; ++i)
  {
    error[i] = fd_fom[i]/fd_fom.norm() - (*fd_rom)[i]/fd_rom->norm();
  }

  printf(" keff_ref %10.5e  rom %10.5e\n", keff_fom, keff_rom);


  std::cout << "the error norm in the fission density =   " << error.norm()<<  "\n";
  std::cout << "The absolute error in the eigenvalue = " << abs(keff_rom - keff_fom) << "\n";

  TEST(soft_equiv(manager.state()->eigenvalue(), ROM.keff(), 1E-6));

 return 0;
}

int test_ROM_EnergyDependent(int argc, char *argv[])
{
  Mesh1D::SP_mesh mesh = get_mesh(1, "core");
  Material::SP_material mat = get_mat();
  InputDB::SP_input input = get_input();


  input->put<std::string>("equation", "dd");
  input->put<std::string>("operator", "energy-dependent");

  int n = mesh->number_cells();
  int r = 7;

  // get the basis
  SP_matrix U;
  U  = std::make_shared<callow::MatrixDense>(2*n, 2*r);
  ROMBasis::GetBasis("../../../source/src/solvers/test/flux_basis_core0_transport_r=7", U);

  // ROM
  ROMSolver<_1D> ROM(input, mesh, mat);
  SP_vector  ROM_flux;
  ROM_flux  = std::make_shared<callow::Vector>(2*n, 0.0);
  ROM.Solve(U, ROM_flux);
  double keff_rom = ROM.keff();

  // FOM
  EigenvalueManager<_1D> manager(input, mat, mesh);
  manager.solve();
  double keff_fom = manager.state()->eigenvalue();

  // error and testing
  callow::Vector phi0_fom(n, 0.0);
  callow::Vector phi1_fom(n, 0.0);

  callow::Vector phi0_rom(n, 0.0);
  callow::Vector phi1_rom(n, 0.0);

  for (int i = 0; i < n; ++i)
  {
   phi0_fom[i] = manager.state()->phi(0)[i];
   phi1_fom[i] = manager.state()->phi(1)[i];
   phi0_rom[i] = (*ROM_flux)[i];
   phi1_rom[i] = (*ROM_flux)[i + n];
  }
  vec_dbl error1 (n, 0);
  vec_dbl error2 (n, 0);

  for (int i = 0; i < n; ++i)
  {
    error1[i] = phi0_fom[i]/phi0_fom.norm() - phi0_rom[i]/phi0_rom.norm();

    error2[i] = phi1_fom[i]/phi1_fom.norm() - phi1_rom[i]/phi1_rom.norm();
  }

  printf(" keff_ref %10.5e  rom %10.5e\n", keff_fom, keff_rom);


  std::cout << "The error in group 1 is " << detran_utilities::norm(error1, "L2") << "\n";
  std::cout << "The error in group 2 is " << detran_utilities::norm(error2, "L2") << "\n";
  std::cout << "The absolute error in the eigenvalue  " << abs(keff_rom - keff_fom) << "\n";


 TEST(soft_equiv(manager.state()->eigenvalue(), ROM.keff(), 1E-6))

 return 0;
}


