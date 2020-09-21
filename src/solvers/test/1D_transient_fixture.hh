/**
 *  @file   test_1D_transient.cc
 *  @author Ruby
 *  @date   March 16, 2020
 *  @brief  1D slab transient
 */
//-

#define TEST_LIST               \
        FUNC(test_1D_transient)

#include "TestDriver.hh"
#include "Mesh1D.hh"
#include "external_source/ConstantSource.hh"
#include "callow/utils/Initialization.hh"
#include "callow/vector/Vector.hh"
#include "boundary/BoundaryTraits.hh"
#include "boundary/BoundarySN.hh"
#include "kinetics/LinearExternalSource.hh"
#include "kinetics/LinearMaterial.hh"
#include "solvers/EigenvalueManager.hh"
#include "kinetics/TimeDependentMaterial.hh"
#include "TimeStepper.hh"


using namespace detran_test;
using namespace detran;
using namespace detran_material;
using namespace detran_external_source;
using namespace detran_geometry;
using namespace detran_utilities;
using namespace std;
using std::cout;
using std::endl;



class SlabMaterial: public TimeDependentMaterial
{
public:
  typedef TimeDependentMaterial Base;
  SlabMaterial(bool transport = false)
    : Base(4, 2, 8, "SlabMaterial")
    , d_transport(transport)
  {
   std::cout << "SlabMaterial **********" << "\n";

    update_impl();
  }
  // Update the materials.

  void update_impl()
  {
  std::cout << "update_impl ********" << "\n";

   double t = time();
   // materials are reflector, fuel, control-1, control-2
   double xs_t[2][2] = {{0.222222,  0.666667}, {0.25641,  0.66667}};
   // double xs_a[4][2] = {{0.0002, 0.01}, {0.0105, 0.1140}, {0.0105, 0.1640}, {0.0105, 0.1640}};
   double xs_f[4][2] = {{0.0, 0.0}, {0.0012, 0.076}, {0.0012, 0.076}, {0.0012, 0.076}};
   double xs_s[4][2][2] = {{{0.190022, 0.0},
	  	      	    	   {0.032, 0.6566667}},
	  	      			   {{0.223910, 0.0},
	  	      			   {0.022, 0.55267}},
	  	      			   {{0.223910, 0.0},
	  	      				{0.022, 0.55267}},
	  						{{0.223910, 0.0},
	  						{0.022, 0.55267}}};


   double diff_coeff[4][2] = {{1.5, 0.5}, {1.3, 0.5}, {1.3, 0.5}, {1.3, 0.5}};

   for (int m = 0; m < 4; ++m)
   {
     for (int g = 0; g < 2; ++g)
     {
	  if (g == 0 & m > 0) set_chi(m, g, 1.0);
	  if (m == 0 || m ==1 )
	  {
	   set_sigma_t(m, g, xs_t[m][g]);
	  }
	  set_nu(m, g, 2.5);
	  set_sigma_f(m, g, xs_f[m][g]);
	  set_diff_coef(m, g, diff_coeff[m][g]);

	  for (int gp = 0; gp < 2; ++gp)
	  {
	   set_sigma_s(m, gp, g, xs_s[m][gp][g]);
   	  }
     }
	};

	// set absorption cross section in control rod cells
	// this is the one that moves
	set_sigma_t(2, 0, 0.256410);
	set_sigma_t(2, 1,  xa_perturbed(t));

	// this one stays in this position
	set_sigma_t(3, 0, 0.256410);
	set_sigma_t(3, 1, compute_xs_t(0.25));

	// kinetics
	set_velocity(0,      2.0e7);
	set_velocity(1,      220000);
	double beta[8] = {2.18e-4, 1.02e-03, 6.05e-4, 1.31e-03, 2.20e-03, 6.00e-4, 5.40e-04, 1.52e-04};
	double lambda[8] = {1.246700e-02, 2.829299e-02, 4.2524000e-02, 1.330420e-01, 2.924670e-01, 6.664888e-01, 1.634781e+00, 3.554601e+00};

	for (int i=0; i < 8; i++)
	{ set_lambda(i, lambda[i]);
	  set_beta(i, beta[i]);

	  for (int m = 1; m < 4; m++)
	  {
	   set_chi_d(m, i, 0, 1);
	  }
	}

	compute_sigma_a();

	finalize();
    }

  float compute_xs_t(float pos)
  {
    float pos1 = 0.0;  // rod all in
    float pos2 = 1.0;  // rod all out
    float xs1 = 0.164; // absorption if all in
    float xs2 = 0.114; // absorption if all out

    // linear interpolation
    float xs_a = -(pos1 - pos)*(xs1 - xs2)/(pos1 - pos2) + xs1;
    // add scattering cross section to compute total
    return xs_a + 0.55267;
  }

  double xa_perturbed(double time)
  {
    float xs_a;
    float initial_pos = 0.25;
    float initial_xa = compute_xs_t(initial_pos);
    double pos; //control position
    float c = 0.05/2; // withdrawal/insertion  rate


    if (2.0 <= time && time < 4.0)
    {
     pos = initial_pos + c*(time - 2);
     xs_a = compute_xs_t(pos);
    }
    else if (4 <= time  && time < 10)
    {
      pos = initial_pos + 0.05;
      xs_a = compute_xs_t(pos);
    }

   else if ( 10 <= time && time <= 12)
   {
    pos = initial_pos + 0.05 - c*(time - 10);
    xs_a = compute_xs_t(pos);
   }
   else
   {
     xs_a = initial_xa;
   }

      return xs_a;
  }

private:
  /// Flag for transport
  bool d_transport;
};

// -------------------------------------------------
Mesh1D::SP_mesh get_mesh(Mesh1D::size_t fmm = 50)
{
    Mesh1D::vec_dbl cm(8);
    cm[0] =  0.0;
    cm[1] = 10.0;
    cm[2] = 20.0;
    cm[3] = 30.0;
    cm[4] = 40.0;
    cm[5] = 50.0;
    cm[6] = 60.0;
    cm[7] = 70.0;
    Mesh1D::vec_int fm(7);
    fm[0] = fmm;
    fm[1] = fmm;
    fm[2] = fmm;
    fm[3] = fmm;
    fm[4] = fmm;
    fm[5] = fmm;
    fm[6] = fmm;

    Mesh1D::vec_int mt(7);
    mt[0] = 0; mt[1] = 1; mt[2] = 2;
    mt[3] = 1; mt[4] = 3; mt[5] = 1;
    mt[6] = 0;
    Mesh1D::SP_mesh mesh = Mesh1D::Create(fm, cm, mt);

    return mesh;
}

InputDB::SP_input get_input()
{
  InputDB::SP_input inp(new InputDB("1D slab"));
  inp->put<int>("dimension",                      1);
  inp->put<int>("number_groups",                  2);
  inp->put<std::string>("equation",               "diffusion");
  inp->put<std::string>("bc_west",                "vacuum");
  inp->put<std::string>("bc_east",                "vacuum");
  inp->put<int>("bc_zero_flux",                   0);
  inp->put<double>("ts_final_time",               62);
  inp->put<double>("ts_step_size",                0.1);
  inp->put<int>("ts_max_steps",                   1000);
  inp->put<int>("ts_output",                      0);
  inp->put<int>("ts_monitor_level",               1);
  inp->put<int>("ts_no_extrapolation",            0);
  inp->put<string>("eigen_solver",                "arnoldi");
  inp->put<string>("outer_solver",          "GMRES");
  inp->put<int>("outer_krylov_group_cutoff",      0);
  inp->put<double>("eigen_solver_tol", 1e-16);
  inp->put<int>("eigen_solver_maxit",  1000);


  //inp->put<int>("compute_boundary_flux",          1);

  inp->put<string>("inner_solver",          "SI");
  inp->put<double>("outer_tolerance",       1e-14);
  inp->put<double>("inner_tolerance",       1e-14);
  inp->put<int>("inner_max_iters",          1e5);

  inp->put<int>("inner_print_level",        0);
  inp->put<int>("outer_print_level",        1);
  // inner gmres parameters
  InputDB::SP_input db(new InputDB("inner_solver_db"));
  db->put<double>("linear_solver_atol",                 1e-16);
  db->put<double>("linear_solver_rtol",                 1e-15);
  db->put<string>("linear_solver_type",                 "jacobi");
  db->put<string>("pc_type",                            "petsc_pc");
  db->put<string>("petsc_pc_type",                      "lu");
  db->put<int>("linear_solver_maxit",                   2000);
  db->put<int>("linear_solver_gmres_restart",           30);
  db->put<int>("linear_solver_monitor_level",           0);
  db->put<string>("eigen_solver_type",                  "slepc");
  db->put<double>("eigen_solver_tol",                   1e-15);
  inp->put<InputDB::SP_input>("inner_solver_db",        db);
  inp->put<InputDB::SP_input>("outer_solver_db",        db);
  inp->put<InputDB::SP_input>("eigen_solver_db",        db);
  inp->put<int>("compute_boundary_flux",                1);
  inp->put<double>("ts_final_time",               60);
  inp->put<double>("ts_step_size",                0.1);
  inp->put<int>("linear_solver_maxit",   1000);
  inp->put<int>("linear_solver_monitor_level", 0);
  inp->put<int>("linear_solver_monitor_diverge", 0);
  if (inp->get<std::string>("equation") != "diffusion")
  {
	inp->put<int>("ts_discrete",              1);
	inp->put<int>("store_angular_flux",       1);
  }

return inp;
}


///////////////////////////////////////////////////////////////

/*
int test_1D_transient(int argc, char *argv[])
{

  typedef TimeStepper<_1D> TS_1D;

  //-------------------------------------------------------------------------//
  // INPUT
  //-------------------------------------------------------------------------//

  InputDB::SP_input inp(new InputDB("1D slab"));
  inp->put<int>("dimension",                      1);
  inp->put<int>("number_groups",                  2);
  inp->put<std::string>("equation",               "diffusion");
  inp->put<std::string>("bc_west",                "vacuum");
  inp->put<std::string>("bc_east",                "vacuum");
  inp->put<int>("bc_zero_flux",                   0);
  inp->put<double>("ts_final_time",               62);
  inp->put<double>("ts_step_size",                0.1);
  inp->put<int>("ts_max_steps",                   1000);
  inp->put<int>("ts_scheme",                      TS_1D::IMP);
  inp->put<int>("ts_output",                      0);
  inp->put<int>("ts_monitor_level",               1);
  inp->put<int>("ts_no_extrapolation",            0);
  inp->put<string>("eigen_solver",                "arnoldi");
  inp->put<string>("outer_solver",          "GMRES");
  inp->put<int>("outer_krylov_group_cutoff",      0);
  //inp->put<int>("compute_boundary_flux",          1);

  inp->put<string>("inner_solver",          "SI");
  inp->put<double>("outer_tolerance",       1e-14);
  inp->put<double>("inner_tolerance",       1e-14);
  inp->put<int>("inner_max_iters",          1e5);

  inp->put<int>("inner_print_level",        0);
  inp->put<int>("outer_print_level",        1);
  // inner gmres parameters
  InputDB::SP_input db(new InputDB("inner_solver_db"));
  db->put<double>("linear_solver_atol",                 1e-12);
  db->put<double>("linear_solver_rtol",                 1e-15);
  db->put<string>("linear_solver_type",                 "petsc");
  db->put<string>("pc_type",                            "petsc_pc");
  db->put<string>("petsc_pc_type",                      "lu");
  db->put<int>("linear_solver_maxit",                   2000);
  db->put<int>("linear_solver_gmres_restart",           30);
  db->put<int>("linear_solver_monitor_level",           0);
  db->put<string>("eigen_solver_type",                  "slepc");
  db->put<double>("eigen_solver_tol",                   1e-15);
  inp->put<InputDB::SP_input>("inner_solver_db",        db);
  inp->put<InputDB::SP_input>("outer_solver_db",        db);
  inp->put<InputDB::SP_input>("eigen_solver_db",        db);
  inp->put<int>("compute_boundary_flux",                1);
  if (inp->get<std::string>("equation") != "diffusion")
  {
    inp->put<int>("ts_discrete",              1);
    inp->put<int>("store_angular_flux",       1);
  }

  //-------------------------------------------------------------------------//
  // MATERIAL
  //-------------------------------------------------------------------------//

  // Create slab material
  bool transport = false;
  if (inp->get<std::string>("equation") != "diffusion") transport = true;
  TS_1D::SP_material mat(new SlabMaterial(transport));

  //-------------------------------------------------------------------------//
  // MESH
  //-------------------------------------------------------------------------//

  TS_1D::SP_mesh mesh = get_mesh(10);

  //-------------------------------------------------------------------------//
  // STEADY STATE
  //-------------------------------------------------------------------------//

  EigenvalueManager<_1D> manager(inp, mat, mesh);
  manager.solve();
  State::SP_state ic = manager.state();
  mat->set_eigenvalue(ic->eigenvalue());
  mat->update(0, 0, 1, false);
  mat->display();
  mesh->display();
  //ic->display();

  // Normalize state.
  double F = 0;
  TS_1D::vec_int matmap = mesh->mesh_map("MATERIAL");
  for (int i = 0; i < mesh->number_cells(); ++i)
  {
    int m = matmap[i];
    F += ic->phi(0)[i] * mat->sigma_f(m, 0) +
         ic->phi(1)[i] * mat->sigma_f(m, 1);
  }
  ic->scale(1.0/F);

  //-------------------------------------------------------------------------//
  // TIME STEPPER
  //-------------------------------------------------------------------------//


  TS_1D stepper(inp, mat, mesh, true);
  stepper.set_monitor(test_monitor);
  stepper.solve(ic);


  printf(" %20.16f %20.16f ", ic->phi(0)[0], ic->phi(0)[1]);
  std::cout << std::endl;

  State::SP_state final = stepper.state();

  printf(" %20.16f %20.16f ", final->phi(0)[0], final->phi(0)[1]);

  std::cout << std::endl;
  std::cout << "*************"<< std::endl;

  return 0;

}

*/
