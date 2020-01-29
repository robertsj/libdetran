//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   test_TWIGL.cc
 *  @author robertsj
 *  @date   Apr 4, 2012
 *  @brief  Implements the TWIGL kinetics benchmark
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST               \
        FUNC(test_TWIGL)

#include "TestDriver.hh"
#include "TimeStepper.hh"
#include "Mesh2D.hh"
#include "external_source/ConstantSource.hh"
#include "callow/utils/Initialization.hh"
#include "callow/vector/Vector.hh"
#include "boundary/BoundaryTraits.hh"
#include "boundary/BoundarySN.hh"
#include "kinetics/LinearExternalSource.hh"
#include "kinetics/LinearMaterial.hh"
#include "solvers/EigenvalueManager.hh"


using namespace detran_test;
using namespace detran;
using namespace detran_material;
using namespace detran_external_source;
using namespace detran_geometry;
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

//---------------------------------------------------------------------------//
// TEST DEFINITIONS
//---------------------------------------------------------------------------//

void test_monitor(void* data, TimeStepper<_2D>* ts, int step, double t,
                  double dt, int it, bool conv)
{
  double F = 0;
  TimeStepper<_2D>::vec_int matmap = ts->mesh()->mesh_map("MATERIAL");
  for (int i = 0; i < ts->mesh()->number_cells(); ++i)
  {
    int m = matmap[i];
    F += ts->state()->phi(0)[i] * ts->material()->sigma_f(m, 0) +
         ts->state()->phi(1)[i] * ts->material()->sigma_f(m, 1);
  }
  //ts->state()->display();
  printf(" %5i  %16.13f  %16.13f %5i \n", step, t, F, it);
}

//---------------------------------------------------------------------------//
class TWIGLMaterial: public TimeDependentMaterial
{
public:
  typedef TimeDependentMaterial Base;
  TWIGLMaterial(int perturbation = 0, bool transport = false)
    : Base(3, 2, 1, "TWIGLMaterial")
    , d_perturbation(perturbation)
    , d_transport(transport)
  {
    update_impl();
  }
  // Update the materials.
  void update_impl()
  {
    double t = time();

    if (d_transport)
    {
      double s_a[3][2] = {{0.01, 0.15}, {0.01, 0.15}, {0.008, 0.05}};
      double s_t[3][2] = {{0.2481, 0.9833}, {0.2481, 0.9833}, {0.2644, 0.7167}};
      double s_f[3][2] = {{0.007, 0.2}, {0.007, 0.2}, {0.003, 0.06}};
      double s_s[3][2][2] =  {  {{0.2281, 0.0000},
    	    	                 {0.0100, 0.8333}},
    			        {{0.2281, 0.0000},
    				 {0.0100, 0.8333}},
    				{{0.2464, 0.0000},
    				 {0.0100, 0.6667}} };
      for (int m = 0; m < 3; ++m)
      {
    	set_chi(m, 0, 1.0);
    	set_chi(m, 1, 0.0);
    	set_chi_d(m, 0, 0, 1.0);
    	set_chi_d(m, 0, 1, 0.0);
    	for (int g = 0; g < 2; ++g)
    	{
    	  set_sigma_t(m, g, s_t[m][g]);
    	  set_sigma_f(m, g, s_f[m][g]);
    	  set_nu(m, g, 1.0);
    	  for (int gp = 0; gp < 2; ++gp)
    	  {
    	    set_sigma_s(m, g, gp, s_s[m][g][gp]);
    	  }
    	}
      }
      set_lambda(0,        0.08);
      set_beta(0,          0.0075);
      set_velocity(0,      1.0e7);
      set_velocity(1,      2.0e5);
      if (d_perturbation == 1)
      {
        double st;
        if (t < 0.1)
          st = s_t[0][1];
        else if (t >= 0.1 && t < 0.3)
          st = s_t[0][1]*(0.3-t)/0.2 + (s_t[0][1]-0.0035)*(t-0.1)/0.2;
        else
          st = s_t[0][1]-0.0035;
        set_sigma_t(0, 1, st);
      }
      else
      {
        if (t >= 0.1)
          set_sigma_t(0, 1, s_t[0][1]-0.0035);
      }
      compute_sigma_a();
      compute_diff_coef();
      finalize();
    }
    else
    {
        int m = 0;
        set_diff_coef(m, 0,  1.400);
        set_diff_coef(m, 1,  0.400);
        double sa = 0.150;
        if (d_perturbation == 1) // RAMP
        {
          if (t >= 0.1 && t < 0.3)
            sa = 0.150 * (1.0 - 0.11667 * (t - 0.1));
          else if (t >= 0.3)
            sa = 0.150 * 0.97666;
        }
        else if (d_perturbation == 2) // STEP
        {
          if (t >= 0.1)
            sa = 0.150 - 0.0035;
        }
        if (d_transport)
        {
          set_sigma_t(m, 0,    1/(3*1.4));
          set_sigma_t(m, 1,    1/(3*0.4));
          set_sigma_s(m, 0, 0, 1/(3*1.4) - 0.01 - 0.01);
          set_sigma_s(m, 1, 1, 1/(3*0.4) - sa);
        }
        else
        {
          set_sigma_t(m, 0,    0.010 + 0.010);
          set_sigma_t(m, 1,    sa);
        }
        set_sigma_s(m, 1, 0, 0.010);
        set_sigma_f(m, 0,    0.007);
        set_sigma_f(m, 1,    0.200);
        set_nu(m, 0,         1.0);
        set_nu(m, 1,         1.0);
        set_chi(m, 0,        1.000);
        //
        m = 1;
        set_diff_coef(m, 0,  1.400);
        set_diff_coef(m, 1,  0.400);
        if (d_transport)
        {
          set_sigma_t(m, 0,    1/(3*1.4));
          set_sigma_t(m, 1,    1/(3*0.4));
          set_sigma_s(m, 0, 0, 1/(3*1.4) - 0.01 - 0.01);
          set_sigma_s(m, 1, 1, 1/(3*0.4) - 0.150);
        }
        else
        {
          set_sigma_t(m, 0,    0.010 + 0.010);
          set_sigma_t(m, 1,    0.150);
        }
        set_sigma_s(m, 1, 0, 0.010);
        set_sigma_f(m, 0,    0.007);
        set_sigma_f(m, 1,    0.200);
        set_nu(m, 0,         1.0);
        set_nu(m, 1,         1.0);
        set_chi(m, 0,        1.000);
        //
        m = 2;
        set_diff_coef(m, 0,  1.300);
        set_diff_coef(m, 1,  0.500);
        if (d_transport)
        {
          set_sigma_t(m, 0,    1/(3*1.3));
          set_sigma_t(m, 1,    1/(3*0.5));
          set_sigma_s(m, 0, 0, 1/(3*1.3) - 0.008 - 0.01);
          set_sigma_s(m, 1, 1, 1/(3*0.5) - 0.050);
        }
        else
        {
          set_sigma_t(m, 0,    0.008 + 0.01);
          set_sigma_t(m, 1,    0.050);
        }
        set_sigma_s(m, 1, 0, 0.010);
        set_sigma_f(m, 0,    0.003);
        set_sigma_f(m, 1,    0.060);
        set_nu(m, 0,         1.0);
        set_nu(m, 1,         1.0);
        set_chi(m, 0,        1.000);
        // kinetics
        set_lambda(0,        0.08);
        set_beta(0,          0.0075);
        set_velocity(0,      1.0e7);
        set_velocity(1,      2.0e5);
        // delayed chi(m, i, g, v)
        set_chi_d(0, 0, 0,    1.00000);
        set_chi_d(1, 0, 0,    1.00000);
        set_chi_d(2, 0, 0,    1.00000);
        //
        finalize();
    }

  }

private:
  /// 0-steady, 1-ramp, 2-step
  int d_perturbation;
  /// Flag for transport
  bool d_transport;
};

//---------------------------------------------------------------------------//
Mesh2D::SP_mesh get_mesh(Mesh2D::size_t fmm = 1)
{
  if (1)
  {
    Mesh2D::vec_dbl cm(4);
    cm[0] =  0.0;
    cm[1] = 24.0;
    cm[2] = 56.0;
    cm[3] = 80.0;
    Mesh2D::vec_int fm(3);
    fm[0] = fmm * 3;
    fm[1] = fmm * 4;
    fm[2] = fmm * 3;
    Mesh2D::vec_int mt(9);
    mt[0] = 2; mt[1] = 1; mt[2] = 2;
    mt[3] = 1; mt[4] = 0; mt[5] = 2;
    mt[6] = 2; mt[7] = 2; mt[8] = 2;
    Mesh2D::SP_mesh mesh = Mesh2D::Create(fm, fm, cm, cm, mt);
    return mesh;
  }
  else
  { // debugging homogeneous problem
    Mesh2D::vec_dbl cm(2);
    cm[1] = 5.0;
    Mesh2D::vec_int fm(1);
    fm[0] = 3;
    Mesh2D::vec_int mt(1);
    mt[0] = 0;
    Mesh2D::SP_mesh mesh = Mesh2D::Create(fm, fm, cm, cm, mt);
    return mesh;
  }
}

//---------------------------------------------------------------------------//
int test_TWIGL(int argc, char *argv[])
{

  typedef TimeStepper<_2D> TS_2D;

  //-------------------------------------------------------------------------//
  // INPUT
  //-------------------------------------------------------------------------//

  InputDB::SP_input inp(new InputDB("TWIGL benchmark"));
  inp->put<int>("dimension",                      2);
  inp->put<int>("number_groups",                  2);
  inp->put<std::string>("equation",               "dd");
  inp->put<std::string>("bc_west",                "reflect");
  inp->put<std::string>("bc_east",                "vacuum");
  inp->put<std::string>("bc_south",               "reflect");
  inp->put<std::string>("bc_north",               "vacuum");
  inp->put<int>("bc_zero_flux",                   0);
  inp->put<double>("ts_final_time",               0.6);
  inp->put<double>("ts_step_size",                0.001);
  inp->put<int>("ts_max_steps",                   1000);
  inp->put<int>("ts_scheme",                      TS_2D::IMP);
  inp->put<int>("ts_output",                      0);
  inp->put<int>("ts_monitor_level",               1);
  inp->put<int>("ts_no_extrapolation",            0);
  inp->put<string>("eigen_solver",                "arnoldi");
  inp->put<int>("quad_number_polar_octant",       8);
  inp->put<int>("quad_number_azimuth_octant",     8);

  inp->put<string>("outer_solver",          "GMRES");
  inp->put<int>("outer_krylov_group_cutoff",      0);
  //inp->put<int>("compute_boundary_flux",          1);

  inp->put<string>("inner_solver",          "SI");
  inp->put<double>("outer_tolerance",       1e-14);
  inp->put<double>("inner_tolerance",       1e-14);
  inp->put<int>("inner_max_iters",          1e5);

  inp->put<int>("inner_print_level",        0);
  inp->put<int>("outer_print_level",        1);
  inp->put<int>("quad_number_azimuth_octant",   1);
  inp->put<int>("quad_number_polar_octant",     1);

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

  // Create a TWIGL material with ramp reactivity
  bool transport = false;
  if (inp->get<std::string>("equation") != "diffusion") transport = true;
  TS_2D::SP_material mat(new TWIGLMaterial(2, transport));

  //-------------------------------------------------------------------------//
  // MESH
  //-------------------------------------------------------------------------//

  TS_2D::SP_mesh mesh = get_mesh(2);

  //-------------------------------------------------------------------------//
  // STEADY STATE
  //-------------------------------------------------------------------------//

  EigenvalueManager<_2D> manager(inp, mat, mesh);
  manager.solve();
  State::SP_state ic = manager.state();
  mat->set_eigenvalue(ic->eigenvalue());
  mat->update(0, 0, 1, false);

  //ic->display();
  //return 0;
  //inp->put<string>("inner_solver",          "GMRES");

  // Normalize state.
  double F = 0;
  TS_2D::vec_int matmap = mesh->mesh_map("MATERIAL");
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


  TS_2D stepper(inp, mat, mesh, true);
  stepper.set_monitor(test_monitor);
  stepper.solve(ic);


  printf(" %20.16f %20.16f ", ic->phi(0)[0], ic->phi(0)[1]);
  std::cout << std::endl;

  State::SP_state final = stepper.state();

  printf(" %20.16f %20.16f ", final->phi(0)[0], final->phi(0)[1]);
  std::cout << std::endl;

  return 0;

}
//---------------------------------------------------------------------------//
//              end of test_TWIGL.cc
//---------------------------------------------------------------------------//
