/**
 *  @file   SlabReactor_fixture.hh
 *  @author Rabab Elzohery
 *  @date   Oct 18, 2021
 *  @brief  1D slab transient
 */
//-

#define TEST_LIST               \
        FUNC(test_1D_transient)

#include "Mesh1D.hh"
#include "TestDriver.hh"
#include "external_source/ConstantSource.hh"
#include "callow/utils/Initialization.hh"
#include "callow/vector/Vector.hh"
#include "boundary/BoundaryTraits.hh"
#include "boundary/BoundarySN.hh"
#include "kinetics/LinearExternalSource.hh"
#include "kinetics/LinearMaterial.hh"
#include "solvers/EigenvalueManager.hh"
#include "kinetics/TimeDependentMaterial.hh"
#include "kinetics/MultiPhysics.hh"
#include "TimeStepper.hh"


using namespace detran;
using namespace detran_material;
using namespace detran_external_source;
using namespace detran_geometry;
using namespace detran_utilities;
using namespace std;
using std::cout;
using std::endl;



// materials are reflector, fuel, control-1, control-2
double xs_t[2][2] = {{0.222222,  0.666667}, {0.25641,  0.66667}};
double xs_a[4][2] = {{0.0002, 0.01}, {0.0105, 0.1140}, {0.0105, 0.1640}, {0.0105, 0.1640}};
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
double beta_[8] = {2.18e-4, 1.02e-03, 6.05e-4, 1.31e-03, 2.20e-03, 6.00e-4, 5.40e-04, 1.52e-04};
double lambda_[8] = {1.246700e-02, 2.829299e-02, 4.2524000e-02, 1.330420e-01, 2.924670e-01, 6.664888e-01, 1.634781e+00, 3.554601e+00};

const double GAMMA   = 3.034e-3;
const double ALPHA   = 3.830e-11;
const double KAPPA   = 3.204e-11;
const double B = 0.001;

//---------------------------------------------------------------------------//

Mesh1D::SP_mesh get_mesh(Mesh1D::size_t fmm = 50)
{
  Mesh1D::vec_dbl cm(8);

  for (int i=0; i<8; i++) cm[i] = i*10;

  Mesh1D::vec_int fm(7, fmm);

  Mesh1D::vec_int mt(7);
  mt[0] = 0; mt[1] = 1; mt[2] = 2;
  mt[3] = 1; mt[4] = 3; mt[5] = 1;
  mt[6] = 0;

  Mesh1D::vec_int assembly(7);
  assembly[0] = 0; assembly[1] = 1; assembly[2] = 2;
  assembly[3] = 3; assembly[4] = 4; assembly[5] = 5;
  assembly[6] = 0;

  Mesh1D::SP_mesh mesh = Mesh1D::Create(fm, cm, mt);
  mesh->add_coarse_mesh_map("COARSEMESH", mt);
  mesh->add_coarse_mesh_map("ASSEMBLY", assembly);

  // fine mesh material map
  Mesh1D::vec_int  mm(mesh->number_cells(), 0);
  for (int i = 0; i < mm.size(); ++i)
    mm[i]  = i;
  mesh->add_mesh_map("MATERIAL", mm);

  return mesh;
}

//---------------------------------------------------------------------------//

InputDB::SP_input get_input()
{
  InputDB::SP_input inp(new InputDB("1D slab"));
  inp->put<int>("dimension",                      1);
  inp->put<int>("number_groups",                  2);
  inp->put<std::string>("equation",               "dd");
  inp->put<std::string>("bc_west",                "vacuum");
  inp->put<std::string>("bc_east",                "vacuum");
  inp->put<int>("bc_zero_flux",                   0);
  inp->put<double>("ts_final_time",              60);
  inp->put<double>("ts_step_size",                0.01);
  inp->put<int>("ts_max_steps",                   10000);
  inp->put<int>("ts_output",                      0);
  inp->put<int>("ts_monitor_level",               1);
  inp->put<int>("ts_no_extrapolation",            0);
  inp->put<int>("ts_max_iters",             1000);
  inp->put<double>("ts_tolerance",          1.0e-8);
  inp->put<string>("eigen_solver",                "arnoldi");
  inp->put<string>("outer_solver",                 "GMRES");
  inp->put<int>("outer_krylov_group_cutoff",      0);
  inp->put<double>("eigen_solver_tol", 1e-16);
  inp->put<int>("eigen_solver_maxit",  1000);
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
  db->put<string>("linear_solver_type",                 "gmres");
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
  inp->put<InputDB::SP_input>("rom_solver_db",          db);
  inp->put<int>("compute_boundary_flux",                1);
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

//---------------------------------------------------------------------------//

class SlabMaterial: public TimeDependentMaterial
{
  public:
  typedef TimeDependentMaterial Base;
  typedef detran_geometry::Mesh::SP_mesh            SP_mesh;
  typedef detran::MultiPhysics::SP_multiphysics     SP_multiphysics;


  SlabMaterial(SP_mesh mesh, bool transport, int perturbation)
    : Base(mesh->number_cells(), 2, 8, "SlabMaterial")
    , d_transport(transport),
	  d_mesh(mesh),
	  d_T(mesh->number_cells(), 300.0),
	  //d_T_old(mesh->number_cells(), 300.0),
	  d_P(mesh->number_cells(), 0.0),
      d_perturbation(perturbation)
  {

    d_unique_mesh_map = d_mesh->mesh_map("COARSEMESH");
	d_assembly_map    = d_mesh->mesh_map("ASSEMBLY");
	vec_int mat_map   = d_mesh->mesh_map("MATERIAL");

	for (int i = 0; i < d_mesh->number_cells(); ++i)
	{
	  // Requiring unique fine mesh materials.
	  Require(mat_map[i] == i);
	}
	// Create physics
	d_physics = new detran::MultiPhysics(1);
	d_physics->variable(0).resize(d_mesh->number_cells(), 300.0);

	initialize_materials();
  }

  //---------------------------------------------------------------------------//

  SP_material Create(SP_mesh mesh, bool flag, bool steady, int perturbation)
  {
    SP_material p(new SlabMaterial(mesh, flag,  perturbation));
    return p;
  }

  //---------------------------------------------------------------------------//

  void initialize_materials()
  {
    double t = 0;

    for (size_t i = 0; i < d_mesh->number_cells(); ++i)
    {
      size_t m = d_unique_mesh_map[i];

      if (m > 0) set_chi(i, 0, 1.0);

      set_nu(i, 0, 2.5);
      set_nu(i, 1, 2.5);

      set_sigma_f(i, 0, xs_f[m][0]);
      set_sigma_f(i, 1, xs_f[m][1]);

      set_diff_coef(i, 0, diff_coeff[m][0]);
      set_diff_coef(i, 1, diff_coeff[m][1]);

      for (int gp = 0; gp < 2; ++gp)
      {
        set_sigma_s(i, gp, 0, xs_s[m][gp][0]);
        set_sigma_s(i, gp, 1, xs_s[m][gp][1]);
      }

      if (m == 0 || m ==1 )
      {
        set_sigma_t(i, 0, xs_t[m][0]);
        set_sigma_t(i, 1, xs_t[m][1]);
      }
      // set absorption cross section in control rod cells
     // this is the one that moves
      else if (m == 2)
     {
       set_sigma_t(i, 0, 0.256410);
       set_sigma_t(i, 1,  xt_perturbed(t));
     }
     // this one stays in this position
     else if (m == 3)
     {
       set_sigma_t(i, 0, 0.256410);
       set_sigma_t(i, 1, compute_xs_t(0.25));
     }

     for (int p=0; p < 8; p++)
     {
      set_chi_d(i, p, 0, 1);
     }
     compute_sigma_a();
    }

    for (int p=0; p < 8; p++)
    {
      set_lambda(p, lambda_[p]);
      set_beta(p, beta_[p]);
    }

    set_velocity(0,      2.0e7);
    set_velocity(1,      220000);

    finalize();

  }

  //---------------------------------------------------------------------------//

  // Update the materials.
  void update_impl()
  {
    initialize_materials();

    double t = time();

    vec_dbl &T = d_physics->variable(0);

    for (int i = 0; i < d_mesh->number_cells(); ++i)
    {
      size_t m =  d_unique_mesh_map[i];

      double sigma_t2 ;

      // update the THERMAL cross section
      if (m == 2) //  material of the control that moves
      {
        sigma_t2 = xt_perturbed(t);
        set_sigma_t(i, 1, sigma_t2);
        set_sigma_a(i, 1, sigma_t2-0.55267); //adjust the absorption XS
      }

     // update the FAST cross section
    double sigma_a1;
    if (m != 0)  // fuel materials
    {
      sigma_a1 = xs_a[m][0] * (1.0 + GAMMA * (std::sqrt(T[i]) - std::sqrt(300.0)));
      double delta_1  = sigma_a1 - xs_a[m][0];

     if (d_transport)
     {
       if (m ==1) set_sigma_t(i, 0, xs_t[m][0] + delta_1);
       if (m==2| m==3) set_sigma_t(i, 0, 0.256410 + delta_1);
       set_sigma_a(i, 0, sigma_a1);

     }

     else
     {
       if (m ==1) set_sigma_t(i, 0, xs_t[m][0] + delta_1);
       if (m==2| m==3) set_sigma_t(i, 0, 0.256410 + delta_1);
       set_sigma_a(i, 0, sigma_a1);
     }
   }

    // chi and fission
    set_chi(i, 0, 1.0);
    set_sigma_f(i, 0, xs_f[m][0]);
    set_sigma_f(i, 1, xs_f[m][1]);

    }
  }

  //---------------------------------------------------------------------------//

  /// compute total cross-section by interpolation
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

  //---------------------------------------------------------------------------//

  /// perturb absorption cross-section based on control movement and compute total
  double xt_perturbed(double time)
  {
    double xs_a;
    double initial_pos = 0.25; // initial control position
    double initial_xs = compute_xs_t(initial_pos);
    double pos; //control position
    double c; // withdrawal/insertion  rate

    if (d_perturbation == 2) // prompt super critical (Ye's problem)
    {
      c = 0.15/2; //control was drawn 15% out in 2 sec

      if (2.0 <= time && time <= 4.0)
      {
        pos = initial_pos + c*(time - 2);
        xs_a = compute_xs_t(pos);
      }

      else if (time > 4.0) xs_a = compute_xs_t(c*2 + initial_pos);
      else if (time < 2)  xs_a = initial_xs;
    }

    else if (d_perturbation == 1) // delayed super critical (Rabab's problem)
    {
      c = 0.05/2; //control was drawn 5% out in 2 sec
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
      xs_a = initial_xs;
    }
  }

  return xs_a;
  }

  //---------------------------------------------------------------------------//

  void update_P_and_T(double t, double dt)
  {
    // Get fluxes
    const detran::State::moments_type &phi0 = d_state->phi(0);
    const detran::State::moments_type &phi1 = d_state->phi(1);

    // Compute power and temperature.  Note, we "unscale" by keff.
    vec_dbl &T = d_physics->variable(0);
    double F = 0;
    for (size_t i = 0; i < d_mesh->number_cells(); ++i)
    {
      F = sigma_f(i, 0) * phi0[i] + sigma_f(i, 1) * phi1[i];

      d_P[i] = KAPPA * F;
      if (t > 0.0)
        T[i] = ALPHA * F;
    }
  }


  //---------------------------------------------------------------------------//


  vec_dbl T() {return d_T;}
  vec_dbl P() {return d_P;}
  SP_multiphysics physics() {return d_physics;}

  //---------------------------------------------------------------------------//

  /// Multiphysics (just T)
   SP_multiphysics d_physics;
   /// Fine mesh temperature
   vec_dbl d_T;
   /// Old fine mesh temperature
   vec_dbl d_T_old;
   /// Fine mesh power density
   vec_dbl d_P;

   vec_int d_unique_mesh_map;

   vec_int d_assembly_map;

   int d_perturbation;

private:
  /// Flag for transport
  bool d_transport;

  SP_mesh d_mesh;


}; // end of the material class


//---------------------------------------------------------------------------//

void update_T_rhs(void* data,
                    detran::TimeStepper<_1D>* step,
                    double t,
                    double dt)
  {
    Require(data);
    Require(step);

    // cast data as LRA
    SlabMaterial* mat = (SlabMaterial*) data;

    // update
    mat->update_P_and_T(t, dt);
  }

// -------------------------------------------------
