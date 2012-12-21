//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   LRA.cc
 *  @author robertsj
 *  @date   Nov 29, 2012
 *  @brief  LRA class definition.
 */
//---------------------------------------------------------------------------//

#include "LRA.hh"

namespace detran_user
{



//---------------------------------------------------------------------------//
LRA::LRA(SP_mesh mesh, bool doingtransport, bool steady)
  : Base(mesh->number_cells(), 2, 2, "LRA_MATERIAL")
  , d_mesh(mesh)
  , d_flag(doingtransport)
  , d_steady(steady)
  , d_T(mesh->number_cells(), 300.0)
  , d_T_old(mesh->number_cells(), 300.0)
  , d_P(mesh->number_cells(), 0.0)
  , d_A(0.0)
  , d_current_time(0.0)
{
  // Preconditions
  Require(d_mesh->mesh_map_exists("COARSEMESH"));
  Require(d_mesh->mesh_map_exists("ASSEMBLY"));

  d_unique_mesh_map = d_mesh->mesh_map("COARSEMESH");
  d_assembly_map    = d_mesh->mesh_map("ASSEMBLY");
  vec_int mat_map   = d_mesh->mesh_map("MATERIAL");
  for (int i = 0; i < d_mesh->number_cells(); ++i)
  {
    // Requiring unique fine mesh materials.
    Require(mat_map[i] == i);
  }

  for (int i = 0; i < d_mesh->number_cells(); ++i)
    if (d_unique_mesh_map[i] != 4) d_A += d_mesh->volume(i);

  // Create physics
  d_physics = new detran::MultiPhysics(1);
  d_physics->variable(0).resize(d_mesh->number_cells(), 300.0);

  initialize_materials();
}

//---------------------------------------------------------------------------//
LRA::SP_material LRA::Create(SP_mesh mesh, bool flag, bool steady)
{
  SP_material p(new LRA(mesh, flag, steady));
  return p;
}

//---------------------------------------------------------------------------//
void LRA::set_state(SP_state state)
{
  Require(state);
  d_state = state;
}

//---------------------------------------------------------------------------//
void LRA::initialize_materials()
{

  for (size_t i = 0; i < d_mesh->number_cells(); ++i)
  {
    size_t m = d_unique_mesh_map[i];

    if (d_flag) // transport
    {
      set_sigma_t(i, 0,     T1[m] + B*D1[m]);
      set_sigma_t(i, 1,     T2[m] + B*D2[m]);
      set_sigma_s(i, 0, 0,  S11[m] / (1.0 + mu0[m]));
      set_sigma_s(i, 1, 1,  S22[m] / (1.0 + mu1[m]));
    }
    else // diffusion (put removal into total; Sgg = 0)
    {
      set_sigma_t(i, 0,     A1[m] + S21[m]  + B*D1[m]);
      set_sigma_t(i, 1,     A2[m]           + B*D2[m]);
    }
    set_sigma_s(i, 1, 0,  S21[m]);
    set_sigma_a(i, 0,     A1[m] + B*D1[m]);
    set_sigma_a(i, 1,     A2[m] + B*D2[m]);
    set_sigma_f(i, 0,     F1[m]);
    set_sigma_f(i, 1,     F2[m]);
    set_nu(i, 0,          NU);
    set_nu(i, 1,          NU);
    set_chi(i, 0,         1.00000);
    set_diff_coef(i, 0,   D1[m]);
    set_diff_coef(i, 1,   D2[m]);
    // delayed chi
    set_chi_d(i, 0, 0,    1.00000);
    set_chi_d(i, 1, 0,    1.00000);
  }
  // beta
  set_beta(0,     BETA0);
  set_beta(1,     BETA1);
  // decay constants
  set_lambda(0,   LAMBDA0);
  set_lambda(1,   LAMBDA1);
  // velocities
  set_velocity(0, VELOCITY0);
  set_velocity(1, VELOCITY1);
  // finalize and return
  finalize();
}

//---------------------------------------------------------------------------//
void LRA::update_impl()
{
  initialize_materials();

  if (d_steady) return;

  // Remember, d_t is the time given to the step routine.  It
  // is the time at which we compute the flux.  d_dt, however,
  // may be a half step if extrapolating.

  // Thermal cross section perturbation
  double sigma_a2 = 0.878763 * A2[ROD];
  if (d_t <= 2.0) sigma_a2 = A2[ROD] * (1.0 - 0.0606184 * d_t);
  double delta_2 = sigma_a2 - A2[ROD];

  vec_dbl &T = d_physics->variable(0);

  for (int i = 0; i < d_mesh->number_cells(); ++i)
  {
    size_t m = d_unique_mesh_map[i];

    // update the THERMAL cross section
    double sa = A2[m];
    double del = 0.0;
    if (m == ROD)
    {
      sa = sigma_a2;
      del = delta_2;
    }
    if (d_flag)
    {
      // transport
      set_sigma_t(i, 1, T2[m] + B * D2[m] + del);
      set_sigma_a(i, 1, sa + B * D2[m]);
    }
    else
    {
      // diffusion
      set_sigma_t(i, 1, sa + B * D2[m]);
      set_sigma_a(i, 1, sa + B * D2[m]);
    }

    // update the FAST cross section
    double sigma_a1 = A1[m];
    if (m != REFLECTOR) // only FUEL has feedback
      sigma_a1 = A1[m] * (1.0 + GAMMA * (std::sqrt(T[i]) - std::sqrt(300.0)));
    double delta_1  = sigma_a1 - A1[m];

    if (d_flag)
    {
      set_sigma_t(i, 0, T1[m] + B * D1[m] + delta_1);
      set_sigma_a(i, 0, sigma_a1 + B * D1[m]);
    }
    else
    {
      set_sigma_t(i, 0, sigma_a1 + B * D1[m] + S21[m]);
      set_sigma_a(i, 0, sigma_a1 + B * D1[m]         );
    }

    // chi and fission
    set_chi(i, 0, 1.0);
    set_sigma_f(i, 0, F1[m]);
    set_sigma_f(i, 1, F2[m]);

  }

}

//---------------------------------------------------------------------------//
void LRA::update_P_and_T(double t, double dt)
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
  std::cout << " T[0]=" << T[0] <<  " F=" << sigma_f(0, 0) * phi0[0] + sigma_f(0, 1) * phi1[0] << std::endl;
}




} // end namespace detran_user

