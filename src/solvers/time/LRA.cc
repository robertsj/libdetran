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
  //, d_T(mesh->number_cells(), 300.0)
  //, d_T_old(mesh->number_cells(), 300.0)
  , d_P(mesh->number_cells(), 0.0)
  , d_A(0.0)
  , d_current_time(0.0)
  , d_temp_0(300.0)
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
  d_physics->variable(0).resize(d_mesh->number_cells(), 310);


  for (int i = 0; i < 6; ++i)
  {
	  d_T1[i] = T1[i];
	  d_T2[i] = T2[i];
	  d_A1[i] = A1[i];
	  d_A2[i] = A2[i];
	  d_S11[i] = S11[i];
	  d_S21[i] = S21[i];
	  d_S22[i] = S22[i];
	  d_F1[i] = F1[i];
	  d_F2[i] = F2[i];
	  d_D1[i] = D1[i];
	  d_D2[i] = D2[i];
	  d_mu0[i] = mu0[i];
	  d_mu1[i] = mu1[i];
  }
  d_NU = NU;
  d_B  = B;
  d_ALPHA   = ALPHA;
  d_GAMMA   = GAMMA;
  d_KAPPA   = KAPPA;
  d_LAMBDA0 = LAMBDA0;
  d_LAMBDA1 = LAMBDA1;
  d_BETA0   = BETA0;
  d_BETA1   = BETA1;
  d_VELOCITY0 = VELOCITY0;
  d_VELOCITY1 = VELOCITY1;


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
void LRA::perturb(double T0, double SigA20, double alpha, double gamma)
{
	Require(T0 > 0.0);
	Require(SigmaA20 > 0.0);
	Require(alpha > 0.0);
	Require(gamma > 0.0);

	// Update SigmaA20, alpha, and gamma
	d_A2[ROD] = SigA20;
	d_ALPHA = alpha;
	d_GAMMA = gamma;

	// update the fast cross sections to reflect this temperature
	for (int m = 0; m < 6; ++m)
	{
		if (m != REFLECTOR) // only FUEL has feedback
		{
			d_A1[m] = d_A1[m] * (1.0 + d_GAMMA * (std::sqrt(T0) - std::sqrt(300)));
		}
	}

	// reset the initial temperature
	for (int i = 0; i < d_physics->variable(0).size(); ++i)
	{
	  d_physics->variable(0)[i] = T0;
	}
	d_temp_0 = T0;

	initialize_materials();
}



//---------------------------------------------------------------------------//
void LRA::initialize_materials()
{

  for (size_t i = 0; i < d_mesh->number_cells(); ++i)
  {
    size_t m = d_unique_mesh_map[i];

    if (d_flag) // transport
    {
      set_sigma_t(i, 0,     d_T1[m] + d_B*d_D1[m]);
      set_sigma_t(i, 1,     d_T2[m] + d_B*d_D2[m]);
      set_sigma_s(i, 0, 0,  d_S11[m] / (1.0 + d_mu0[m]));
      set_sigma_s(i, 1, 1,  d_S22[m] / (1.0 + d_mu1[m]));
    }
    else // diffusion (put removal into total; Sgg = 0)
    {
      set_sigma_t(i, 0,     d_A1[m] + d_S21[m]  + d_B*d_D1[m]);
      set_sigma_t(i, 1,     d_A2[m]             + d_B*d_D2[m]);
    }
    set_sigma_s(i, 1, 0,  d_S21[m]);
    set_sigma_a(i, 0,     d_A1[m] + d_B*d_D1[m]);
    set_sigma_a(i, 1,     d_A2[m] + d_B*d_D2[m]);
    set_sigma_f(i, 0,     d_F1[m]);
    set_sigma_f(i, 1,     d_F2[m]);
    set_nu(i, 0,          d_NU);
    set_nu(i, 1,          d_NU);
    set_chi(i, 0,         1.00000);
    set_diff_coef(i, 0,   d_D1[m]);
    set_diff_coef(i, 1,   d_D2[m]);
    // delayed chi
    set_chi_d(i, 0, 0,    1.00000);
    set_chi_d(i, 1, 0,    1.00000);
  }
  // beta
  set_beta(0,     d_BETA0);
  set_beta(1,     d_BETA1);
  // decay constants
  set_lambda(0,   d_LAMBDA0);
  set_lambda(1,   d_LAMBDA1);
  // velocities
  set_velocity(0, d_VELOCITY0);
  set_velocity(1, d_VELOCITY1);
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
  double sigma_a2 = 0.878763 * d_A2[ROD];
  if (d_t <= 2.0)
  {
	sigma_a2 = d_A2[ROD] * (1.0 - 0.0606184 * d_t);
  }
  double delta_2 = sigma_a2 - d_A2[ROD];

  vec_dbl &T = d_physics->variable(0);

  for (int i = 0; i < d_mesh->number_cells(); ++i)
  {
    size_t m = d_unique_mesh_map[i];

    // update the THERMAL cross section
    double sa = d_A2[m];
    double del = 0.0;
    if (m == ROD)
    {
      sa = sigma_a2;
      del = delta_2;
    }
    if (d_flag)
    {
      // transport
      set_sigma_t(i, 1, d_T2[m] + d_B * d_D2[m] + del);
      set_sigma_a(i, 1, sa + d_B * d_D2[m]);
    }
    else
    {
      // diffusion
      set_sigma_t(i, 1, sa + d_B * d_D2[m]);
      set_sigma_a(i, 1, sa + d_B * d_D2[m]);
    }

    // update the FAST cross section
    double sigma_a1 = d_A1[m];
    if (m != REFLECTOR) // only FUEL has feedback
      sigma_a1 = d_A1[m] * (1.0 + d_GAMMA * (std::sqrt(T[i]) - std::sqrt(300)));
    double delta_1  = sigma_a1 - d_A1[m];

    if (d_flag)
    {
      set_sigma_t(i, 0, d_T1[m] + d_B * d_D1[m] + delta_1);
      set_sigma_a(i, 0, sigma_a1 + d_B * d_D1[m]);
    }
    else
    {
      set_sigma_t(i, 0, sigma_a1 + d_B * d_D1[m] + d_S21[m]);
      set_sigma_a(i, 0, sigma_a1 + d_B * d_D1[m]         );
    }

    // chi and fission
    set_chi(i, 0, 1.0);
    set_sigma_f(i, 0, d_F1[m]);
    set_sigma_f(i, 1, d_F2[m]);

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

    d_P[i] = d_KAPPA * F;
    if (t > 0.0)
      T[i] = d_ALPHA * F;
  }
  // std::cout << " T[0]=" << T[0] <<  " F="
	//	    << sigma_f(0, 0) * phi0[0] + sigma_f(0, 1) * phi1[0]
	//		<< std::endl;
}




} // end namespace detran_user

