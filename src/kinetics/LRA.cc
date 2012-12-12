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

// Raw material data
//                      fuel 1, blade in    fuel 1, blade out   fuel 2, blade in    fuel 2, blade out   reflector           cr (= fuel 2, blade in @ t = 0)
const double T1[]  =  { 0.265604249667995 , 0.262881177707676 , 0.264760391845380 , 0.264760391845380 , 0.265181649429859 , 0.264760391845380 };
const double T2[]  =  { 1.579778830963665 , 1.752541184717841 , 1.594133588394707 , 1.594133588394707 , 2.093802345058626 , 1.594133588394707 };
const double A1[]  =  { 0.008252000000000 , 0.007181000000000 , 0.008002000000000 , 0.008002000000000 , 0.000602400000000 , 0.008002000000000 };
const double A2[]  =  { 0.100300000000000 , 0.070470000000000 , 0.083440000000000 , 0.073324000000000 , 0.019110000000000 , 0.083440000000000 };
const double F1[]  =  { 0.001893827160494 , 0.001896707818930 , 0.001918930041152 , 0.001918930041152 , 0.000000000000000 , 0.001918930041152 };
const double F2[]  =  { 0.044897119341564 , 0.035699588477366 , 0.042016460905350 , 0.042016460905350 , 0.000000000000000 , 0.042016460905350 };
const double S11[] =  { 0.232022249667995 , 0.228030177707676 , 0.230588391845380 , 0.230588391845380 , 0.21703924942985  , 0.230588391845380 };
const double S21[] =  { 0.025330000000000 , 0.027670000000000 , 0.026170000000000 , 0.026170000000000 , 0.047540000000000 , 0.026170000000000 };
const double S22[] =  { 1.479478830963665 , 1.682071184717841 , 1.510693588394707 , 1.520809588394707 , 2.074692345058626 , 1.510693588394707 };
const double D1[]  =  { 1.255000000000000 , 1.268000000000000 , 1.259000000000000 , 1.259000000000000 , 1.257000000000000 , 1.259000000000000 };
const double D2[]  =  { 0.211000000000000 , 0.190200000000000 , 0.209100000000000 , 0.209100000000000 , 0.159200000000000 , 0.209100000000000 };
const double mu0[] =  { 0.000000000000000 , 0.000000000000000 , 0.000000000000000 , 0.000000000000000 , 0.000000000000000 , 0.000000000000000 };
const double mu1[] =  { 0.000000000000000 , 0.000000000000000 , 0.000000000000000 , 0.000000000000000 , 0.000000000000000 , 0.000000000000000 };
const double NU = 2.43;
const double B  = 0.0001;
const double ALPHA = 3.830e-11;
const double GAMMA = 2.034e-3;
const double KAPPA = 3.204e-11;
const int ROD = 5;
const int REFLECTOR = 4;

//---------------------------------------------------------------------------//
LRA::LRA(SP_mesh mesh, bool flag, bool steady)
  : Base(mesh->number_cells(), 2, 2, "LRA_MATERIAL")
  , d_mesh(mesh)
  , d_flag(flag)
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

    if (d_flag)
    {
      set_sigma_t(i, 0,     T1[m] + B*D1[m]);
      set_sigma_t(i, 1,     T2[m] + B*D2[m]);
      set_sigma_s(i, 0, 0,  S11[m] / (1.0 + mu0[m]));
      set_sigma_s(i, 1, 1,  S22[m] / (1.0 + mu1[m]));
    }
    else
    {
      set_sigma_t(i, 0,     A1[m] + S21[m]  + B*D1[m]);
      set_sigma_t(i, 1,     A2[m] + B*D2[m]          );
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
  set_beta(0,     0.0054);
  set_beta(1,     0.0010873);
  // decay constants
  set_lambda(0,   0.00654);
  set_lambda(1,   1.35);
  // velocities
  set_velocity(0, 3.0e7);
  set_velocity(1, 3.0e5);
  // finalize and return
  finalize();
}

//---------------------------------------------------------------------------//
void LRA::update_impl()
{
  initialize_materials();

  if (d_steady) return;

  // Thermal cross section perturbation
  double sigma_a2 = 0.878763 * A2[ROD];
  if (d_t <= 2.0) sigma_a2 = A2[ROD] * (1.0 - 0.0606184 * d_t);
  double delta_2 = sigma_a2 - A2[ROD];

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
      set_sigma_t(i, 1, T2[m] + B * D2[m] + delta_2);
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
      sigma_a1 = A1[m] * (1.0 + GAMMA * (std::sqrt(d_T[i]) - std::sqrt(300.0)));
    double delta_1  = sigma_a1 - A1[m];

    if (d_flag)
    {
      set_sigma_t(i, 0, sigma_a1 + B * D1[m] + delta_1);
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

//  if (time() > 0)
//  {
//    //display();
//    std::cout << " KEFF = " << kcrit() << std::endl;
//    THROW("Done");
//  }

}

//---------------------------------------------------------------------------//
void LRA::update_P_and_T(double t)
{
  bool step = false;
  if (t > d_current_time)
  {
    d_current_time = t;
    step = true;
  }

  // Get fluxes
  const detran::State::moments_type &phi0 = d_state->phi(0);
  const detran::State::moments_type &phi1 = d_state->phi(1);

  // Compute power and temperature.  Note, we "unscale" by keff.
  for (size_t i = 0; i < d_mesh->number_cells(); ++i)
  {
    double F = sigma_f(i, 0) * phi0[i] + sigma_f(i, 1) * phi1[i];
    d_P[i] = KAPPA * F;
    if (d_t > 0.0) d_T[i] = d_T_old[i] + d_dt * ALPHA * F;
    if (step) d_T_old[i] = d_T[i];
  }


}

} // end namespace detran_user

