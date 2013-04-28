//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   TimeDependentMaterial.cc
 *  @brief  TimeDependentMaterial
 *  @author Jeremy Roberts
 *  @date   Nov 14, 2012
 */
//---------------------------------------------------------------------------//

#include "TimeDependentMaterial.hh"
#include "BDFCoefficients.hh"

namespace detran
{

//---------------------------------------------------------------------------//
TimeDependentMaterial::
TimeDependentMaterial(const size_t number_materials,
                      const size_t number_energy_groups,
                      const size_t number_precursor_groups,
                      std::string  name)
  : Base(number_materials, number_energy_groups, number_precursor_groups, name)
  , d_t(0.0)
  , d_dt(0.0)
  , d_kcrit(1.0)
{
  /* ... */
}

//---------------------------------------------------------------------------//
void TimeDependentMaterial::update(const double t,
                                   const double dt,
                                   const size_t order,
                                   const bool   synthetic)
{
  // Store time step data
  d_t     = t;
  d_dt    = dt;
  d_order = order;

  // Update the base materials.
  update_impl();

  // BDF coefficient
  Assert(order >= 1 && order <= 6);
  double a_0 = bdf_coefs[order - 1][0];

  // Add synthetic components
  if (synthetic)
  {
    for (size_t m = 0; m < number_materials(); ++m)
    {
      for (size_t g = 0; g < number_groups(); ++g)
      {
        // synthetic total cross section
        d_sigma_t[g][m] += a_0 / (d_velocity[g] * dt);

        // synthetic prompt chi spectrum
        double chi = (1.0 - beta_total(m)) * d_chi[g][m];
        for (size_t i = 0; i < d_number_precursor_groups; ++i)
        {
          double den = a_0 + dt * d_lambda[i];
          chi += d_lambda[i] * d_chi_d[g][i][m] * dt * d_beta[i][m] / den;
        }
        d_chi[g][m] = chi;
        // synthetic fission cross section
        d_sigma_f[g][m] /= d_kcrit;
      } // end groups
    } // end materials
  }
  else
  {
    // Just scale by the eigenvalue.  This must be done
    // to get the right precursor concentration for initial
    // conditions.
    for (size_t m = 0; m < number_materials(); ++m)
    {
      for (size_t g = 0; g < number_groups(); ++g)
      {
        d_sigma_f[g][m] /= d_kcrit;
        //d_chi[g][m] *= (1.0 - beta_total(m));
      }
    }
  }

  // Finalize
  finalize();

}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file TimeDependentMaterial.cc
//---------------------------------------------------------------------------//
