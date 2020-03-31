//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   TimeIndependentMaterial.cc
 *  @brief  TimeInependentMaterial
 *  @author Rabab Elzohery
 *  @date   March 28, 2020
 */
//---------------------------------------------------------------------------//

#include "TimeIndependentMaterial.hh"
#include "KineticsMaterial.hh"

namespace detran
{

//---------------------------------------------------------------------------//
TimeIndependentMaterial::
TimeIndependentMaterial(KineticsMaterial::SP_material km)
    :Base(km->number_materials(), km->number_groups(), km->number_precursor_groups()),
	  d_km(km)
{
	 update_impl();
}

TimeIndependentMaterial::SP_material
TimeIndependentMaterial::Create(KineticsMaterial::SP_material KM)
{
  SP_material p(new TimeIndependentMaterial(KM));

  return p;
}

void TimeIndependentMaterial::update_impl()
{
  size_t m = d_number_materials - 1;
  for (int i=0; i < d_number_precursor_groups; i++)
    {
      d_lambda[i] = d_km->lambda(i);
      for (int m = 0; m < d_number_materials; m++)
        {
          d_beta[i][m] = d_km->beta(m, i);
	  for (int g=0; g < d_number_groups; g++)
          {
	    d_sigma_t[g][m] = d_km->sigma_t(m, g);
	    d_sigma_a[g][m] = d_km->sigma_a(m, g);
	    d_nu_sigma_f[g][m] = d_km->nu_sigma_f(m, g);
	    d_sigma_f[g][m] = d_km->sigma_f(m, g);
	    d_nu[g][m] = d_km->nu(m, g);
	    d_chi[g][m] = d_km->chi(m, g);
	    d_diff_coef[g][m] = d_km->diff_coef(m, g);
	    d_chi_d[g][i][m] = d_km->chi_d(m, i, g);
	    d_velocity[g] = d_km->velocity(g);
	    for (int gp = 0; gp < d_number_groups; gp++)
	    {
	      d_sigma_s[g][gp][m] = d_km->sigma_s(m, g, gp);
	    }
        }
     }
  }
}
}; // end namespace detran

//---------------------------------------------------------------------------//
//              end of file TimeIndependentMaterial.cc
//---------------------------------------------------------------------------//
