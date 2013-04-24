//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   LinearMaterial.cc
 *  @brief  LinearMaterial
 *  @author Jeremy Roberts
 *  @date   Nov 14, 2012
 */
//---------------------------------------------------------------------------//

#include "LinearMaterial.hh"

namespace detran
{

//---------------------------------------------------------------------------//
LinearMaterial::LinearMaterial(const vec_dbl      &times,
                               const vec_material &materials,
                               std::string        name)
  : Base(materials[0]->number_materials(),
         materials[0]->number_groups(),
         materials[0]->number_precursor_groups(),
         name)
  , d_times(times)
  , d_materials(materials)
{
  // Preconditions
  Require(d_times.size() > 0);
  Require(d_times.size() == d_materials.size());
  Require(d_materials[0]);
  for (size_t i = 0; i < d_times.size(); ++i)
  {
    // Strictly monotonic increasing
    if (i > 0)
    {
      Require(d_times[i - 1] < d_times[i]);
    }
    // Material better be there
    Require(d_materials[i]);
    // Materials have same properties
    Require(d_materials[i]->number_materials() ==
            d_number_materials);
    Require(d_materials[i]->number_groups() ==
            d_number_groups);
    Require(d_materials[i]->number_precursor_groups() ==
            d_number_precursor_groups);
  }

  for (size_t g = 0; g < d_number_groups; ++g)
    d_velocity[g] = d_materials[0]->velocity(g);

  for (size_t i = 0; i < d_number_precursor_groups; ++i)
    d_lambda[i] = d_materials[0]->lambda(i);

  d_number_times = d_times.size();
}

//---------------------------------------------------------------------------//
LinearMaterial::SP_material
LinearMaterial::Create(const vec_dbl      &times,
                       const vec_material &materials,
                       std::string        name)
{
  SP_material p(new LinearMaterial(times, materials, name));
  return p;
}


//---------------------------------------------------------------------------//
// IMPLEMENTATION
//---------------------------------------------------------------------------//


//---------------------------------------------------------------------------//
void LinearMaterial::update_impl()
{

  // Determine which materials to interpolate and the interpolation,
  // factors, i.e. sigma = A*sigma_0 + B*sigma_1
  size_t index_A = 0;
  size_t index_B = 0;
  double A = 1.0;
  double B = 0.0;
  if (d_t <= d_times[0])
  {
    /* ... */
  }
  else if (d_t > d_times[d_number_times-1])
  {
    index_A = d_number_times - 1;
    index_B = index_A;
  }
  else
  {
    for (size_t i = 1; i < d_number_times; ++i)
      if (d_t > d_times[i-1] && d_t <= d_times[i])
        index_B = i;
    // f = fA + slope * delta = fA + (fB-fA)/(tB-tA) * delta
    //   = fA*(1 - delta/(tB-tA)) + fB*delta/(tB-tA)
    //   = (1-B)*fA + B*fB
    index_A = index_B - 1;
    B = (d_t - d_times[index_A]) /
        (d_times[index_B] - d_times[index_A]);
    A = 1.0 - B;
  }

  // Fill internal material
  update_material(index_A, index_B, A, B);

}

//---------------------------------------------------------------------------//
void LinearMaterial::update_material(const size_t iA,
                                     const size_t iB,
                                     const double cA,
                                     const double cB)
{
  SP_kineticsmaterial A = d_materials[iA];
  SP_kineticsmaterial B = d_materials[iB];

  for (size_t m = 0; m < A->number_materials(); ++m)
  {
    for (size_t g = 0; g < A->number_groups(); ++g)
    {

      d_sigma_t[g][m]    = cA*A->sigma_t(m, g) + cB*B->sigma_t(m, g);
      d_sigma_a[g][m]    = cA*A->sigma_a(m, g) + cB*B->sigma_a(m, g);
      d_sigma_f[g][m]    = cA*A->sigma_f(m, g) + cB*B->sigma_f(m, g);
      d_nu[g][m]         = cA*A->nu(m, g) + cB*B->nu(m, g);
      d_diff_coef[g][m]  = cA*A->diff_coef(m, g)  + cB*B->diff_coef(m, g);

      for (size_t gp = 0; gp < A->number_groups(); ++gp)
        d_sigma_s[g][gp][m] = cA*A->sigma_s(m, g, gp) + cB*B->sigma_s(m, g, gp);

      for (size_t i = 0; i < d_number_precursor_groups; ++i)
        d_chi_d[g][i][m] = cA*A->chi_d(m, i, g)  + cB*B->chi_d(m, i, g);

      d_chi[g][m] = cA * A->chi(m, g) + cB * B->chi(m, g);

    } // end groups

    for (size_t i = 0; i < d_number_precursor_groups; ++i)
      d_beta[i][m] = cA*A->beta(m, i)  + cB*B->beta(m, i);

  } // end materials

}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file LinearMaterial.cc
//---------------------------------------------------------------------------//
