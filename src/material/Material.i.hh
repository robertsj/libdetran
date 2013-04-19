//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Material.i.hh
 *  @author Jeremy Roberts
 *  @brief  Material class inline member definitions.
 */
//---------------------------------------------------------------------------//

#ifndef detran_material_MATERIAL_I_HH_
#define detran_material_MATERIAL_I_HH_

#include "utilities/DBC.hh"

namespace detran_material
{

//---------------------------------------------------------------------------//
inline double Material::sigma_t(size_t m, size_t g) const
{
  Require(m < d_number_materials);
  Require(g < d_number_groups);
  return d_sigma_t[g][m];
}

//---------------------------------------------------------------------------//
inline double Material::sigma_a(size_t m, size_t g) const
{
  Require(m < d_number_materials);
  Require(g < d_number_groups);
  return d_sigma_a[g][m];
}

//---------------------------------------------------------------------------//
inline double Material::nu_sigma_f(size_t m, size_t g) const
{
  Require(m < d_number_materials);
  Require(g < d_number_groups);
  return d_nu_sigma_f[g][m];
}

//---------------------------------------------------------------------------//
inline double Material::sigma_f(size_t m, size_t g) const
{
  Require(m < d_number_materials);
  Require(g < d_number_groups);
  return d_sigma_f[g][m];
}

//---------------------------------------------------------------------------//
inline double Material::nu(size_t m, size_t g) const
{
  Require(m < d_number_materials);
  Require(g < d_number_groups);
  return d_nu[g][m];
}

//---------------------------------------------------------------------------//
inline double Material::chi(size_t m, size_t g) const
{
  Require(m < d_number_materials);
  Require(g < d_number_groups);
  return d_chi[g][m];
}

//---------------------------------------------------------------------------//
inline double Material::sigma_s(size_t m, size_t g, size_t gp) const
{
  Require(m < d_number_materials);
  Require(g < d_number_groups);
  Require(gp < d_number_groups);
  return d_sigma_s[g][gp][m];
}

//---------------------------------------------------------------------------//
inline double Material::diff_coef(size_t m, size_t g) const
{
  Require(m < d_number_materials);
  Require(g < d_number_groups);
  return d_diff_coef[g][m];
}

// Vectorized.  Since we store with material on the
// outer index, we actual have to copy the materials
// to a new vector.  These are largely implemented to
// help in I/O of material data.

//---------------------------------------------------------------------------//
inline Material::vec_dbl Material::sigma_t(size_t m) const
{
  Require(m < d_number_materials);
  vec_dbl v(d_number_groups, 0.0);
  for (size_t g = 0; g < d_number_groups; g++)
    v[g] = d_sigma_t[g][m];
  return v;
}

//---------------------------------------------------------------------------//
inline Material::vec_dbl Material::sigma_a(size_t m) const
{
  Require(m < d_number_materials);
  vec_dbl v(d_number_groups, 0.0);
  for (size_t g = 0; g < d_number_groups; g++)
    v[g] = d_sigma_a[g][m];
  return v;
}

//---------------------------------------------------------------------------//
inline Material::vec_dbl Material::nu_sigma_f(size_t m) const
{
  Require(m < d_number_materials);
  vec_dbl v(d_number_groups, 0.0);
  for (size_t g = 0; g < d_number_groups; g++)
    v[g] = d_nu_sigma_f[g][m];
  return v;
}

//---------------------------------------------------------------------------//
inline Material::vec_dbl Material::sigma_f(size_t m) const
{
  Require(m < d_number_materials);
  vec_dbl v(d_number_groups, 0.0);
  for (size_t g = 0; g < d_number_groups; g++)
    v[g] = d_sigma_f[g][m];
  return v;
}

//---------------------------------------------------------------------------//
inline Material::vec_dbl Material::nu(size_t m) const
{
  Require(m < d_number_materials);
  vec_dbl v(d_number_groups, 0.0);
  for (size_t g = 0; g < d_number_groups; g++)
    v[g] = d_nu[g][m];
  return v;
}

//---------------------------------------------------------------------------//
inline Material::vec_dbl Material::chi(size_t m) const
{
  Require(m < d_number_materials);
  vec_dbl v(d_number_groups, 0.0);
  for (size_t g = 0; g < d_number_groups; g++)
    v[g] = d_chi[g][m];
  return v;
}

//---------------------------------------------------------------------------//
inline Material::vec2_dbl Material::sigma_s(size_t m) const
{
  Require(m < d_number_materials);
  vec2_dbl v(d_number_groups, vec_dbl(d_number_groups, 0.0));
  for (size_t g = 0; g < d_number_groups; g++)
    for (size_t gp = 0; gp < d_number_groups; gp++)
      v[g][gp] = d_sigma_s[g][gp][m];
  return v;
}

//---------------------------------------------------------------------------//
inline Material::vec_dbl Material::diff_coef(size_t m) const
{
  Require(m < d_number_materials);
  vec_dbl v(d_number_groups, 0.0);
  for (size_t g = 0; g < d_number_groups; g++)
    v[g] = d_diff_coef[g][m];
  return v;
}

} // end namespace detran

#endif /* detran_material_MATERIAL_I_HH_ */
