/*
 * Material.i.hh
 *
 *  Created on: Mar 20, 2012
 *      Author: robertsj
 */

#ifndef MATERIAL_I_HH_
#define MATERIAL_I_HH_

namespace detran
{

double Material::sigma_t(int m, int g) const
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(g >= 0);
  Require(g < d_number_groups);
  return d_sigma_t[g][m];
}

double Material::sigma_a(int m, int g) const
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(g >= 0);
  Require(g < d_number_groups);
  return d_sigma_a[g][m];
}

double Material::nu_sigma_f(int m, int g) const
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(g >= 0);
  Require(g < d_number_groups);
  return d_nu_sigma_f[g][m];
}

double Material::sigma_f(int m, int g) const
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(g >= 0);
  Require(g < d_number_groups);
  return d_sigma_f[g][m];
}

double Material::nu(int m, int g) const
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(g >= 0);
  Require(g < d_number_groups);
  return d_nu[g][m];
}

double Material::chi(int m, int g) const
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(g >= 0);
  Require(g < d_number_groups);
  return d_chi[g][m];
}

double Material::sigma_s(int m, int g, int gp) const
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(g >= 0);
  Require(g < d_number_groups);
  Require(gp >= 0);
  Require(gp < d_number_groups);
  return d_sigma_s[g][gp][m];
}

double Material::diff_coef(int m, int g) const
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(g >= 0);
  Require(g < d_number_groups);
  return d_diff_coef[g][m];
}

// Vectorized

vec_dbl Material::sigma_t(int m) const
{
  Require(m >= 0);
  Require(m < d_number_materials);
  vec_dbl v(d_number_groups, 0.0);
  for (int g = 0; g < d_number_groups; g++)
    v[g] = d_sigma_t[g][m];
  return v;
}

vec_dbl Material::sigma_a(int m) const
{
  Require(m >= 0);
  Require(m < d_number_materials);
  vec_dbl v(d_number_groups, 0.0);
  for (int g = 0; g < d_number_groups; g++)
    v[g] = d_sigma_a[g][m];
  return v;
}

vec_dbl Material::nu_sigma_f(int m) const
{
  Require(m >= 0);
  Require(m < d_number_materials);
  vec_dbl v(d_number_groups, 0.0);
  for (int g = 0; g < d_number_groups; g++)
    v[g] = d_nu_sigma_f[g][m];
  return v;
}

vec_dbl Material::sigma_f(int m) const
{
  Require(m >= 0);
  Require(m < d_number_materials);
  vec_dbl v(d_number_groups, 0.0);
  for (int g = 0; g < d_number_groups; g++)
    v[g] = d_sigma_f[g][m];
  return v;
}

vec_dbl Material::nu(int m) const
{
  Require(m >= 0);
  Require(m < d_number_materials);
  vec_dbl v(d_number_groups, 0.0);
  for (int g = 0; g < d_number_groups; g++)
    v[g] = d_nu[g][m];
  return v;
}

vec_dbl Material::chi(int m) const
{
  Require(m >= 0);
  Require(m < d_number_materials);
  vec_dbl v(d_number_groups, 0.0);
  for (int g = 0; g < d_number_groups; g++)
    v[g] = d_chi[g][m];
  return v;
}

vec2_dbl Material::sigma_s(int m) const
{
  Require(m >= 0);
  Require(m < d_number_materials);
  vec2_dbl v(d_number_groups, vec_dbl(d_number_groups, 0.0));
  for (int g = 0; g < d_number_groups; g++)
    for (int gp = 0; gp < d_number_groups; gp++)
      v[g][gp] = d_sigma_s[g][gp][m];
  return v;
}

vec_dbl Material::diff_coef(int m) const
{
  Require(m >= 0);
  Require(m < d_number_materials);
  vec_dbl v(d_number_groups, 0.0);
  for (int g = 0; g < d_number_groups; g++)
    v[g] = d_diff_coef[g][m];
  return v;
}

} // end namespace detran

#endif /* MATERIAL_I_HH_ */
