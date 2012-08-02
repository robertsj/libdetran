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
  return d_sigma_t[m][g];
}

double Material::sigma_a(int m, int g) const
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(g >= 0);
  Require(g < d_number_groups);
  return d_sigma_a[m][g];
}

double Material::nu_sigma_f(int m, int g) const
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(g >= 0);
  Require(g < d_number_groups);
  return d_nu_sigma_f[m][g];
}

double Material::sigma_f(int m, int g) const
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(g >= 0);
  Require(g < d_number_groups);
  return d_sigma_f[m][g];
}

double Material::nu(int m, int g) const
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(g >= 0);
  Require(g < d_number_groups);
  return d_nu[m][g];
}

double Material::chi(int m, int g) const
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(g >= 0);
  Require(g < d_number_groups);
  return d_chi[m][g];
}

double Material::sigma_s(int m, int g, int gp) const
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(g >= 0);
  Require(g < d_number_groups);
  Require(gp >= 0);
  Require(gp < d_number_groups);
  return d_sigma_s[m][g][gp];
}

double Material::diff_coef(int m, int g) const
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(g >= 0);
  Require(g < d_number_groups);
  return d_diff_coef[m][g];
}

// Vectorized

vec_dbl Material::sigma_t(int m) const
{
  Require(m >= 0);
  Require(m < d_number_materials);
  return d_sigma_t[m];
}

vec_dbl Material::sigma_a(int m) const
{
  Require(m >= 0);
  Require(m < d_number_materials);
  return d_sigma_a[m];
}

vec_dbl Material::nu_sigma_f(int m) const
{
  Require(m >= 0);
  Require(m < d_number_materials);
  return d_nu_sigma_f[m];
}

vec_dbl Material::sigma_f(int m) const
{
  Require(m >= 0);
  Require(m < d_number_materials);
  return d_sigma_f[m];
}

vec_dbl Material::nu(int m) const
{
  Require(m >= 0);
  Require(m < d_number_materials);
  return d_nu[m];
}

vec_dbl Material::chi(int m) const
{
  Require(m >= 0);
  Require(m < d_number_materials);
  return d_chi[m];
}

vec2_dbl Material::sigma_s(int m) const
{
  Require(m >= 0);
  Require(m < d_number_materials);
  return d_sigma_s[m];
}

vec_dbl Material::diff_coef(int m) const
{
  Require(m >= 0);
  Require(m < d_number_materials);
  return d_diff_coef[m];
}

} // end namespace detran

#endif /* MATERIAL_I_HH_ */
