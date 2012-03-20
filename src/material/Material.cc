/*
 * Material.cc
 *
 *  Created on: Mar 20, 2012
 *      Author: robertsj
 */

#include "Material.hh"

namespace detran
{

/*
 * For now, do the simplest implementation, i.e. allocate all things
 * at once.  For some problems, everything (e.g. fission or diffusion)
 * might be unused, but if it's an issue, write another constructor.
 */
Material::Material(int number_groups,
                   int number_materials,
                   bool downscatter)
 : d_number_groups(number_groups)
 , d_number_materials(number_materials)
 , d_downscatter(downscatter)
 , d_sigma_t(number_materials, vec_dbl(number_groups, 0.0))
 , d_sigma_a(number_materials, vec_dbl(number_groups, 0.0))
 , d_nu_sigma_f(number_materials, vec_dbl(number_groups, 0.0))
 , d_chi(number_materials, vec_dbl(number_groups, 0.0))
 , d_diff_coef(number_materials, vec_dbl(number_groups, 0.0))
 , d_sigma_s(number_materials, vec2_dbl(number_groups, vec_dbl(number_groups, 0.0)))
 , d_scatter_bounds(number_groups, vec_int(2, 0))
{
  // Preconditions
  Require(number_groups    >= 0);
  Require(number_materials >= 0);

  // Postconditions
  Ensure(d_sigma_t.size() == number_materials);
  Ensure(d_sigma_t[0].size() == number_groups);
}

//----------------------------------------------------------------------------//
// Setters
//----------------------------------------------------------------------------//

void Material::set_sigma_t(int m, int g, double v)
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(g >= 0);
  Require(g < d_number_groups);
  Require(v >= 0.0);
  d_sigma_t[m][g] = v;
}

void Material::set_sigma_a(int m, int g, double v)
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(g >= 0);
  Require(g < d_number_groups);
  Require(v >= 0.0);
  d_sigma_a[m][g] = v;
}


void Material::set_nu_sigma_f(int m, int g, double v)
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(g >= 0);
  Require(g < d_number_groups);
  Require(v >= 0.0);
  d_nu_sigma_f[m][g] = v;
}

void Material::set_chi(int m, int g, double v)
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(g >= 0);
  Require(g < d_number_groups);
  Require(v >= 0.0);
  d_chi[m][g] = v;
}

// Note, not including anisotropic for now...
void Material::set_sigma_s(int m, int g, int gp, double v)
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(g >= 0);
  Require(g < d_number_groups);
  Require(gp >= 0);
  Require(gp < d_number_groups);
  Require(v >= 0.0);
  d_sigma_s[m][g][gp] = v;
}

void Material::set_diff_coef(int m, int g, double v)
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(g >= 0);
  Require(g < d_number_groups);
  Require(v >= 0.0);
  d_diff_coef[m][g] = v;
}

void Material::set_sigma_t(int m, vec_dbl &v)
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(v.size() == d_number_groups);
  d_sigma_t[m] = v;
}

void Material::set_nu_sigma_f(int m, vec_dbl &v)
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(v.size() == d_number_groups);
  d_nu_sigma_f[m] = v;
}

void Material::set_chi(int m, vec_dbl &v)
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(v.size() == d_number_groups);
  d_chi[m] = v;
}

void Material::set_sigma_s(int m, int g, vec_dbl &v)
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(g >= 0);
  Require(g < d_number_groups);
  Require(v.size() == d_number_groups);
  d_sigma_s[m][g] = v;
}

void Material::set_diff_coef(int m, vec_dbl &v)
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(v.size() == d_number_groups);
  d_diff_coef[m] = v;
}

//----------------------------------------------------------------------------//
// Getters
//----------------------------------------------------------------------------//

int Material::lower(int g)
{
  Require(d_finalized);
  Require(g >= 0);
  Require(g < d_number_groups);
  return d_scatter_bounds[g][0];
}

int Material::upper(int g)
{
  Require(d_finalized);
  Require(g >= 0);
  Require(g < d_number_groups);
  return d_scatter_bounds[g][1];
}

void Material::finalize()
{

  /*
   * Set the scatter group bounds.  For each group, we compute the
   * lowest index (highest energy) that leads to downscatter.  We also
   * compute the highest index (lowest energy) that upscatters into the
   * group.  Knowing these bounds eliminates a bit of computation in
   * computing the scattering source.
   *
   */
  for (int g = 0; g < d_number_groups; g++)
  {
    int lower = g;
    int upper = g;

    for (int m = 0; m < d_number_materials; m++)
    {
      // Downscatter from gp to g
      for (int gp = 0; gp < g; gp++)
      {
        if (d_sigma_s[m][g][gp] > 0.0) lower = std::min(gp, lower);
      }

      // Upscatter from gp to g
      for (int gp = 0; gp < d_number_groups; gp++)
      {
        if (d_sigma_s[m][g][gp] > 0.0) upper = std::max(gp, upper);
      }

    }
    d_scatter_bounds[g][0] = lower;
    d_scatter_bounds[g][1] = upper;
  }

  /*
   * Go through the scatter bounds for each g.  If for some g, the
   * upper scatter bound is larger than g, then upscatter exists
   * from that lower energy group.  The first group g for which
   * this occurs is the upscatter cutoff.
   *
   */
  d_upscatter_cutoff = d_number_groups;
  for (int g = 0; g < d_number_groups; g++)
  {
    if (d_scatter_bounds[g][1] > g)
    {
      d_upscatter_cutoff = g;
      break;
    }
  }

  // If our materials have no upscatter, then we set the
  // downscatter-only flag.
  if (d_upscatter_cutoff == d_number_groups)
  {
    if (d_downscatter == false)
    {
      detran_utils::warning(detran_utils::USER_INPUT,
        "Upscatter is being turned off since no upscatter exists in the data.");
    }
    d_downscatter = true;
  }

  d_finalized = true;
}


} // end namespace detran
