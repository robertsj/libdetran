//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Material.cc
 * \author Jeremy Roberts
 * \brief  Material class member definitions.
 */
//---------------------------------------------------------------------------//

#include "Material.hh"
#include "utilities/Warning.hh"
#include <iostream>
#include <cstdio>
#include <cmath>

namespace detran_material
{

/*
 * For now, do the simplest implementation, i.e. allocate all things
 * at once.  For some problems, everything (e.g. fission or diffusion)
 * might be unused, but if it's an issue, write another constructor.
 */
Material::Material(const size_t number_materials,
                   const size_t number_groups,
                   const bool downscatter)
 : d_number_materials(number_materials)
 , d_number_groups(number_groups)
 , d_downscatter(downscatter)
 , d_sigma_t(number_groups, vec_dbl(number_materials, 0.0))
 , d_sigma_a(number_groups, vec_dbl(number_materials, 0.0))
 , d_nu_sigma_f(number_groups, vec_dbl(number_materials, 0.0))
 , d_sigma_f(number_groups, vec_dbl(number_materials, 0.0))
 , d_nu(number_groups, vec_dbl(number_materials, 1.0))
 , d_chi(number_groups, vec_dbl(number_materials, 0.0))
 , d_diff_coef(number_groups, vec_dbl(number_materials, 0.0))
 , d_sigma_s(number_groups, vec2_dbl(number_groups, vec_dbl(number_materials, 0.0)))
 , d_scatter_bounds(number_groups, vec_int(2, 0))
{
  // Preconditions
  Require(number_groups    >= 0);
  Require(number_materials >= 0);

  // Postconditions
  Ensure(d_sigma_t.size() == number_groups);
  Ensure(d_sigma_t[0].size() == number_materials);
}

//----------------------------------------------------------------------------//
// Setters
//----------------------------------------------------------------------------//

void Material::set_sigma_t(size_t m, size_t g, double v)
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(g >= 0);
  Require(g < d_number_groups);
  Require(v >= 0.0);
  d_sigma_t[g][m] = v;
}

void Material::set_sigma_a(size_t m, size_t g, double v)
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(g >= 0);
  Require(g < d_number_groups);
  Require(v >= 0.0);
  d_sigma_a[g][m] = v;
}


void Material::set_nu_sigma_f(size_t m, size_t g, double v)
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(g >= 0);
  Require(g < d_number_groups);
  Require(v >= 0.0);
  d_nu_sigma_f[g][m] = v;
}

void Material::set_sigma_f(size_t m, size_t g, double v)
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(g >= 0);
  Require(g < d_number_groups);
  Require(v >= 0.0);
  d_sigma_f[g][m] = v;
}

void Material::set_nu(size_t m, size_t g, double v)
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(g >= 0);
  Require(g < d_number_groups);
  Require(v >= 0.0);
  d_nu[g][m] = v;
}

void Material::set_chi(size_t m, size_t g, double v)
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(g >= 0);
  Require(g < d_number_groups);
  Require(v >= 0.0);
  d_chi[g][m] = v;
}

// Note, not including anisotropic for now...
void Material::set_sigma_s(size_t m, size_t g, size_t gp, double v)
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(g >= 0);
  Require(g < d_number_groups);
  Require(gp >= 0);
  Require(gp < d_number_groups);
  Require(v >= 0.0);
  d_sigma_s[g][gp][m] = v;
}

void Material::set_diff_coef(size_t m, size_t g, double v)
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(g >= 0);
  Require(g < d_number_groups);
  Require(v >= 0.0);
  d_diff_coef[g][m] = v;
}

void Material::set_sigma_t(size_t m, vec_dbl &v)
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(v.size() == d_number_groups);
  for (size_t g = 0; g < d_number_groups; g++)
    d_sigma_t[g][m] = v[g];
}

void Material::set_sigma_a(size_t m, vec_dbl &v)
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(v.size() == d_number_groups);
  for (size_t g = 0; g < d_number_groups; g++)
    d_sigma_a[g][m] = v[g];
}


void Material::set_nu_sigma_f(size_t m, vec_dbl &v)
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(v.size() == d_number_groups);
  for (size_t g = 0; g < d_number_groups; g++)
    d_nu_sigma_f[g][m] = v[g];
}

void Material::set_sigma_f(size_t m, vec_dbl &v)
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(v.size() == d_number_groups);
  for (size_t g = 0; g < d_number_groups; g++)
    d_sigma_f[g][m] = v[g];
}

void Material::set_nu(size_t m, vec_dbl &v)
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(v.size() == d_number_groups);
  for (size_t g = 0; g < d_number_groups; g++)
    d_nu[g][m] = v[g];
}

void Material::set_chi(size_t m, vec_dbl &v)
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(v.size() == d_number_groups);
  for (size_t g = 0; g < d_number_groups; g++)
    d_chi[g][m] = v[g];
}

void Material::set_sigma_s(size_t m, size_t g, vec_dbl &v)
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(g >= 0);
  Require(g < d_number_groups);
  Require(v.size() == d_number_groups);
  d_sigma_s[m][g] = v;
  for (size_t gp = 0; gp < d_number_groups; gp++)
    d_sigma_s[g][gp][m] = v[gp];
}

void Material::set_diff_coef(size_t m, vec_dbl &v)
{
  Require(m >= 0);
  Require(m < d_number_materials);
  Require(v.size() == d_number_groups);
  for (size_t g = 0; g < d_number_groups; g++)
    d_diff_coef[g][m] = v[g];
}

//----------------------------------------------------------------------------//
// Getters
//----------------------------------------------------------------------------//

Material::size_t Material::lower(size_t g) const
{
  Require(d_finalized);
  Require(g < d_number_groups);
  return d_scatter_bounds[g][0];
}

Material::size_t Material::upper(size_t g) const
{
  Require(d_finalized);
  Require(g < d_number_groups);
  return d_scatter_bounds[g][1];
}

void Material::compute_sigma_a()
{
  for (size_t m = 0; m < d_number_materials; m++)
  {
    for (size_t g = 0; g < d_number_groups; g++)
    {
      double sa = d_sigma_t[g][m];
      for (size_t gp = 0; gp < d_number_groups; gp++)
        sa -= d_sigma_s[gp][g][m];
      d_sigma_a[g][m] = sa;
    }
  }
}

void Material::compute_diff_coef()
{
  double coef = 1.0 / 3.0;
  for (size_t m = 0; m < d_number_materials; m++)
  {
    for (size_t g = 0; g < d_number_groups; g++)
    {
      d_diff_coef[g][m] =  coef * d_sigma_t[g][m];
    }
  }
}


void Material::finalize()
{

  /*
   * Set the scatter group bounds.  For each group, we compute the
   * lowest index (highest energy) that leads to downscatter.  We also
   * compute the highest index (lowest energy) that upscatters size_to the
   * group.  Knowing these bounds eliminates a bit of computation in
   * computing the scattering source.
   *
   */
  for (size_t g = 0; g < d_number_groups; g++)
  {
    size_t lower = g;
    size_t upper = g;

    for (size_t m = 0; m < d_number_materials; m++)
    {
      // Downscatter from gp to g
      for (size_t gp = 0; gp < g; gp++)
      {
        if (d_sigma_s[g][gp][m] > 0.0) lower = std::min(gp, lower);
      }

      // Upscatter from gp to g
      for (size_t gp = 0; gp < d_number_groups; gp++)
      {
        if (d_sigma_s[g][gp][m] > 0.0) upper = std::max(gp, upper);
      }

      // Compute nu*sigma_f
      d_nu_sigma_f[g][m] = d_nu[g][m] * d_sigma_f[g][m];

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
  for (size_t g = 0; g < d_number_groups; g++)
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
      detran_utilities::warning(detran_utilities::USER_INPUT,
        "Upscatter is being turned off since no upscatter exists in the data.");
    }
    d_downscatter = true;
  }

  d_finalized = true;
}

void Material::display()
{

  using std::cout;
  using std::endl;
  using std::printf;

  /*
   *   "Material 1 Description"
   *
   *    0               1               2               3
   *    sigma_t1        sigma_t2        ...             ...
   *    nu_sigma_f1     nu_sigmaf2      ...             ...
   *    chi1            chi2            ...             ...
   *    sigma_s1<-1     sigma_s1<-2     ...
   *    sigma_s2<-1     ...
   *
   *    "Material 2 ..."
   */
  printf("\n");
  printf("Detran Material");
  printf("\n");
  printf("Number of materials: %5i\n", d_number_materials);
  printf("   Number of groups: %5i\n", d_number_groups);
  printf("\n");
  for (size_t m = 0; m < d_number_materials; m++)
  {
    printf("material      %5i\n", m);
    printf("     gp: ");
    for (size_t g = 0; g < d_number_groups; g++)
    {
      printf("%13i ", g);
    }
    printf("\n  total  ");
    // Total
    for (size_t g = 0; g < d_number_groups; g++)
    {
      printf("%13.10f ", sigma_t(m, g));
    }
    printf("\n  absorption  ");
    // Total
    for (size_t g = 0; g < d_number_groups; g++)
    {
      printf("%13.10f ", sigma_a(m, g));
    }
    printf("\n  nufis  ");
    // nu*Fission
    for (size_t g = 0; g < d_number_groups; g++)
    {
      printf("%13.10f ", nu_sigma_f(m, g));
    }
    printf("\n  fiss   ");
    // Fission
    for (size_t g = 0; g < d_number_groups; g++)
    {
      printf("%13.10f ", sigma_f(m, g));
    }
    printf("\n  nu     ");
    // nu
    for (size_t g = 0; g < d_number_groups; g++)
    {
      printf("%13.10f ", nu(m, g));
    }
    printf("\n  chi    ");
    // Chi
    for (size_t g = 0; g < d_number_groups; g++)
    {
      printf("%13.10f ", chi(m, g));
    }
    printf("\n");
    // Scatter
    for (size_t gp = 0; gp < d_number_groups; gp++)
    {
      printf("%3i<-gp  ", gp);
      for (size_t g = 0; g < d_number_groups; g++)
      {
        printf("%13.10f ", sigma_s(m, gp, g));
      }
      printf("\n");
    }
    printf("\n");
  } // end material loop

  // Other info
  printf("scattering bounds: \n");
  for (size_t g = 0; g < d_number_groups; g++)
  {
    printf("  group %3i  lower = %3i  upper %3i \n", g, lower(g), upper(g));
  }
  printf("\n");
  printf("upscatter cutoff %4i : \n", d_upscatter_cutoff);
}



} // end namespace detran
