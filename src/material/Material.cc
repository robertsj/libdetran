//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  Material.cc
 *  @brief Material class member definitions
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#include "Material.hh"
#include "utilities/Warning.hh"
#include <iostream>
#include <cstdio>
#include <cmath>

namespace detran_material
{

//---------------------------------------------------------------------------//
/*
 * For now, do the simplest implementation, i.e. allocate all things
 * at once.  For some problems, everything (e.g. fission or diffusion)
 * might be unused, but if it's an issue, write another constructor.
 */
Material::Material(const size_t number_materials,
                   const size_t number_groups,
                   std::string  name)
 : d_name(name)
 , d_number_groups(number_groups)
 , d_number_materials(number_materials)
 , d_sigma_t(number_groups, vec_dbl(number_materials, 0.0))
 , d_sigma_a(number_groups, vec_dbl(number_materials, 0.0))
 , d_nu_sigma_f(number_groups, vec_dbl(number_materials, 0.0))
 , d_sigma_f(number_groups, vec_dbl(number_materials, 0.0))
 , d_nu(number_groups, vec_dbl(number_materials, 1.0))
 , d_chi(number_groups, vec_dbl(number_materials, 0.0))
 , d_sigma_s(number_groups,
		         vec2_dbl(number_groups,
		    		          vec_dbl(number_materials, 0.0)))
 , d_diff_coef(number_groups, vec_dbl(number_materials, 0.0))
 , d_scatter_bounds(number_groups, vec_size_t(4, 0))
 , d_finalized(false)
{
  Ensure(d_sigma_t.size() == number_groups);
  Ensure(d_sigma_t[0].size() == number_materials);
  d_downscatter[0] = false;
  d_downscatter[1] = false;
  d_upscatter_cutoff[0] = d_number_groups;
  d_upscatter_cutoff[1] = 0;
}

//---------------------------------------------------------------------------//
Material::SP_material
Material::Create(const size_t number_materials,
                 const size_t number_groups,
                 std::string  name)
{
  SP_material p(new Material(number_materials, number_groups, name));
  return p;
}

//----------------------------------------------------------------------------//
// Setters
//----------------------------------------------------------------------------//

// Single Values

//---------------------------------------------------------------------------//
void Material::set_sigma_t(size_t m, size_t g, double v)
{
  Require(m < d_number_materials);
  Require(g < d_number_groups);
  Require(v >= 0.0);
  d_sigma_t[g][m] = v;
}

//---------------------------------------------------------------------------//
void Material::set_sigma_a(size_t m, size_t g, double v)
{
  Require(m < d_number_materials);
  Require(g < d_number_groups);
  Require(v >= 0.0);
  d_sigma_a[g][m] = v;
}

//---------------------------------------------------------------------------//
void Material::set_nu_sigma_f(size_t m, size_t g, double v)
{
  Require(m < d_number_materials);
  Require(g < d_number_groups);
  Require(v >= 0.0);
  d_nu_sigma_f[g][m] = v;
}

//---------------------------------------------------------------------------//
void Material::set_sigma_f(size_t m, size_t g, double v)
{
  Require(m < d_number_materials);
  Require(g < d_number_groups);
  Require(v >= 0.0);
  d_sigma_f[g][m] = v;
}

//---------------------------------------------------------------------------//
void Material::set_nu(size_t m, size_t g, double v)
{
  Require(m < d_number_materials);
  Require(g < d_number_groups);
  Require(v >= 0.0);
  d_nu[g][m] = v;
}

//---------------------------------------------------------------------------//
void Material::set_chi(size_t m, size_t g, double v)
{
  Require(m < d_number_materials);
  Require(g < d_number_groups);
  Require(v >= 0.0);
  d_chi[g][m] = v;
}

//---------------------------------------------------------------------------//
// Note, not including anisotropic for now...
void Material::set_sigma_s(size_t m, size_t g, size_t gp, double v)
{
  Require(m < d_number_materials);
  Require(g < d_number_groups);
  Require(gp < d_number_groups);
  Require(v >= 0.0);
  d_sigma_s[g][gp][m] = v;
}

//---------------------------------------------------------------------------//
void Material::set_diff_coef(size_t m, size_t g, double v)
{
  Require(m < d_number_materials);
  Require(g < d_number_groups);
  Require(v >= 0.0);
  d_diff_coef[g][m] = v;
}

// Vectorized

//---------------------------------------------------------------------------//
void Material::set_sigma_t(size_t m, vec_dbl &v)
{
  Require(m < d_number_materials);
  Require(v.size() == d_number_groups);
  for (size_t g = 0; g < d_number_groups; g++)
    d_sigma_t[g][m] = v[g];
}

//---------------------------------------------------------------------------//
void Material::set_sigma_a(size_t m, vec_dbl &v)
{
  Require(m < d_number_materials);
  Require(v.size() == d_number_groups);
  for (size_t g = 0; g < d_number_groups; g++)
    d_sigma_a[g][m] = v[g];
}

//---------------------------------------------------------------------------//
void Material::set_nu_sigma_f(size_t m, vec_dbl &v)
{
  Require(m < d_number_materials);
  Require(v.size() == d_number_groups);
  for (size_t g = 0; g < d_number_groups; g++)
    d_nu_sigma_f[g][m] = v[g];
}

//---------------------------------------------------------------------------//
void Material::set_sigma_f(size_t m, vec_dbl &v)
{
  Require(m < d_number_materials);
  Require(v.size() == d_number_groups);
  for (size_t g = 0; g < d_number_groups; g++)
    d_sigma_f[g][m] = v[g];
}

//---------------------------------------------------------------------------//
void Material::set_nu(size_t m, vec_dbl &v)
{
  Require(m < d_number_materials);
  Require(v.size() == d_number_groups);
  for (size_t g = 0; g < d_number_groups; g++)
    d_nu[g][m] = v[g];
}

//---------------------------------------------------------------------------//
void Material::set_chi(size_t m, vec_dbl &v)
{
  Require(m < d_number_materials);
  Require(v.size() == d_number_groups);
  for (size_t g = 0; g < d_number_groups; g++)
    d_chi[g][m] = v[g];
}

//---------------------------------------------------------------------------//
void Material::set_sigma_s(size_t m, size_t g, vec_dbl &v)
{
  Require(m < d_number_materials);
  Require(g < d_number_groups);
  Require(v.size() == d_number_groups);
  for (size_t gp = 0; gp < d_number_groups; gp++)
    d_sigma_s[g][gp][m] = v[gp];
}

//---------------------------------------------------------------------------//
void Material::set_diff_coef(size_t m, vec_dbl &v)
{
  Require(m < d_number_materials);
  Require(v.size() == d_number_groups);
  for (size_t g = 0; g < d_number_groups; g++)
    d_diff_coef[g][m] = v[g];
}

//----------------------------------------------------------------------------//
// Getters
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
Material::size_t Material::lower(size_t g, bool tran) const
{
  Require(d_finalized);
  Require(g < d_number_groups);
  if (tran) return d_scatter_bounds[g][2];
  return d_scatter_bounds[g][0];
}

//----------------------------------------------------------------------------//
Material::size_t Material::upper(size_t g, bool tran) const
{
  Require(d_finalized);
  Require(g < d_number_groups);
  if (tran) return d_scatter_bounds[g][3];
  return d_scatter_bounds[g][1];
}

//----------------------------------------------------------------------------//
bool Material::downscatter(bool tran) const
{
  if (tran) return d_downscatter[1];
  return d_downscatter[0];
}

//----------------------------------------------------------------------------//
Material::size_t Material::upscatter_cutoff(bool tran) const
{
  if (tran) return d_upscatter_cutoff[1];
  return d_upscatter_cutoff[0];
}

//----------------------------------------------------------------------------//
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

//----------------------------------------------------------------------------//
void Material::compute_diff_coef()
{
  double coef = 1.0 / 3.0;
  for (size_t m = 0; m < d_number_materials; m++)
  {
    for (size_t g = 0; g < d_number_groups; g++)
    {
      d_diff_coef[g][m] =  coef / d_sigma_t[g][m];
    }
  }
}

//----------------------------------------------------------------------------//
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
  for (size_t g = 0; g < d_number_groups; g++)
  {
    size_t lower = g;
    size_t upper = g;
    for (size_t m = 0; m < d_number_materials; m++)
    {
      // Downscatter from gp to g
      for (size_t gp = 0; gp < g; gp++)
        if (d_sigma_s[g][gp][m] > 0.0) lower = std::min(gp, lower);
      // Upscatter from gp to g
      for (size_t gp = 0; gp < d_number_groups; gp++)
        if (d_sigma_s[g][gp][m] > 0.0) upper = std::max(gp, upper);
    }
    d_scatter_bounds[g][0] = lower;           // lowest gp into a row g
    d_scatter_bounds[g][1] = upper;           // highest gp into a row g
  }
  // Now for the transpose
  for (size_t gp = 0; gp < d_number_groups; gp++)
  {
    size_t lower = gp;
    size_t upper = gp;
    for (size_t m = 0; m < d_number_materials; m++)
    {
      // Downscatter from g to gp
      for (size_t g = 0; g < gp; g++)
        if (d_sigma_s[g][gp][m] > 0.0) lower = std::min(g, lower);
      // Upscatter from gp to g
      for (size_t g = 0; g < d_number_groups; g++)
        if (d_sigma_s[g][gp][m] > 0.0) upper = std::max(g, upper);
    }
    d_scatter_bounds[gp][2] = upper; // Numerically, the bounds are
    d_scatter_bounds[gp][3] = lower; // in reverse order for adjoint.
  }

  /*
   * Go through the scatter bounds for each g.  If for some g, the
   * upper scatter bound is larger than g, then upscatter exists
   * from that lower energy group.  The first group g for which
   * this occurs is the upscatter cutoff.
   *
   */
  for (size_t g = 0; g < d_number_groups; ++g)
  {
    if (d_scatter_bounds[g][1] > g)
    {
      d_upscatter_cutoff[0] = g;
      break;
    }
  }
  // Now for transpose.  Here, we want the highest energy group from which
  // there is no downscatter.
  for (int g = d_number_groups - 1; g >= 0; --g)
  {
    if (d_scatter_bounds[g][3] < g)
    {
      d_upscatter_cutoff[1] = g;
      break;
    }
  }

  // If our materials have no upscatter, then we set the downscatter-only flag.
  if (d_upscatter_cutoff[0] == d_number_groups)
    d_downscatter[0] = true;
  if (d_upscatter_cutoff[1] == d_number_groups)
    d_downscatter[1] = true;

  // Compute nu*sigma_f
  for (size_t g = 0; g < d_number_groups; g++)
    for (size_t m = 0; m < d_number_materials; m++)
      d_nu_sigma_f[g][m] = d_nu[g][m] * d_sigma_f[g][m];

  d_finalized = true;
}

void Material::display()
{
  material_display();
}

//----------------------------------------------------------------------------//
// IMPLEMENTATION
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
void Material::material_display()
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
  std::cout << "Detran Material: " << d_name << std::endl;
  printf("\n");
  printf("Number of materials: %5i\n", d_number_materials);
  printf("   Number of groups: %5i\n", d_number_groups);
  printf("\n");
  for (size_t m = 0; m < d_number_materials; m++)
  {
    printf("material      %5i\n", m);
    printf("     gp: ");
    for (size_t g = 0; g < d_number_groups; g++)
      printf("%13i ", g);

    printf("\n  total  ");
    for (size_t g = 0; g < d_number_groups; g++)
      printf("%13.6e ", sigma_t(m, g));

    printf("\n  absrb  ");
    for (size_t g = 0; g < d_number_groups; g++)
      printf("%13.6e ", sigma_a(m, g));

    printf("\n  nufis  ");
    for (size_t g = 0; g < d_number_groups; g++)
      printf("%13.6e ", nu_sigma_f(m, g));

    printf("\n  fiss   ");
    for (size_t g = 0; g < d_number_groups; g++)
      printf("%13.6e ", sigma_f(m, g));

    printf("\n  nu     ");
    for (size_t g = 0; g < d_number_groups; g++)
      printf("%13.6e ", nu(m, g));

    printf("\n  chi    ");
    for (size_t g = 0; g < d_number_groups; g++)
      printf("%13.6e ", chi(m, g));

    printf("\n  dc     ");
    for (size_t g = 0; g < d_number_groups; g++)
      printf("%13.6e ", diff_coef(m, g));

    printf("\n");
    for (size_t gp = 0; gp < d_number_groups; gp++)
    {
      printf("%3i<-gp  ", gp);
      for (size_t g = 0; g < d_number_groups; g++)
        printf("%13.6e ", sigma_s(m, gp, g));
      printf("\n");
    }
    printf("\n");

  } // end material loop

  // Other info
  printf("scattering bounds: \n");
  for (size_t g = 0; g < d_number_groups; g++)
    printf("  group %3i  lower = %3i  upper %3i \n", g, lower(g), upper(g));
  printf("\n");
  printf("upscatter cutoff %4i : \n\n", d_upscatter_cutoff[0]);

  printf("adjoint scattering bounds: \n");
  for (size_t g = 0; g < d_number_groups; g++)
    printf("  group %3i  lower = %3i  upper %3i \n", g, lower(g, true), upper(g, true));
  printf("\n");
  printf("downscatter cutoff %4i : \n", d_upscatter_cutoff[1]);
}

} // end namespace detran_material
