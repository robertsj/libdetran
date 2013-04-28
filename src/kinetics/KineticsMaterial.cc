//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   KineticsMaterial.cc
 *  @author robertsj
 *  @date   Oct 2, 2012
 *  @brief  KineticsMaterial class definition.
 */
//---------------------------------------------------------------------------//

#include "KineticsMaterial.hh"
#include "utilities/SoftEquivalence.hh"
#include <cstdio>

namespace detran
{

//---------------------------------------------------------------------------//
KineticsMaterial::KineticsMaterial(const size_t num_materials,
                                   const size_t num_energy_groups,
                                   const size_t num_precursor_groups,
                                   std::string  name)
  : Material(num_materials, num_energy_groups, name)
  , d_number_precursor_groups(num_precursor_groups)
  , d_velocity(number_groups(), 0.0)
  , d_lambda(d_number_precursor_groups, 0.0)
  , d_beta(d_number_precursor_groups, vec_dbl(number_materials(), 0.0))
  , d_chi_d(d_number_groups,
            vec2_dbl(d_number_precursor_groups,
                     vec_dbl(number_materials(), 0.0)))
  , d_verified(false)
{
  /* ... */
}

//---------------------------------------------------------------------------//
KineticsMaterial::SP_material
KineticsMaterial::Create(const size_t number_materials,
                         const size_t number_energy_groups,
                         const size_t number_precursor_groups,
                         std::string  name)
{
  SP_material p(new KineticsMaterial(number_materials,
                                     number_energy_groups,
                                     number_precursor_groups,
                                     name));
  // Postconditions
  Ensure(p);
  return p;
}

//---------------------------------------------------------------------------//
void KineticsMaterial::set_velocity(const size_t g, const double v)
{
  Require(g < d_number_groups);
  Require(v > 0.0);
  d_velocity[g] = v;
}

//---------------------------------------------------------------------------//
void KineticsMaterial::set_lambda(const size_t i, const double v)
{
  Require(i < d_number_precursor_groups);
  Require(v > 0.0);
  d_lambda[i] = v;
}

//---------------------------------------------------------------------------//
void KineticsMaterial::set_beta(const size_t m, const size_t i, const double v)
{
  Require(m < number_materials());
  Require(i < d_number_precursor_groups);
  //Require(v > 0.0);
  d_beta[i][m] = v;
}

//---------------------------------------------------------------------------//
void KineticsMaterial::set_beta(const size_t i, const double v)
{
  Require(i < d_number_precursor_groups);
  //Require(v > 0.0);
  for (size_t m = 0; m < number_materials(); ++m)
    d_beta[i][m] = v;
}

//---------------------------------------------------------------------------//
void KineticsMaterial::set_chi_d(const size_t m,
                                 const size_t i,
                                 const size_t g,
                                 const double v)
{
  Require(m < number_materials());
  Require(i < d_number_precursor_groups);
  Require(g < number_groups());
  d_chi_d[g][i][m] = v;
}

//---------------------------------------------------------------------------//
void KineticsMaterial::finalize()
{
  Material::finalize();
  for (size_t i = 0; i < d_number_precursor_groups; ++i)
  {
    Insist(d_lambda[i] > 0.0, "Decay constants must be greater than zero.");
  }
  for (size_t g = 0; g < number_groups(); ++g)
  {
    Insist(d_velocity[g] > 0.0, "Velocities must be greater than zero.");
  }

  // Check that chi's add to 1
  for (size_t m = 0; m < number_materials(); ++m)
  {
    bool has_fission = false;
    double v = 0;
    for (size_t g = 0; g < number_groups(); ++g)
    {
      if (d_sigma_f[g][m] > 0.0) has_fission = true;
      v += d_chi[g][m];
    }
    Assert(v > 0.0 || !has_fission);

    for (size_t i = 0; i < d_number_precursor_groups; ++i)
    {
      v = 0;
      for (size_t g = 0; g < number_groups(); ++g)
      {
        v += d_chi_d[g][i][m];
      }
      Assert(v > 0.0 || !has_fission);
    }
  }

  // Note, beta is not checked, since we don't do anything special for
  // non fissile fuel.
  d_finalized = true;
}

//---------------------------------------------------------------------------//
void KineticsMaterial::display()
{
  // Display base data
  material_display();
  // Display kinetics data
  kinetics_display();
}


//---------------------------------------------------------------------------//
// IMPLEMENTATION
//---------------------------------------------------------------------------//

void KineticsMaterial::kinetics_display()
{
  printf("\n");
  printf("Number of precursor groups: %5i\n", d_number_precursor_groups);
  printf("\n");
  printf("Decay constants [1/s]:\n");
  for (size_t i = 0; i < d_number_precursor_groups; ++i)
  {
    printf("%5i  %16.8f\n", i, d_lambda[i]);
  }
  printf("\n");
  printf("Velocities [cm/s]:\n");
  for (size_t i = 0; i < d_number_groups; ++i)
  {
    printf("%5i  %16.8f\n", i, d_velocity[i]);
  }
  printf("\n");
  printf("Delayed fractions:\n");
  for (size_t m = 0; m < d_number_materials; m++)
  {
    printf("material      %5i\n", m);
    for (size_t i = 0; i < d_number_precursor_groups; i++)
    {
      printf("%5i %16.8f \n", i, d_beta[i][m]);
    }
    printf("TOTAL %16.8f \n", beta_total(m));
    printf("\n");
  } // end material loop
}

} // end detran
