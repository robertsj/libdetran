//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   KineticsMaterial.i.hh
 *  @author robertsj
 *  @date   Oct 2, 2012
 *  @brief  KineticsMaterial inline member definitions
 */
//---------------------------------------------------------------------------//

#ifndef detran_KINETICSMATERIAL_I_HH_
#define detran_KINETICSMATERIAL_I_HH_

namespace detran
{

//---------------------------------------------------------------------------//
inline double KineticsMaterial::velocity(const size_t g) const
{
  Require(g < number_groups());
  return d_velocity[g];
}

//---------------------------------------------------------------------------//
inline double KineticsMaterial::lambda(const size_t i) const
{
  Require(i < d_number_precursor_groups);
  return d_lambda[i];
}

//---------------------------------------------------------------------------//
inline double KineticsMaterial::beta(const size_t m, const size_t i) const
{
  Require(m < number_materials());
  Require(i < d_number_precursor_groups);
  return d_beta[i][m];
}

//---------------------------------------------------------------------------//
inline double KineticsMaterial::beta_total(const size_t m) const
{
  Require(m < number_materials());
  double value = 0;
  for (size_t i = 0; i < d_number_precursor_groups; ++i)
    value += d_beta[i][m];
  return value;
}

//---------------------------------------------------------------------------//
inline double KineticsMaterial::chi_d(const size_t m,
                                      const size_t i,
                                      const size_t g) const
{
  Require(m < number_materials());
  Require(i < d_number_precursor_groups);
  Require(g < number_groups());
  return d_chi_d[g][i][m];
}

//---------------------------------------------------------------------------//
inline KineticsMaterial::size_t
KineticsMaterial::number_precursor_groups() const
{
  return d_number_precursor_groups;
}

} // end namespace detran


#endif /* detran_KINETICSMATERIAL_I_HH_ */
