//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   TimeDependentMaterial.cc
 *  @brief  TimeDependentMaterial
 *  @author Jeremy Roberts
 *  @date   Nov 14, 2012
 */
//---------------------------------------------------------------------------//

#include "TimeDependentMaterial.hh"

namespace detran
{

//---------------------------------------------------------------------------//
TimeDependentMaterial::
TimeDependentMaterial(const size_t number_materials,
                      const size_t number_energy_groups,
                      const size_t number_precursor_groups,
                      std::string  name)
  : Base(number_materials, number_energy_groups, number_precursor_groups, name)
  , d_t(0.0)
  , d_dt(0.0)
  , d_kcrit(1.0)
{
  /* ... */
}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file TimeDependentMaterial.cc
//---------------------------------------------------------------------------//
