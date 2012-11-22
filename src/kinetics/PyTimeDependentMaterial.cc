//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   PyTimeDependentMaterial.cc
 *  @author robertsj
 *  @date   Nov 21, 2012
 *  @brief  PyTimeDependentMaterial class definition.
 */
//---------------------------------------------------------------------------//

#include "PyTimeDependentMaterial.hh"

namespace detran
{

//---------------------------------------------------------------------------//
PyTimeDependentMaterial::
PyTimeDependentMaterial(const size_t number_materials,
                        const size_t number_energy_groups,
                        const size_t number_precursor_groups,
                        std::string  name)
  : Base(number_materials, number_energy_groups, number_precursor_groups, name)
  , d_update_impl(NULL)
  , d_update_impl_data(NULL)
{
  /* ... */
}

//---------------------------------------------------------------------------//
PyTimeDependentMaterial::SP_material
PyTimeDependentMaterial::Create(const size_t number_materials,
                                const size_t number_energy_groups,
                                const size_t number_precursor_groups,
                                std::string  name)
{
  SP_material p(new PyTimeDependentMaterial(number_materials,
                                            number_energy_groups,
                                            number_precursor_groups,
                                            name));
  return p;
}


//---------------------------------------------------------------------------//
void PyTimeDependentMaterial::set_update_impl(callback_ptr f, void* data)
{
  std::cout << "C++ setting implementation" << std::endl;
  d_update_impl = f;
  d_update_impl_data = data;
}

//---------------------------------------------------------------------------//
void PyTimeDependentMaterial::update_impl()
{
  Require(d_update_impl);
  d_update_impl(d_update_impl_data);
}

} // end namespace detran


