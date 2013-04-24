//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   PyPyTimeDependentMaterial.hh
 *  @author robertsj
 *  @date   Nov 21, 2012
 *  @brief  PyPyTimeDependentMaterial class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_PYTIMEDEPENDENTMATERIAL_HH_
#define detran_PYTIMEDEPENDENTMATERIAL_HH_

#include "TimeDependentMaterial.hh"

namespace detran
{

/**
 *  @class PyTimeDependentMaterial
 *  @brief Base class for time-dependent materials defined in Python
 */
class KINETICS_EXPORT PyTimeDependentMaterial: public TimeDependentMaterial
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef TimeDependentMaterial             Base;
  typedef void (*callback_ptr)(void *);

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param number_materials           Number of materials
   *  @param number_energy_groups       Number of energy groups
   *  @param number_precursor_groups    Number of precursor groups
   *  @param state                      State vector
   */
  PyTimeDependentMaterial(const size_t number_materials,
                          const size_t number_energy_groups,
                          const size_t number_precursor_groups,
                          std::string  name = "PyTimeDependentMaterial");

  /// SP constructor
  static SP_material Create(const size_t number_materials,
                            const size_t number_energy_groups,
                            const size_t number_precursor_groups,
                            std::string  name = "PyTimeDependentMaterial");

  /// Virtual destructor
  virtual ~PyTimeDependentMaterial(){}

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /// Set the user defined material update.
  void set_update_impl(callback_ptr f, void* data);

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  // Pointer to user-defined update function
  callback_ptr d_update_impl;

  // Data for callback (usually unused)
  void* d_update_impl_data;

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL TIME DEPENDENT MATERIALS MUST IMPLEMENT THESE
  //-------------------------------------------------------------------------//

  /// Call the callback function
  void update_impl();

};

} // end namespace detran

#endif /* detran_PYTIMEDEPENDENTMATERIAL_HH_ */
