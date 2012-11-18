//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   TimeDependentMaterial.hh
 *  @author robertsj
 *  @date   Oct 2, 2012
 *  @brief  TimeDependentMaterial class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_TIMEDEPENDENTMATERIAL_HH_
#define detran_TIMEDEPENDENTMATERIAL_HH_

#include "KineticsMaterial.hh"
#include "transport/State.hh"

namespace detran
{

/**
 *  @class TimeDependentMaterial
 *  @brief Base class for time-dependent materials
 *
 *  This implementation allows for arbitrary time-dependent materials
 *  to be defined, possibly as a function of state.
 *
 *  The client must implement the private implementation of update,
 *  which defines the normal cross section values.  This private
 *  function is called in update, which then adjusts the cross sections
 *  for the time step (i.e. it produces synthetic cross sections that
 *  are dependent on the time increment)
 *
 */
class TimeDependentMaterial: public KineticsMaterial
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef KineticsMaterial                              Base;
  typedef Base::SP_material                             SP_base;
  typedef detran_utilities::SP<TimeDependentMaterial>   SP_material;
  typedef detran::State::SP_state                       SP_state;

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
  TimeDependentMaterial(const size_t number_materials,
                        const size_t number_energy_groups,
                        const size_t number_precursor_groups,
                        SP_state     state,
                        std::string  name = "TimeDependentMaterial");

  /// Virtual destructor
  virtual ~TimeDependentMaterial(){}

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  void set_eigenvalue(const double keff)
  {
    Require(keff > 0.0);
    d_kcrit = keff;
  }

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL TIME DEPENDENT MATERIALS MUST IMPLEMENT THESE
  //-------------------------------------------------------------------------//

  /**
   *  @brief Update the materials
   *  @param t      New time (in seconds)
   *  @param dt     Time step (in seconds)
   */
  virtual void update(const double t, const double dt) = 0;

protected:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// State
  SP_state d_state;
  /// Time
  double d_t;
  /// Time step
  double d_dt;
  /// Eigenvalue (for scaling fission)
  double d_kcrit;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  /* ... */

};

} // end namespace detran

#endif /* detran_TIMEDEPENDENTMATERIAL_HH_ */
