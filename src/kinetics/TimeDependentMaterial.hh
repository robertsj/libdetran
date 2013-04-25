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
class KINETICS_EXPORT TimeDependentMaterial: public KineticsMaterial
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

  /// Set the state vector.
  void set_state(SP_state s)
  {
    d_state = s;
  }

  /// Get the state vector.
  SP_state state()
  {
    return d_state;
  }

  double time() const
  {
    return d_t;
  }

  double dt() const
  {
    return d_dt;
  }

  size_t order() const
  {
    return d_order;
  }

  double kcrit() const
  {
    return d_kcrit;
  }

  /**
   *  @brief Update the materials
   *  @param t          New time (in seconds)
   *  @param dt         Time step (in seconds)
   *  @param order      BDF order
   *  @paral synthetic  Flag for creating a synthetic material
   */
  void update(const double t,
              const double dt,
              const size_t order,
              const bool flag = true);

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
  /// BDF order
  size_t d_order;
  /// Eigenvalue (for scaling fission)
  double d_kcrit;

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL TIME DEPENDENT MATERIALS MUST IMPLEMENT THESE
  //-------------------------------------------------------------------------//

  /*
   *  @brief User-defined material update
   *
   *  This routine is called from within update.  This fills the
   *  internal cross section vectors with their actual values.  The
   *  update function then adds the synthetic components.
   *
   */
  virtual void update_impl() = 0;

};

KINETICS_TEMPLATE_EXPORT(detran_utilities::SP<TimeDependentMaterial>)

} // end namespace detran

#endif /* detran_TIMEDEPENDENTMATERIAL_HH_ */
