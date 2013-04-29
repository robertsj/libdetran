//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   TimeStepper.hh
 *  @author robertsj
 *  @date   Oct 2, 2012
 *  @brief  TimeStepper class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_TIMESTEPPER_HH_
#define detran_TIMESTEPPER_HH_

#include "FixedSourceManager.hh"
#include "angle/Quadrature.hh"
#include "ioutils/SiloOutput.hh"
#include "kinetics/BDFCoefficients.hh"
#include "kinetics/MultiPhysics.hh"
#include "kinetics/TimeDependentMaterial.hh"
#include "kinetics/TimeDependentExternalSource.hh"
#include "kinetics/SyntheticSource.hh"
#include "utilities/DBC.hh"
#include "utilities/Definitions.hh"
#include "utilities/SP.hh"
#include "transport/State.hh"
#include <cstdio>

namespace detran
{

/**
 *  @class TimeStepper
 *  @brief Solve a time dependent problem by stepping forward in time.
 *
 *  Available time stepping options are the BDF methods of
 *  order 1 through 6 and the implicit midpoint rule.
 *
 */

template <class D>
class TimeStepper
{

public:

  //-------------------------------------------------------------------------//
  // ENUMERATION
  //-------------------------------------------------------------------------//

  enum time_schemes
  {
    IMP, BDF1, BDF2, BDF3, BDF4, BDF5, BDF6, END_TIME_SCHEMES
  };

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<TimeStepper<D> >           SP_timestepper;
  typedef detran_utilities::InputDB::SP_input             SP_input;
  typedef TimeDependentMaterial::SP_material              SP_material;
  typedef detran_geometry::Mesh::SP_mesh                  SP_mesh;
  typedef detran_angle::Quadrature::SP_quadrature         SP_quadrature;
  typedef State::SP_state                                 SP_state;
  typedef detran_utilities::size_t                        size_t;
  typedef FixedSourceManager<D>                           Fixed_T;
  typedef typename Fixed_T::SP_manager                    SP_mg_solver;
  typedef TimeDependentExternalSource::SP_tdsource        SP_tdsource;
  typedef std::vector<SP_tdsource>                        vec_source;
  typedef SyntheticSource::SP_externalsource              SP_syntheticsource;
  typedef SyntheticSource::vec_states                     vec_states;
  typedef Precursors::SP_precursors                       SP_precursors;
  typedef SyntheticSource::vec_precursors                 vec_precursors;
  typedef MultiPhysics::SP_multiphysics                   SP_multiphysics;
  typedef std::vector<SP_multiphysics>                    vec_multiphysics;
  typedef FissionSource::SP_fissionsource                 SP_fissionsource;
  typedef detran_utilities::vec_int                       vec_int;
  typedef detran_ioutils::SiloOutput::SP_silooutput       SP_silooutput;
  /// Pointer to callback function for monitoring
  typedef void (*monitor_pointer)
               (void*, TimeStepper<D>*, int, double, double, int, bool);
  /// Pointer to callback function for updating physics
  typedef void (*multiphysics_pointer)
               (void*, TimeStepper<D>*, double, double);
  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param input
   *  @param material
   */
  TimeStepper(SP_input      input,
              SP_material   material,
              SP_mesh       mesh,
              bool          multiply);

  /// Virtual destructor
  virtual ~TimeStepper(){}

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /// Add a time-dependent external source
  void add_source(SP_tdsource source);

  /// Solve
  void solve(SP_state initial_state);

  /// Getters
  SP_state state() {return d_state;}
  SP_mesh mesh() {return d_mesh;}
  SP_material material() {return d_material;}
  SP_quadrature quadrature() {return d_quadrature;}
  size_t monitor_level() const {return d_monitor_level;}
  SP_precursors precursor() {return d_precursor;}
  SP_multiphysics multiphysics() {return d_multiphysics;}
  SP_fissionsource fissionsource() {return d_fissionsource;}
  double residual_norm() {return d_residual_norm;}

  /// Set a user-defined monitor function.
  void set_monitor(monitor_pointer monitor, void* monitor_data = NULL)
  {
    Require(monitor);
    d_monitor = monitor;
    d_monitor_data = monitor_data;
  }

  /**
   *  @brief Set the multiphysics
   *
   *  The user must initialize the multiphysics state vector.  The
   *  user also specifies the multiphysics update function as well
   *  as any data required.
   */
  void set_multiphysics(SP_multiphysics ic,
                        multiphysics_pointer update_multiphysics_rhs,
                        void* multiphysics_data = NULL);

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Input database
  SP_input d_input;
  /// Time-dependent material
  SP_material d_material;
  /// Mesh
  SP_mesh d_mesh;
  /// Quadrature
  SP_quadrature d_quadrature;
  /// Flag for moment or discrete problem
  bool d_discrete;
  /// Number of groups
  size_t d_number_groups;
  /// Include fission
  bool d_multiply;
  /// Synthetic source
  SP_syntheticsource d_syntheticsource;
  /// External sources.  The first external source is the synthetic source.
  vec_source d_sources;
  /// Fission source
  SP_fissionsource d_fissionsource;
  /// Working state.  This is initially assigned the initial condition.
  SP_state d_state;
  /// Working precursor vector.
  SP_precursors d_precursor;
  /// Working multiphysics vector
  SP_multiphysics d_multiphysics;
  /// Previous state iterate.
  SP_state d_state_0;
  /// Previous precursor iterate.
  SP_precursors d_precursor_0;
  /// Previous multiphysics iterate
  SP_multiphysics d_multiphysics_0;
  /// Time step size
  double d_dt;
  /// Step factor (for IMP)
  double d_step_factor;
  /// Final time
  double d_final_time;
  /// Number of time steps
  size_t d_number_steps;
  /// Integration scheme
  size_t d_scheme;
  /// Order of method.
  size_t d_order;
  /// Vector of previous states
  vec_states d_states;
  /// Vector of previous precursor concentrations
  vec_precursors d_precursors;
  /// Vector of previous physics iterates
  vec_multiphysics d_vec_multiphysics;
  /// Flag to write out time-dependent fluxes
  bool d_do_output;
  /// SILO output
  SP_silooutput d_silooutput;
  /// Fixed source solver
  SP_mg_solver d_solver;
  /// IMP negative flux fixup.
  bool d_fixup;
  /// No extrapolation flag (for first steps of high order BDF)
  bool d_no_extrapolation;
  /// Iteration count for nonlinear problems
  size_t d_iteration;
  /// Residual norm
  double d_residual_norm;
  /// Time step monitor
  monitor_pointer d_monitor;
  /// Monitor data
  void* d_monitor_data;
  /// Monitor level (only on or off right now)
  size_t d_monitor_level;
  /// Tolerance for nonlinear iterations
  double d_tolerance;
  /// Maximum nonlinear iterations
  size_t d_maximum_iterations;
  /// Multiphysics callback
  multiphysics_pointer d_update_multiphysics_rhs;
  /// Multiphysics data
  void* d_multiphysics_data;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  /// Given the initial fluxes, compute the precursor concentration.
  void initialize_precursors();

  /// Given the new state and previous precursor vectors, compute new precursor
  void update_precursors(const double t, const double dt, const size_t order);

  /// Given the new state and previous precursor vectors, compute new precursor
  void update_multiphysics(const double t, const double dt, const size_t order);


  /**
   *  @brief Cycle the states and precursors.
   *
   *  The first term gets a copy of the initial state.
   */
  void cycle_states_precursors(const size_t order);

  /// Update the synthetic and any external sources
  void update_sources(const double t, const double dt, const size_t order);

  /// Extrapolate for IMP
  void extrapolate();

  /// Perform a step
  void step(const double t,
            const double dt,
            const size_t order,
            const bool   flag);

  /// Check convergence
  bool check_convergence();

};

/**
 *  @brief Default time step monitor
 *  @param  data  Pointer to arbitrary user data
 *  @param  ts    TimeStepper pointer
 *  @param  step  Current time step
 *  @param  t     Time at step
 *  @param  dt    Step size
 *  @param  it    Iteration count (for nonlinear problems)
 */
template <class D>
void ts_default_monitor(void* data,
                        TimeStepper<D>* ts,
                        int step,
                        double t,
                        double dt,
                        int it,
                        bool converged);

} // end namespace detran

#endif /* detran_TIMESTEPPER_HH_ */
