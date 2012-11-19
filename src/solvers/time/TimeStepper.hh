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
#include "kinetics/TimeDependentMaterial.hh"
#include "kinetics/TimeDependentExternalSource.hh"
#include "kinetics/SyntheticSource.hh"
#include "utilities/DBC.hh"
#include "utilities/Definitions.hh"
#include "utilities/SP.hh"
#include "transport/State.hh"
#include "kinetics/BDFCoefficients.hh"

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
  typedef SyntheticSource::SP_precursors                  SP_precursors;
  typedef SyntheticSource::vec_precursors                 vec_precursors;
  typedef FissionSource::SP_fissionsource                 SP_fissionsource;
  typedef detran_utilities::vec_int                       vec_int;
  typedef detran_ioutils::SiloOutput::SP_silooutput       SP_silooutput;

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

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /// Add a time-dependent external source
  void add_source(SP_tdsource source);

  /// Solve
  void solve(SP_state initial_state);

  SP_state state() {return d_state;}

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
  /// Number of previous steps kept
  size_t d_number_previous;
  /// Order of method.
  size_t d_order;

  vec_states d_states;
  vec_precursors d_precursors;

  /// Flag to write out time-dependent fluxes
  bool d_do_output;
  SP_silooutput d_silooutput;

  SP_mg_solver d_solver;

  /// IMP negative flux fixup.
  bool d_fixup;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  /// Given the initial fluxes, compute the precursor concentration.
  void initialize_precursors();

  /// Given the new state and previous precursor vectors, compute new precursor
  void update_precursors();

  /**
   *  @brief Cycle the states and precursors.
   *
   *  The first term gets a copy of the initial state.
   */
  void cycle_states_precursors();

  /// Update the synthetic and any external sources
  void update_sources(const double t, const double dt);

  /// Extrapolate for IMP
  void extrapolate();

};

} // end namespace detran_kinetics

#endif /* detran_TIMESTEPPER_HH_ */
