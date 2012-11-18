//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   TimeDependentManager.hh
 *  @brief  TimeDependentManager
 *  @author Jeremy Roberts
 *  @date   Nov 16, 2012
 */
//---------------------------------------------------------------------------//

#ifndef detran_TIMEDEPENDENTMANAGER_HH_
#define detran_TIMEDEPENDENTMANAGER_HH_

namespace detran
{

/**
 *  @class TimeDependentManager
 *  @brief Manage solution of time-dependnent problems
 */
template <class D>
class TimeDependentManager
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<TimeDependentManager<D> >  SP_manager;
  typedef detran_utilities::InputDB::SP_input             SP_input;
  typedef State::SP_state                                 SP_state;
  typedef detran_geometry::Mesh::SP_mesh                  SP_mesh;
  typedef detran_material::Material::SP_material          SP_material;
  typedef detran_angle::Quadrature::SP_quadrature         SP_quadrature;
  typedef BoundaryBase<D>                                 Boundary_T;
  typedef typename Boundary_T::SP_boundary                SP_boundary;
  typedef detran_external_source::
          ExternalSource::SP_externalsource               SP_source;
  typedef detran_external_source::
          ExternalSource::vec_externalsource              vec_source;
  typedef FissionSource::SP_fissionsource                 SP_fissionsource;
  typedef State::moments_type                             moments_type;
  typedef typename MGSolver<D>::SP_solver                 SP_solver;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param input      parameter database
   *  @param material   material database
   *  @param mesh       mesh definition
   */
  TimeDependentManager(SP_input    input,
                       SP_material material,
                       SP_mesh     mesh);

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /**
   *  @brief Sets up a problem to be solved
   *
   *  This sets up the discretization, quadrature, boundary, and
   *  state.  By changing the appropriate database parameters,
   *  a call to setup will rebuild the problem using a different
   *  discretization, etc.
   *
   *  Source defining external sources, the quadrature is often needed.
   *  By calling setup first, the client can extract the quadrature for
   *  building an external source.
   *
   */
  void setup();

  /**
   *  @brief Add a new external source.
   *
   *  This is called after setup, since external sources often
   *  need the quadrature set produced by the manager.
   */
  void set_source(SP_source q);

  /**
   *  @brief Set the solver based on parameter database.
   *
   *  This must be called after setup is called and after all sources
   *  have been set.  By changing the appropriate parameters, a
   *  problem can be re-run with a different solver.
   */
  bool set_solver();

  /**
   *  @brief Solve the system
   *  @param keff   Scaling factor for multiplying problems
   */
  bool solve(const double keff = 1.0);

  /// @name Getters
  /// @{
  SP_input input() const { return d_input; }
  SP_material material() const { return d_material; }
  SP_mesh mesh() const { return d_mesh; }
  SP_state state() const { return d_state; }
  SP_boundary boundary() const { return d_boundary; }
  SP_quadrature quadrature() const { return d_quadrature; }
  SP_fissionsource fissionsource() const { return d_fissionsource; }
  int discretization() const { return d_discretization; }
  SP_solver solver() const { return d_solver; }
  /// @}

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// User input
  SP_input d_input;
  /// Material database
  SP_material d_material;
  /// Mesh
  SP_mesh d_mesh;
  /// State vector
  SP_state d_state;
  /// Quadrature
  SP_quadrature d_quadrature;
  /// Boundary container
  SP_boundary d_boundary;
  /// Fission source
  SP_fissionsource d_fissionsource;
  /// External volume sources
  vec_source d_sources;
  /// Multigroup solver
  SP_solver d_solver;
  /// Adjoint mode flag
  bool d_adjoint;
  /// Discretization general type (SN, MOC, or diffusion)
  int d_discretization;
  /// Flag for multiplying fixed source problem
  bool d_multiply;
  /// Flag for including fission source
  bool d_fission;
  /// Problem setup status flag
  bool d_is_setup;
  /// Solver setup status flag
  bool d_is_ready;
  /// Iteration generation
  int d_generation;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//


};

} // end namespace detran

#endif // detran_TIMEDEPENDENTMANAGER_HH_

//---------------------------------------------------------------------------//
//              end of file TimeDependentManager.hh
//---------------------------------------------------------------------------//
