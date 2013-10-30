//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  FixedSourceManager.hh
 *  @brief FixedSourceManager class definition
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_FIXEDSOURCEMANAGER_HH_
#define detran_FIXEDSOURCEMANAGER_HH_

#include "solvers/solvers_export.hh"
#include "TransportManager.hh"
#include "solvers/mg/MGSolver.hh"
#include "angle/Quadrature.hh"

namespace detran
{

/**
 *  @class FixedSourceManager
 *  @brief Manage solution of a multigroup fixed source problem
 *
 *  Fixed source problems are classified generically as "fixed" or
 *  "multiplying".  For the latter, the fission source is included
 *  based on an assumed k-eigenvalue (defaulted to unity).  A
 *  multiplying problem can be solved with fission treated like
 *  scatter, or it can be added via iteration outside the normal
 *  multigroup solver.  The latter is useful when we want to
 *  expand solutions in fission generation series.
 */
template <class D>
class FixedSourceManager: TransportManager
{

public:

  //--------------------------------------------------------------------------//
  // ENUMERATION
  //--------------------------------------------------------------------------//

  // basic spatial discretization categories
  enum EQTYPES
  {
    SN, MOC, DIFF, END_EQTYPES
  };

  // various fixed source problem approaches
  enum FIXEDTYPE
  {
    FIXED, MULTIPLY, END_FIXEDTYPE
  };

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<FixedSourceManager<D> >  SP_manager;
  typedef detran_utilities::InputDB::SP_input           SP_input;
  typedef State::SP_state                               SP_state;
  typedef detran_geometry::Mesh::SP_mesh                SP_mesh;
  typedef detran_material::Material::SP_material        SP_material;
  typedef detran_angle::Quadrature::SP_quadrature       SP_quadrature;
  typedef BoundaryBase<D>                               Boundary_T;
  typedef typename Boundary_T::SP_boundary              SP_boundary;
  typedef detran_external_source::
          ExternalSource::SP_externalsource             SP_source;
  typedef detran_external_source::
          ExternalSource::vec_externalsource            vec_source;
  typedef FissionSource::SP_fissionsource               SP_fissionsource;
  typedef State::moments_type                           moments_type;
  typedef typename MGSolver<D>::SP_solver               SP_solver;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param argc       command line count
   *  @param argv       command line values
   *  @param input      parameter database
   *  @param material   material database
   *  @param mesh       mesh definition
   *  @param multiply   flag for multiplying fixed source problem
   *  @param fission    ensure a fission source is built (e.g. for eigen)
   */
  FixedSourceManager(int         argc,
                     char       *argv[],
                     SP_input    input,
                     SP_material material,
                     SP_mesh     mesh,
                     bool        multiply = false,
                     bool        fission = false);

  /// Constructor (without command line)
  FixedSourceManager(SP_input    input,
                     SP_material material,
                     SP_mesh     mesh,
                     bool        multiply = false,
                     bool        fission = false);

  /// Virtual destructor
  virtual ~FixedSourceManager(){}

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

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

  /**
   *  @brief Perform a fission iteration
   *
   *  This function provides a way for a client to perform expansion
   *  of a fixed source problem solution in terms of a series in
   *  fission generations.
   *
   *  To use this function, the client sets up a problem in the typical
   *  way by setting a source, options, etc.  For each generation, the
   *  client iterates.  The function returns the norm of the difference
   *  between fluxes of successive generations.  Following each iteration,
   *  the state is filled with the m-th iteration contribution to the
   *  flux.  In other words, the solution is represented as
   *  @f[
   *      \phi = \phi_0 + \phi_1 + \ldots
   *  @f]
   *  where @f$ \phi_0 @f$ is the flux due to the initial source,
   *        @f$ \phi_1 @f$ is the "once-fissioned" flux, and so on.
   *  Thus, for arbitrary eigenvalue @f$ k @f$, we have
   *  @f[
   *      \phi = \phi_0 + \frac{1}{k}\phi_1 + \frac{1}{k^2}\phi_2 + \ldots
   *  @f]
   *  Moreover, it turns out that the ratio of successive terms tends
   *  toward @f$ k @f$ of the domain in a vacuum, @f$ k_v @f$, so that the
   *  series can be approximately summed as
   *  @f[
   *      \phi = \phi_0 + \frac{1}{k}\phi_1 + \frac{1}{k^2(1-k_v/k)}\phi_2
   *  @f]
   *
   *  Note, after generation zero, the fixed sources are eliminated.
   */
  double iterate(const int generation);

  /// Update operators, etc.
  void update()
  {
    Require(d_solver);

    // Refresh the solver
    d_solver->refresh();
  }

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
  int number_sweeps() const { return d_solver->number_sweeps(); }
  bool adjoint() const {return d_adjoint;}
  /// @}

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

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

};

} // end namespace detran

#endif /* detran_FIXEDSOURCEMANAGER_HH_ */

//----------------------------------------------------------------------------//
//              end of FixedSourceManager.hh
//----------------------------------------------------------------------------//
