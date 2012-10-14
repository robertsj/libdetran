//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   FixedSourceManager.hh
 *  @author robertsj
 *  @date   Sep 25, 2012
 *  @brief  FixedSourceManager class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_FIXEDSOURCEMANAGER_HH_
#define detran_FIXEDSOURCEMANAGER_HH_

#include "angle/Quadrature.hh"
#include "boundary/BoundaryBase.hh"
#include "callow/solver/LinearSolver.hh"
#include "external_source/ExternalSource.hh"
#include "geometry/Mesh.hh"
#include "material/Material.hh"
#include "utilities/Definitions.hh"
#include "utilities/InputDB.hh"
#include "transport/FissionSource.hh"
#include "transport/State.hh"
#include <vector>

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
class FixedSourceManager
{

public:

  //-------------------------------------------------------------------------//
  // ENUMERATION
  //-------------------------------------------------------------------------//

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

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::InputDB::SP_input         SP_input;
  typedef State::SP_state                             SP_state;
  typedef detran_geometry::Mesh::SP_mesh              SP_mesh;
  typedef detran_material::Material::SP_material      SP_material;
  typedef detran_angle::Quadrature::SP_quadrature     SP_quadrature;
  //
  typedef BoundaryBase<D>                             Boundary_T;
  typedef typename Boundary_T::SP_boundary            SP_boundary;
  //
  typedef detran_external_source::
          ExternalSource::SP_externalsource           SP_source;
  typedef FissionSource::SP_fissionsource             SP_fissionsource;
  //
  typedef State::moments_type                         moments_type;
  //
  typedef callow::LinearSolver::SP_solver             SP_solver;
  //

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param input      parameter database
   *  @param material   material database
   *  @param mesh       mesh definition
   */
  FixedSourceManager(SP_input    input,
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
   */
  void setup();

  /// Add a new external source
  void set_source(SP_source q);

  /**
   *  @brief Solve the system
   *
   *  By changing the appropriate database parameters, a problem already
   *  set up can be solved by a different method by calling this
   *  method.
   *
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
   *  user iterates.  The function returns the norm of the difference
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

  /// @name Getters
  /// @{
  SP_input input() const { return d_input; }
  SP_material material() const { return d_material; }
  SP_mesh mesh() const { return d_mesh; }
  SP_state state() const { return d_state; }
  SP_boundary boundary() const { return d_boundary; }
  SP_quadrature quadrature() const { return d_quadrature; }
  SP_fissionsource fissionsource() const { return d_fissionsource; }
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
  std::vector<SP_source> d_sources;
  /// Linear solver
  SP_solver d_solver;
  /// Adjoint mode flag
  bool d_adjoint;
  /// Discretization general type (SN, MOC, or diffusion)
  int d_discretization;
  /// Fixed type
  int d_fixed_type;
  /// Setup status flag
  bool d_is_setup;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//


};

} // end namespace detran

#endif /* detran_FIXEDSOURCEMANAGER_HH_ */
