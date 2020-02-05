//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Solver.hh
 *  @brief  Solver class definition
 *  @author Jeremy Roberts
 *  @date   Sep 8, 2012
 */
//---------------------------------------------------------------------------//

#ifndef detran_SOLVER_HH_
#define detran_SOLVER_HH_

#include "boundary/BoundaryBase.hh"
#include "external_source/ExternalSource.hh"
#include "geometry/Mesh.hh"
#include "material/Material.hh"
#include "transport/FissionSource.hh"
#include "utilities/DBC.hh"
#include "utilities/Definitions.hh"
#include "utilities/InputDB.hh"

namespace detran
{

/**
 *  @class Solver
 *  @brief Boilerplate shared by all solvers.
 */
template <class D>
class Solver
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::InputDB::SP_input         SP_input;
  typedef State::SP_state                             SP_state;
  typedef detran_geometry::Mesh::SP_mesh              SP_mesh;
  typedef detran_material::Material::SP_material      SP_material;
  typedef typename BoundaryBase<D>::SP_boundary       SP_boundary;
  typedef detran_external_source::
          ExternalSource::SP_externalsource           SP_externalsource;
  typedef detran_external_source::
          ExternalSource::vec_externalsource          vec_externalsource;
  typedef FissionSource::SP_fissionsource             SP_fissionsource;
  typedef detran_utilities::size_t                    size_t;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *
   *  @param input             Input database.
   *  @param state             State vectors, etc.
   *  @param mesh              Problem mesh.
   *  @param material          Material definitions.
   *  @param boundary          Boundary fluxes.
   *  @param q_e               Vector of user-defined external sources
   *  @param q_f               Fission source.
   *  @param multiply          Flag for a multiplying fixed source problem
   */
  Solver(SP_state                   state,
         SP_material                material,
         SP_boundary                boundary,
         const vec_externalsource  &q_e,
         SP_fissionsource           q_f);

  /// Pure virtual destructor
  virtual ~Solver() = 0;

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /// Reset the tolerance.
  void set_tolerance(double tol)
  {
    Require(tol > 0.0);
    d_tolerance = tol;
  }

  /// Reset the maximum iterations.
  void set_max_iters(int max_iters)
  {
    Require(max_iters > 0);
    d_maximum_iterations = max_iters;
  }

  /**
   *  @brief  Refresh the solver
   *
   *  This can be reimplemented by a subclass if internal structures must
   *  be updated due to some change (e.g. material perturbation)
   */
  virtual void refresh()
  {
    /* ... */
  }

  /// Getters
  //@}
  SP_input input() {return d_input;}
  SP_state state() {return d_state;}
  SP_mesh mesh() {return d_mesh;}
  SP_material material() {return d_material;}
  SP_boundary boundary() {return d_boundary;}
  bool adjoint() {return d_adjoint;}
  //@}

protected:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// User input.
  SP_input d_input;
  /// State vectors.
  SP_state d_state;
  /// Problem mesh
  SP_mesh d_mesh;
  /// Materials definitions.
  SP_material d_material;
  /// Boundary fluxes.
  SP_boundary d_boundary;
  /// Vector of user-defined external sources
  vec_externalsource d_externalsources;
  /// Fission source, if used
  SP_fissionsource d_fissionsource;
  /// Number of groups
  size_t d_number_groups;
  /// Maximum iterations
  size_t d_maximum_iterations;
  /// Convergence tolerance
  double d_tolerance;
  /// Print out flag
  int d_print_level;
  /// Interval for print out
  int d_print_interval;
  /// Adjoint flag
  bool d_adjoint;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  // Default constructor for inherited class use
  Solver()
    : d_maximum_iterations(100)
    , d_tolerance(1e-5)
    , d_print_level(2)
    , d_print_interval(10)
    , d_adjoint(false)
  {/* ... */}

};

} // end namespace detran

#endif // SOLVER_HH_ 

//---------------------------------------------------------------------------//
//              end of file Solver.hh
//---------------------------------------------------------------------------//
