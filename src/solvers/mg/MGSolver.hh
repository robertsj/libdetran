//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  MGSolver.hh
 *  @brief MGSolver class definition
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_MGSOLVER_HH_
#define detran_MGSOLVER_HH_

#include "solvers/Solver.hh"
#include "external_source/ExternalSource.hh"
#include "geometry/Mesh.hh"
#include "material/Material.hh"
#include "transport/FissionSource.hh"
#include "utilities/DBC.hh"
#include "utilities/InputDB.hh"
#include "utilities/SP.hh"

namespace detran
{

/**
 *  @class MGSolver
 *  @brief Base class for multigroup solvers.
 */
template <class D>
class MGSolver: public Solver<D>
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<MGSolver<D> >        SP_solver;
  typedef Solver<D>                                 Base;
  typedef typename Base::SP_input                   SP_input;
  typedef typename Base::SP_state                   SP_state;
  typedef typename Base::SP_mesh                    SP_mesh;
  typedef typename Base::SP_material                SP_material;
  typedef typename Base::SP_boundary                SP_boundary;
  typedef typename Base::SP_externalsource          SP_externalsource;
  typedef typename Base::vec_externalsource         vec_externalsource;
  typedef typename Base::SP_fissionsource           SP_fissionsource;
  typedef typename Base::size_t                     size_t;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param state             State vectors, etc.
   *  @param material          Material definitions.
   *  @param boundary          Boundary fluxes.
   *  @param q_e               Vector of user-defined external sources
   *  @param q_f               Fission source.
   *  @param multiply          Flag for a multiplying fixed source problem
   */
  MGSolver(SP_state                  state,
           SP_material               material,
           SP_boundary               boundary,
           const vec_externalsource &q_e,
           SP_fissionsource          q_f,
           bool                      multiply = false);

  /// Virtual destructor
  virtual ~MGSolver(){};

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL MULTIGROUP SOLVERS MUST IMPLEMENT
  //--------------------------------------------------------------------------//

  /// Solve the multigroup equations.
  virtual void solve(const double keff = 1.0) = 0;

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /// Return number of sweeps
  virtual int number_sweeps() const
  {
    return 0;
  }

protected:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Expose base members
  using Base::d_input;
  using Base::d_state;
  using Base::d_mesh;
  using Base::d_material;
  using Base::d_boundary;
  using Base::d_externalsources;
  using Base::d_fissionsource;
  using Base::d_number_groups;
  using Base::d_maximum_iterations;
  using Base::d_tolerance;
  using Base::d_print_level;
  using Base::d_print_interval;
  using Base::d_adjoint;

  /// Turn off upscatter
  bool d_downscatter;
  /// Flag for multiplying fixed source problem
  bool d_multiply;

};

} // namespace detran

#endif /* detran_MGSOLVER_HH_ */
