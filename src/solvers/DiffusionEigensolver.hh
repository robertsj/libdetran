//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   DiffusionEigensolver.hh
 *  @author robertsj
 *  @date   Jul 26, 2012
 *  @brief  DiffusionEigensolver class definition.
 */
//---------------------------------------------------------------------------//

#ifndef DIFFUSIONEIGENSOLVER_HH_
#define DIFFUSIONEIGENSOLVER_HH_

#include "DiffusionLossOperator.hh"
#include "DiffusionGainOperator.hh"
#include "callow/vector/Vector.hh"
#include "boundary/BoundaryDiffusion.hh"
#include "external_source/ExternalSource.hh"
#include "transport/State.hh"
#include "utilities/Definitions.hh"
#include "callow/solver/EigenSolver.hh"
#include "callow/solver/EigenSolverCreator.hh"

namespace detran
{

/**
 *  @class DiffusionEigensolver
 *  @brief Eigensolver for multigroup diffusion equations.
 */
template <class D>
class DiffusionEigensolver
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::InputDB::SP_input         SP_input;
  typedef State::SP_state                             SP_state;
  typedef detran_geometry::Mesh::SP_mesh              SP_mesh;
  typedef detran_material::Material::SP_material      SP_material;
  //
  typedef BoundaryDiffusion<D>                        Boundary_T;
  typedef typename Boundary_T::SP_boundary            SP_boundary;
  //
  typedef detran_external_source::
          ExternalSource::SP_externalsource           SP_source;
  //
  typedef State::moments_type                         moments_type;
  //
  typedef DiffusionLossOperator::SP_lossoperator        SP_lossoperator;
  typedef DiffusionGainOperator::SP_gainoperator        SP_gainoperator;
  typedef callow::Vector<double>                        Vector_T;
  typedef Vector_T::SP_vector                           SP_vector;
  typedef callow::EigenSolverCreator<double>            Creator_T;
  typedef Creator_T::SP_solver                          SP_solver;
  typedef callow::EigenSolver<double>::SP_linearsolver  SP_linearsolver;
  //
  typedef detran_utilities::size_t                    size_t;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param input      parameter database
   *  @param material   material database
   *  @param mesh       mesh definition
   *  @param state      state vector
   */
  DiffusionEigensolver(SP_input    input,
                       SP_material material,
                       SP_mesh     mesh,
                       SP_state    state);

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /// Solve the system
  void solve();

  /// @name Getters
  /// @{
  SP_input input() const { return d_input; }
  SP_material material() const { return d_material; }
  SP_mesh mesh() const { return d_mesh; }
  SP_state state() const { return d_state; }
  SP_boundary boundary() const { return d_boundary; }
  SP_lossoperator lossoperator() const { return d_M; }
  SP_gainoperator gainoperator() const { return d_F; }
  SP_solver eigensolver() const { return d_solver; }
  SP_solver linearsolver() const { return d_solver->linearsolver(); }
  SP_vector phi() { return d_phi; }
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
  /// Boundary container
  SP_boundary d_boundary;
  /// Loss operator
  SP_lossoperator d_M;
  /// Fission operator
  SP_gainoperator d_F;
  /// Unknown vector
  SP_vector d_phi;
  /// Working vector
  SP_vector d_work;
  /// Number of unknowns
  size_t d_problem_size;
  /// Solver type (e.g. power)
  std::string d_solver_type;
  /// Maximum linear solver iterations
  size_t d_maximum_iterations;
  /// Tolerance for linear solve
  double d_tolerance;
  /// Pointer to eigensolver
  SP_solver d_solver;
  /// Solve in adjoint mode?
  bool d_adjoint;

};

} // end namespace detran

#endif /* DIFFUSIONEIGENSOLVER_HH_ */
