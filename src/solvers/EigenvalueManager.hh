//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   EigenvalueManager.hh
 *  @author robertsj
 *  @date   Oct 24, 2012
 *  @brief  EigenvalueManager class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_EIGENVALUEMANAGER_HH_
#define detran_EIGENVALUEMANAGER_HH_

#include "Eigensolver.hh"

namespace detran
{

/**
 *  @class EigenvalueManager
 *  @brief Manage solution of a multigroup eigenvalue problem
 */
template <class D>
class EigenvalueManager
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef typename Eigensolver<D>::Fixed_T      Fixed_T;
  typedef typename Eigensolver<D>::SP_solver    SP_solver;
  typedef typename Fixed_T::SP_manager          SP_mg_solver;
  typedef typename Fixed_T::SP_input            SP_input;
  typedef typename Fixed_T::SP_state            SP_state;
  typedef typename Fixed_T::SP_mesh             SP_mesh;
  typedef typename Fixed_T::SP_material         SP_material;
  typedef typename Fixed_T::SP_quadrature       SP_quadrature;
  typedef typename Fixed_T::SP_boundary         SP_boundary;
  typedef typename Fixed_T::SP_fissionsource    SP_fissionsource;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param input      parameter database
   *  @param material   material database
   *  @param mesh       mesh definition
   */
  EigenvalueManager(SP_input    input,
                    SP_material material,
                    SP_mesh     mesh);

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /**
   *  @brief Solve the system
   *
   *  By changing the appropriate database parameters, a problem already
   *  set up can be solved by a different method by calling this
   *  method.
   *
   *  @param keff   Scaling factor for multiplying problems
   */
  bool solve();

  /// @name Getters
  /// @{
  SP_input input() const { return d_mg_solver->input(); }
  SP_material material() const { return d_mg_solver->material(); }
  SP_mesh mesh() const { return d_mg_solver->mesh(); }
  SP_state state() const { return d_mg_solver->state(); }
  SP_boundary boundary() const { return d_mg_solver->boundary(); }
  SP_quadrature quadrature() const { return d_mg_solver->quadrature(); }
  SP_fissionsource fissionsource() const { return d_mg_solver->fissionsource(); }
  /// @}

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Eigensolver
  SP_solver d_solver;
  /// Multigroup solver
  SP_mg_solver d_mg_solver;
  /// Adjoint mode flag
  bool d_adjoint;
  /// Discretization general type (SN, MOC, or diffusion)
  int d_discretization;
  /// Problem setup status flag
  bool d_is_setup;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//


};

} // end namespace detran


#endif /* detran_EIGENVALUEMANAGER_HH_ */
