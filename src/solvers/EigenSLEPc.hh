//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   EigenSLEPc.hh
 * \author robertsj
 * \date   Jun 18, 2012
 * \brief  EigenSLEPc class definition.
 */
//---------------------------------------------------------------------------//

#ifndef EIGENSLEPC_HH_
#define EIGENSLEPC_HH_

#include "Eigensolver.hh"
#include "slepceps.h"

namespace detran
{

//---------------------------------------------------------------------------//
/*!
 * \class EigenSLEPc
 * \brief Solves the eigenvalue problem using SLEPc.
 *
 * The eigenvalue problem can be cast in the form
 * \f[
 *     \mathbf{A}d = kd \,
 * \f]
 * where \f$ d \f$ is the fission density and \f$ k \f$ is the
 * eigenvalue.  See Eigensolver for more details on this formulation.
 *
 * SLEPc is a package for solving eigenvalue problems and is built on
 * top of PETSc.  SLEPc offers a number of built in solvers, including
 * the explicitly-restarted Arnoldi method (ERAM), the Krylov-Schur (KS)
 * method, and power iteration (PI).  The default is KS, and other
 * solvers and options are available from the command line.
 *
 * Note, as mentioned in Eigensolver, Krylov methods require an
 * adequately-converged multigroup solve.  The paper by Warsa et al.
 * talks more about this, and it seems the multigroup
 * convergence is most
 * important in the first several iterations.  Thus, it might
 * be worth implementing a dynamic tolerance at some point.
 *
 */
//---------------------------------------------------------------------------//

template <class D>
class EigenSLEPc: public Eigensolver<D>
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef Eigensolver<D>                            Base;
  typedef typename Base::SP_solver                  SP_solver;
  typedef typename Base::SP_mg_solver               SP_mg_solver;
  typedef typename Base::SP_input                   SP_input;
  typedef typename Base::SP_state                   SP_state;
  typedef typename Base::SP_mesh                    SP_mesh;
  typedef typename Base::SP_material                SP_material;
  typedef typename Base::SP_boundary                SP_boundary;
  typedef typename Base::SP_fissionsource           SP_fissionsource;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param mg_solver         Multigroup solver
   */
  EigenSLEPc(SP_mg_solver mg_solver);

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL EIGENSOLVERS MUST IMPLEMENT
  //-------------------------------------------------------------------------//

  /// Solve the eigenvalue problem.
  void solve();

protected:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Expose base members
  using Base::d_input;
  using Base::d_state;
  using Base::d_mesh;
  using Base::d_material;
  using Base::d_fissionsource;
  using Base::d_number_groups;
  using Base::d_maximum_iterations;
  using Base::d_tolerance;
  using Base::d_print_level;
  using Base::d_print_interval;
  using Base::d_adjoint;
  using Base::d_mg_solver;

  /// SLEPc eigensolver.
  EPS d_solver;
  /// Operator.
  Mat d_operator;
  /// PETSc vector for fission density.
  Vec d_density;
  /// System size
  int d_size;
  /// Number of multigroup solves.
  int d_mg_solves;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  /// Set the templated operator function.
  PetscErrorCode set_operation();

  //---------------------------------------------------------------------------//
  /*!
   * \brief A matrix-vector shell for the eigenvalue problem.
   *
   * This is called by thin wrappers, since PETSc needs the matrix-vector
   * operation as a function pointer, which precludes a member function.
   *
   * \note This is public, since I haven't figured out a way to use friends
   *       in this way---maybe there is no way.
   *
   * \param   A       PETSc shell matrix
   * \param   x       Incoming PETSc vector
   * \param   y       Outgoing PETSc vector
   */
public:
  PetscErrorCode apply_eigen(Mat A, Vec x, Vec y);

};

} // namespace detran

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "EigenSLEPc.i.hh"

//---------------------------------------------------------------------------//
// EXTERNAL WRAPPER FUNCTIONS
//---------------------------------------------------------------------------//

PetscErrorCode apply_eigen_1D(Mat A, Vec x, Vec y);
PetscErrorCode apply_eigen_2D(Mat A, Vec x, Vec y);
PetscErrorCode apply_eigen_3D(Mat A, Vec x, Vec y);

#endif /* EIGENSLEPc_HH_ */
