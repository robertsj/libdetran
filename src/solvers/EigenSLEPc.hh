//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   EigenSLEPc.hh
 * \author robertsj
 * \date   Jun 18, 2012
 * \brief  EigenSLEPc class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

#ifndef EIGENSLEPC_HH_
#define EIGENSLEPC_HH_

// Detran
#include "Eigensolver.hh"

// System
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

  typedef SP<EigenSLEPc<D> >                    SP_solver;
  typedef Eigensolver<D>                        Base;
  typedef typename Base::SP_solver              SP_base;
  typedef typename
      MultigroupSolver<D>::SP_solver            SP_mg_solver;
  // basic objects
  typedef InputDB::SP_input                     SP_input;
  typedef State::SP_state                       SP_state;
  typedef Mesh::SP_mesh                         SP_mesh;
  typedef Material::SP_material                 SP_material;
  typedef Quadrature::SP_quadrature             SP_quadrature;
  typedef typename BoundaryBase<D>::SP_boundary SP_boundary;
  typedef FissionSource::SP_source              SP_fissionsource;

  /*!
   *  \brief Constructor
   *
   *  \param input             Input database.
   *  \param state             State vectors, etc.
   *  \param mesh              Problem mesh.
   *  \param mat               Material definitions.
   *  \param quadrature        Angular mesh.
   *  \param boundary          Boundary fluxes.
   *  \param fission_source    Fission source.
   */
  EigenSLEPc(SP_input           input,
             SP_state           state,
             SP_mesh            mesh,
             SP_material        material,
             SP_quadrature      quadrature,
             SP_boundary        boundary,
             SP_fissionsource   q_f);

  /// SP Constructor
  static SP<EigenSLEPc<D> >
  Create(SP<detran::InputDB>          input,
         SP<detran::State>            state,
         SP<detran::Mesh>             mesh,
         SP<detran::Material>         material,
         SP<detran::Quadrature>       quadrature,
         SP<detran::BoundaryBase<D> > boundary,
         SP<detran::FissionSource>    q_f)
  {
    SP_solver p(new EigenSLEPc(input, state, mesh, material,
                               quadrature, boundary, q_f));
    return p;
  }

  /// Solve the eigenvalue problem.
  void solve();

  /// Unimplemented DBC function.
  bool is_valid() const
  {
    return true;
  }


protected:

  using Base::b_input;
  using Base::b_state;
  using Base::b_mesh;
  using Base::b_material;
  using Base::b_quadrature;
  using Base::b_boundary;
  using Base::b_fissionsource;
  using Base::b_mg_solver;
  using Base::b_max_iters;
  using Base::b_tolerance;
  using Base::b_print_out;
  using Base::b_print_interval;

  /// \name Protected Data
  /// \{

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

  /// \}

  /// \name Implementation
  /// \{

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

  /// \}

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
