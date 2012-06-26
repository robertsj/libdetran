//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   InnerGMRES.hh
 * \author robertsj
 * \date   Apr 4, 2012
 * \brief  InnerGMRES class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

#ifndef INNERGMRES_HH_
#define INNERGMRES_HH_

// Detran
#include "InnerIteration.hh"

// System
#include "petsc.h"

namespace detran
{

//---------------------------------------------------------------------------//
/*!
 *  \class InnerGMRES
 *  \brief Solve the within-group problem with GMRES.
 *
 *  The within group equation can be written
 *  \f[
 *      (\mathbf{I} - \mathbf{D}\mathbf{L}^{-1}\mathbf{MS})\phi
 *      = \mathbf{D} \mathbf{L}^{-1} Q \, .
 *  \f]
 *  This class couples with PETSc to make available its set of applicable
 *  solvers, the default being GMRES.  Other solvers are selected
 *  by command line flags, e.g. -ksp_type bcgs uses a BiCongugate Gradient
 *  Stabilized algorithm.  Past experience suggests GMRES works best.
 *
 *  Note, all iterative solvers do better with preconditioning, which is
 *  current not implemented.
 *
 */
//---------------------------------------------------------------------------//

template <class D>
class InnerGMRES: public InnerIteration<D>
{

public:

  typedef SP<InnerGMRES>                        SP_inner;
  typedef InnerIteration<D>                     Base;
  typedef typename InnerIteration<D>::SP_inner  SP_base;
  // basic objects
  typedef InputDB::SP_input                     SP_input;
  typedef State::SP_state                       SP_state;
  typedef Mesh::SP_mesh                         SP_mesh;
  typedef Material::SP_material                 SP_material;
  typedef Quadrature::SP_quadrature             SP_quadrature;
  typedef typename BoundaryBase<D>::SP_boundary SP_boundary;
  typedef typename MomentToDiscrete<D>::SP_MtoD SP_MtoD;
  // source typedefs
  typedef ExternalSource::SP_source             SP_externalsource;
  typedef FissionSource::SP_source              SP_fissionsource;
  //
  typedef typename Sweeper<D>::SP_sweeper       SP_sweeper;
  typedef typename
      SweepSource<D>::SP_sweepsource            SP_sweepsource;
  //
  typedef State::moments_type                   moments_type;

  /*!
   *  \brief Constructor
   *
   *  \param input             Input database.
   *  \param state             State vectors, etc.
   *  \param mesh              Problem mesh.
   *  \param mat               Material definitions.
   *  \param quadrature        Angular mesh.
   *  \param boundary          Boundary fluxes.
   *  \param external_source   User-defined external source.
   *  \param fission_source    Fission source.
   */
  InnerGMRES(SP_input           input,
             SP_state           state,
             SP_mesh            mesh,
             SP_material        material,
             SP_quadrature      quadrature,
             SP_boundary        boundary,
             SP_externalsource  q_e,
             SP_fissionsource   q_f);

  /// SP Constructor
  static SP<InnerGMRES<D> >
  Create(SP<detran::InputDB>          input,
         SP<detran::State>            state,
         SP<detran::Mesh>             mesh,
         SP<detran::Material>         material,
         SP<detran::Quadrature>       quadrature,
         SP<detran::BoundaryBase<D> > boundary,
         SP<detran::ExternalSource>   q_e,
         SP<detran::FissionSource>    q_f)
  {
    SP_inner p(new InnerGMRES(input, state, mesh, material,
                              quadrature, boundary, q_e, q_f));
    return p;
  }

  /// Solve the within group equation.
  void solve(int g);

  // Friend functions for applying dimension-specific operator.
  friend PetscErrorCode apply_WGTO_1D(Mat A, Vec x, Vec y);
  friend PetscErrorCode apply_WGTO_2D(Mat A, Vec x, Vec y);
  friend PetscErrorCode apply_WGTO_3D(Mat A, Vec x, Vec y);

private:

  // Make inherited data visible
  using Base::d_input;
  using Base::d_state;
  using Base::d_mesh;
  using Base::d_material;
  using Base::d_quadrature;
  using Base::d_boundary;
  using Base::d_sweeper;
  using Base::d_sweepsource;
  using Base::d_tolerance;
  using Base::d_max_iters;
  using Base::d_print_out;
  using Base::d_print_interval;
  using Base::d_g;
  //using Base::b_acceleration;

  /// \name Private Data
  /// \{

  /// Main linear solver
  KSP d_solver;

  /// Operator "A" in "Ax = b"
  Mat d_operator;

  /// Solution vector
  Vec d_X;

  /// Right hand side
  Vec d_B;

  /// Size of the moments portion of d_X
  int d_moments_size;

  /// Size of the boundary portion of d_X
  int d_boundary_size;

  int d_reflective_solve_iterations;

  /// \}

  /// \name Implementation
  /// \{

  /// Set the templated operator function.
  PetscErrorCode set_operation();

  /// Build the right hand side.
  void build_rhs(State::moments_type &B);

  //---------------------------------------------------------------------------//
  /*!
   * \brief A matrix-vector shell for the within-group transport operator.
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
  PetscErrorCode apply_WGTO(Mat A, Vec x, Vec y);

  /// \}

};

} // namespace detran

//---------------------------------------------------------------------------//
// EXTERNAL WRAPPER FUNCTIONS
//---------------------------------------------------------------------------//

PetscErrorCode apply_WGTO_1D(Mat A, Vec x, Vec y);
PetscErrorCode apply_WGTO_2D(Mat A, Vec x, Vec y);
PetscErrorCode apply_WGTO_3D(Mat A, Vec x, Vec y);

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "InnerGMRES.i.hh"

#endif /* INNERGMRES_HH_ */

//---------------------------------------------------------------------------//
//              end of InnerGMRES.hh
//---------------------------------------------------------------------------//
