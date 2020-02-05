//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Eigensolver.hh
 *  @brief Eigensolver class definition
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_EIGENSOLVER_HH_
#define detran_EIGENSOLVER_HH_

#include "solvers/FixedSourceManager.hh"

namespace detran
{

//----------------------------------------------------------------------------//
/**
 *  @class Eigensolver
 *  @brief Base class for solving the eigenvalue problem
 *
 *  @section theeigenvalueproblem The Eigenvalue Problem
 *
 *  The steady-state balance of neutrons in a fissile system is
 *  characterized by the multigroup eigenvalue problem
 *
 *  @f[
 *      \hat{\Omega} \cdot \nabla  \psi +
 *        \Sigma_t(\vec{r},E_g) \psi(\vec{r},\hat{\Omega},E_g) =
 *      \frac{1}{4\pi} \sum^G_{g'} \int_{4\pi} d\Omega'
 *        \Sigma_S(\vec{r},\hat{\Omega}'\cdot \hat{\Omega}, E_g'\to E_g)
 *          \psi(\vec{r},\hat{\Omega}',E_{g'}) +
 *      \frac{\chi_g}{4\pi k} \sum^G_{g'} \int_{4\pi} d\Omega'
 *        \nu \Sigma_f(\vec{r}, E_{g'}) \psi(\vec{r}, \hat{\Omega}', E_{g'})
 *  @f]
 *
 *  where \f$ k \f$ is the eigenvalue, which physically represents the
 *  ratio of the number of neutrons in successive generations.
 *
 *  @section operatorform Operator Form
 *
 *  In operator form, the eigenvalue problem is
 *
 *  @f[
 *      (\mathbf{I} - \mathbf{DL}^{-1}\mathbf{MS})\phi =
 *        \frac{1}{k} \mathbf{DL}^{-1}\mathbf{M}\chi \mathbf{F}^T \phi \, .
 *  @f]
 *
 *  where the eigenvector consists of the angular flux moments.  In the more
 *  standard form \f$ \mathbf{A}x = k x \f$, we have two options.  We
 *  can keep the problem in terms of the multigroup fluxes so that
 *
 *  @f[
 *      \overbrace{(\mathbf{I} - \mathbf{DL}^{-1}\mathbf{MS})^{-1}
 *        \mathbf{DL}^{-1}\mathbf{M}\chi \mathbf{F}^T}^{\mathbf{A}} \phi
 *       = k \phi \, .
 *  @f]
 *
 *  Alternatively, we can solve in terms of the fission density,
 *  @f[
 *      d = \mathbf{F}^T \phi = \sum^G_g \nu \Sigma_{fg} \phi_g
 *  @f]
 *  so that
 *  @f[
 *      \overbrace{\mathbf{F}^T (\mathbf{I} - \mathbf{DL}^{-1}\mathbf{MS})^{-1}
 *        \mathbf{DL}^{-1}\mathbf{M}\chi }^{\mathbf{A}_{\text{ei}}} d =
 *        k d \, ,
 *  @f]
 *  where @f$ \mathbf{A}_{\text{ei}} @f$ is the
 *  @ref EnergyIndependentEigenOperator.
 *
 *  We only implement the second approach, since it represents a smaller
 *  system independent of energy and represents the more canonical
 *  formulation of the problem.  The former approach would be
 *  of value when parallelizing over energy, which is not our focus here.
 *
 *  @section anoteaboutthemultigroupsolve A Note About the Multigroup Solve
 *
 *  Note that in both cases, we require the action
 *  @f$ (\mathbf{I} - \mathbf{DL}^{-1}\mathbf{MS})\f$.  This is equivalent
 *  to an \e exact multigroup solve.
 *
 *  For the power method (see \ref EigenPI), an
 *  exact inversion is not really necessary, since the iteration
 *  will eventually converge (though in many more power iterations
 *  than if an exact multigroup solve were used).  In many
 *  cases, this is actually more efficient, as the number of
 *  sweeps to get within some distance to the true eigenvector
 *  and eigenvalue is less than if one converges on the
 *  multigroup problem (and the within-group problems).
 *
 *  The only
 *  problem with this approach is that the convergence criterion is not
 *  very well defined, since the iteration is not truly a power iteration
 *  when exact inversions are not used (which actually means
 *  we never really use power iteration in practice).
 *  Hence, there is no real
 *  way to map the difference between successive iterates in
 *  the approximate case to what the difference between successive
 *  iterates means mathematically when full inversions are used,
 *  i.e. when a full power iteration is used.  Thus, when one
 *  does \e not converge the multigroup problem, care must be
 *  taken to ensure the convergence criterion at least somewhat
 *  represents the desired outcome.
 *
 *  For Krylov methods such as the Arnoldi method, stricter
 *  convergence on the multigroup problem is required.  This is
 *  because inexact (or simply no) convergence on the multigroup
 *  problem represents an entirely different operator \f$ \mathbf{A} \f$,
 *  and hence will yield an entirely different spectrum.
 *
 *  @section eigensolverinput Eigensolver Parameters
 *
 *  @f[
    \begin{tabular}{llll}
      Key                     & Type  & Description                         & Default \\
    \hline
      eigen\_max\_iters       & int   & Maximum iterations allowed          & 100     \\
      eigen\_tolerance        & dbl   & Tolerance on residual               & 1e-5    \\
      eigen\_print\_out       & int   & Diagnostic print level              & 0       \\
      eigen\_print\_interval  & int   & Iteration interval for diagonostics & 10      \\
    \end{tabular}
    @f]
 *
 */
//----------------------------------------------------------------------------//

template <class D>
class Eigensolver: public Solver<D>
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef FixedSourceManager<D>                     Fixed_T;
  typedef typename Fixed_T::SP_manager              SP_mg_solver;
  typedef detran_utilities::SP<Eigensolver<D> >     SP_solver;
  typedef Solver<D>                                 Base;
  typedef typename Base::SP_input                   SP_input;
  typedef typename Base::SP_state                   SP_state;
  typedef typename Base::SP_mesh                    SP_mesh;
  typedef typename Base::SP_material                SP_material;
  typedef typename Base::SP_boundary                SP_boundary;
  typedef typename Base::SP_fissionsource           SP_fissionsource;


  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param mg_solver  Multigroup solver
   */
  Eigensolver(SP_mg_solver mg_solver);

  /// Virtual destructor
  virtual ~Eigensolver() = 0;

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL EIGENSOLVERS MUST IMPLEMENT
  //--------------------------------------------------------------------------//

  /// Solve the eigenvalue problem.
  virtual void solve() = 0;

protected:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

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

  // Multigroup solver
  SP_mg_solver d_mg_solver;

};

} // namespace detran

#endif /* detran_EIGENSOLVER_HH_ */

//----------------------------------------------------------------------------//
//              end of Eigensolver.hh
//----------------------------------------------------------------------------//
