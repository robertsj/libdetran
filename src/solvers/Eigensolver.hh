//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Eigensolver.hh
 * \author robertsj
 * \date   Jun 18, 2012
 * \brief  Eigensolver class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

#ifndef EIGENSOLVER_HH_
#define EIGENSOLVER_HH_

// Detran
#include "ExternalSource.hh"
#include "FissionSource.hh"
#include "GaussSeidel.hh"
#include "MultigroupSolver.hh"
#include "Material.hh"
#include "Mesh.hh"
#include "InnerIteration.hh"
#include "detran_config.h"
#ifdef DETRAN_ENABLE_PETSC
#include "KrylovMG.hh"
#endif

// Utilities
#include "DBC.hh"
#include "InputDB.hh"
#include "SP.hh"

// System
#ifdef DETRAN_ENABLE_PETSC
#include "petsc.h"
#endif

namespace detran
{

//---------------------------------------------------------------------------//
/*!
 *  \class Eigensolver
 *  \brief Base class for solving the eigenvalue problem
 *
 *  \section theeigenvalueproblem The Eigenvalue Problem
 *
 *  The steady-state balance of neutrons in a fissile system is
 *  characterized by the multigroup eigenvalue problem
 *
 *  \f[
 *      \hat{\Omega} \cdot \nabla  \psi +
 *        \Sigma_t(\vec{r},E_g) \psi(\vec{r},\hat{\Omega},E_g) =
 *      \frac{1}{4\pi} \sum^G_{g'} \int_{4\pi} d\Omega'
 *        \Sigma_S(\vec{r},\hat{\Omega}'\cdot \hat{\Omega}, E_g'\to E_g)
 *          \psi(\vec{r},\hat{\Omega}',E_{g'}) +
 *      \frac{\chi_g}{4\pi k} \sum^G_{g'} \int_{4\pi} d\Omega'
 *        \nu \Sigma_f(\vec{r}, E_{g'}) \psi(\vec{r}, \hat{\Omega}', E_{g'})
 *  \f]
 *
 *  where \f$ k \f$ is the eigenvalue, which physically represents the
 *  ratio of the number of neutrons in successive generations.
 *
 *  \section operatorform Operator Form
 *
 *  In operator form, the eigenvalue problem is
 *
 *  \f[
 *      (\mathbf{I} - \mathbf{DL}^{-1}\mathbf{MS})\phi =
 *        \frac{1}{k} \mathbf{DL}^{-1}\mathbf{M}\chi \mathbf{F}^T \phi \, .
 *  \f]
 *
 *  where the eigenvector consists of the angular flux moments.  In the more
 *  standard form \f$ \mathbf{A}x = k x \f$, we have two options.  We
 *  can keep the problem in terms of the multigroup fluxes so that
 *
 *  \f[
 *      \overbrace{(\mathbf{I} - \mathbf{DL}^{-1}\mathbf{MS})^{-1}
 *        \mathbf{DL}^{-1}\mathbf{M}\chi \mathbf{F}^T}^{\mathbf{A}} \phi =
 *        k \phi \, .
 *  \f]
 *
 *  Alternatively, we can solve in terms of the fission density,
 *  \f[
 *      d = \mathbf{F}^T \phi = \sum^G_g \nu \Sigma_{fg} \phi_g
 *  \f]
 *  so that
 *  \f[
 *      \overbrace{(\mathbf{I} - \mathbf{DL}^{-1}\mathbf{MS})^{-1}
 *        \mathbf{DL}^{-1}\mathbf{M}\chi }^{\mathbf{A}} d =
 *        k d \, .
 *  \f]
 *
 *  We choose the second approach, since it represents a smaller
 *  system independent of energy and represents the more canonical
 *  formulation of the problem.  The former approach would be
 *  of value when parallelizing over energy, as in Denovo, but
 *  that is not our focus here.
 *
 *  \section anoteaboutthemultigroupsolve A Note About the Multigroup Solve
 *
 *  Note that in both cases, we require the action
 *  \f$ (\mathbf{I} - \mathbf{DL}^{-1} )\f$.  This is equivalent
 *  to an \e exact multigroup solve.
 *
 *  For the power method (see \ref PowerIteration), an
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
 *  \section eigensolverinput Eigensolver Input Entries
 *
 *
 *  \f[
    \begin{tabular}{llll}
      Key                     & Type  & Description                         & Default \\
    \hline
      eigen\_max\_iters       & int   & Maximum iterations allowed          & 100     \\
      eigen\_tolerance        & dbl   & Tolerance on residual               & 1e-5    \\
      eigen\_print\_out       & int   & Diagnostic print level              & 0       \\
      eigen\_print\_interval  & int   & Iteration interval for diagonostics & 10      \\
    \end{tabular}
    \f]
 *
 */
//---------------------------------------------------------------------------//

template <class D>
class Eigensolver: public Object
{

public:

  typedef SP<Eigensolver<D> >                   SP_solver;
  typedef typename
      MultigroupSolver<D>::SP_solver            SP_mg_solver;
  // basic objects
  typedef InputDB::SP_input                     SP_input;
  typedef State::SP_state                       SP_state;
  typedef Mesh::SP_mesh                         SP_mesh;
  typedef Material::SP_material                 SP_material;
  typedef Quadrature::SP_quadrature             SP_quadrature;
  typedef typename BoundaryBase<D>::SP_boundary SP_boundary;
  // source typedefs
  typedef ExternalSource::SP_source             SP_externalsource;
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
   *  \param external_source   User-defined external source.
   *  \param fission_source    Fission source.
   */
  Eigensolver(SP_input           input,
              SP_state           state,
              SP_mesh            mesh,
              SP_material        material,
              SP_quadrature      quadrature,
              SP_boundary        boundary,
              SP_fissionsource   q_f);

  /// Virtual destructor
  virtual ~Eigensolver(){}

  /// Solve the eigenvalue problem.
  virtual void solve() = 0;

  /// DBC function.
  virtual bool is_valid() const
  {
    return true;
  }

protected:

  /// \name Protected Data
  /// \{

  /// User input.
  SP_input b_input;

  /// State vectors.
  SP_state b_state;

  /// Problem mesh.
  SP_mesh b_mesh;

  /// Materials definitions.
  SP_material b_material;

  /// Angular mesh.
  SP_quadrature b_quadrature;

  /// Boundary fluxes.
  SP_boundary b_boundary;

  /// Fission source.
  SP_fissionsource b_fissionsource;

  /// multigroup solver
  SP_mg_solver b_mg_solver;

  /// Maximum eigensolver iterations
  int b_max_iters;

  /// Eigensolver tolerance, \f$ ||\mathbf{A} x-\lambda x|| < tol \f$.
  double b_tolerance;

  /// Print out flag
  int b_print_out;

  /// Interval for print out
  int b_print_interval;

  /// \}

};

// Constructor
template <class D>
Eigensolver<D>::Eigensolver(SP_input          input,
                            SP_state          state,
                            SP_mesh           mesh,
                            SP_material       material,
                            SP_quadrature     quadrature,
                            SP_boundary       boundary,
                            SP_fissionsource  q_f)
  : b_input(input)
  , b_state(state)
  , b_mesh(mesh)
  , b_material(material)
  , b_quadrature(quadrature)
  , b_boundary(boundary)
  , b_fissionsource(q_f)
  , b_max_iters(100)
  , b_tolerance(1e-5)
  , b_print_out(2)
  , b_print_interval(10)
{
  Require(b_input);
  Require(b_state);
  Require(b_mesh);
  Require(b_material);
  Require(b_quadrature);
  Require(b_boundary);
  Require(b_fissionsource);

  // Get relevant input parameters.
  if (input->check("eigen_max_iters"))
    b_max_iters = input->get<int>("eigen_max_iters");

  if (input->check("eigen_tolerance"))
    b_tolerance = input->get<double>("eigen_tolerance");

  if (input->check("eigen_print_out"))
    b_print_out = input->get<int>("eigen_print_out");

  if (input->check("eigen_print_interval"))
    b_print_interval = input->get<int>("eigen_print_interval");

  // Get the multigroup solver type and create.
  std::string outer_solver = "GS";
  if (input->check("outer_solver"))
  {
    outer_solver = input->get<std::string>("outer_solver");
  }
  if (outer_solver == "GS")
  {
    b_mg_solver = new GaussSeidel<D>(input, state, mesh, material,
                                     quadrature, boundary,
                                     SP_externalsource(), q_f);
  }
  else if (outer_solver == "KrylovMG")
  {
#ifdef DETRAN_ENABLE_PETSC
    b_mg_solver = new KrylovMG<D>(input, state, mesh, material,
                                  quadrature, boundary,
                                  SP_externalsource(), q_f);
#else
    THROW("KrylovMG is not available because PETSc is not enabled.");
#endif
  }
  else
  {
    THROW("Unsupported outer solver type selected: "+outer_solver);
  }

}

} // namespace detran

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

//#include "Eigensolver.i.hh"

#endif /* EIGENSOLVER_HH_ */
