//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   DiffusionFixedSourceSolver.hh
 * \brief  DiffusionFixedSourceSolver class definition
 * \author Jeremy Roberts
 * \date   Sep 11, 2012
 */
//---------------------------------------------------------------------------//

#ifndef DIFFUSIONFIXEDSOURCESOLVER_HH_
#define DIFFUSIONFIXEDSOURCESOLVER_HH_

#include "DiffusionLossOperator.hh"
#include "DiffusionGainOperator.hh"
#include "Vector.hh"
#include "boundary/BoundaryDiffusion.hh"
#include "external_source/ExternalSource.hh"
#include "transport/State.hh"
#include "utilities/Definitions.hh"
#include "petsc.h"

namespace detran
{

/*!
 *  \class DiffusionFixedSourceSolver
 *  \brief Solve a fixed source problem using diffusion
 *
 *  Fixed source problems can be solved in the following modalities:
 *    - fixed, no multiplication
 *    - fixed, multiplying via implicit fission (i.e. in the loss matrix)
 *    - fixed, multiplying via fission iteration
 *
 *  Mathematically, consider the diffusion equation in
 *  operator form:
 *
 *  \f[
 *      T \phi = \frac{1}{k} F \phi + Q \, .
 *  \f]
 *  The first case assumes \f$ F = 0 \f$.  The second case
 *  uses the modified operator \f$ T' = T -  \frac{1}{k} F \phi \f$.
 *  The third case performs the iteration
 *  \f[
 *      T \phi^{n+1} = \frac{1}{k}F\phi^{n} + Q \, .
 *  \f]
 *  This third case is less efficient than the second case, but it
 *  allows use to pull out the solution following each fission
 *  iteration.  The number of such fission iterations can be limited
 *  by the user.
 *
 *  These three cases are selected via diffusion_fixed_type 0,1,2
 */
template <class D>
class DiffusionFixedSourceSolver
{

public:

  //-------------------------------------------------------------------------//
  // ENUMERATIONS
  //-------------------------------------------------------------------------//

  enum fixed_types
  {
    FIXED, MULTIPLY, ITERATE, END_FIXED_TYPES
  };

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::InputDB::SP_input         SP_input;
  typedef State::SP_state                             SP_state;
  typedef detran_geometry::Mesh::SP_mesh              SP_mesh;
  typedef detran_material::Material::SP_material      SP_material;
  typedef BoundaryDiffusion<D>                        Boundary_T;
  typedef typename Boundary_T::SP_boundary            SP_boundary;
  typedef detran_external_source::
          ExternalSource::SP_externalsource           SP_source;
  typedef State::moments_type                         moments_type;
  typedef DiffusionLossOperator::SP_lossoperator      SP_lossoperator;
  typedef DiffusionGainOperator::SP_gainoperator      SP_gainoperator;
  typedef Vector::SP_vector                           SP_vector;
  typedef detran_utilities::size_t                    size_t;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  DiffusionFixedSourceSolver(SP_input input,
                             SP_material material,
                             SP_mesh mesh,
                             SP_state state);

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Build the right hand side
   *
   *  The source consists of external volume and boundary sources.
   *  Currently, the client must set the boundary current directly for
   *  a boundary source.
   *
   *  Note, the volume source is optional.  If it is not present, only
   *  the boundary current contributes.
   *
   *  \param q  Pointer to volume source
   */
  void build_source(SP_source q = SP_source(0));

  /// Solve the fixed source diffusion problem
  void solve();

  // Getters
  SP_input input() const { return d_input; }
  SP_material material() const { return d_material; }
  SP_mesh mesh() const { return d_mesh; }
  SP_state state() const { return d_state; }
  SP_boundary boundary() const { return d_boundary; }

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
  /// Saved unknown vector
  SP_vector d_phi_old;
  /// Fixed source
  SP_vector d_Q;
  /// Fixed source plus fission
  SP_vector d_Q_total;
  size_t d_problem_size;
  size_t d_fixed_type;
  size_t d_maximum_iterations;
  size_t d_maximum_fission_iterations;
  double d_tolerance;
  double d_fission_tolerance;
  double d_fission_scaling;
  KSP d_solver;
  int d_print_out;
  int d_print_interval;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  // Solve directly via M*phi = Q
  void solve_fixed();
  // Solve with iteration on the fission source, M*phi(n+1) = (F/k)*phi(n) + Q
  void solve_iterate();

  // Build the volume source
  void build_volume_source(SP_source q);
  // Build the boundary source
  void build_boundary_source();

  // Fill the boundary with outgoing current
  void fill_boundary();

};

} // end namespace detran

#endif // DIFFUSIONFIXEDSOURCESOLVER_HH_ 

//---------------------------------------------------------------------------//
//              end of file DiffusionFixedSourceSolver.hh
//---------------------------------------------------------------------------//
