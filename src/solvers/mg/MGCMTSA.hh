//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  MGCMTSA.hh
 *  @brief MGCMTSA class definition.
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#ifndef detran_MGCMTSA_HH_
#define detran_MGCMTSA_HH_

#include "MGPreconditioner.hh"
#include "DiffusionLossOperator.hh"
#include "transport/ScatterSource.hh"
#include "callow/solver/LinearSolver.hh"
#include "callow/preconditioner/Preconditioner.hh"
#include "callow/matrix/MatrixShell.hh"

namespace detran
{

/**
 *  @class MGCMTSA.hh
 *  @brief Multigroup coarse mesh transport synthetic acceleration
 *
 *  Recall the multigroup transport equation can be written generically as
 *  @f[
 *     \mathbf{A} \phi
 *       = (\mathbf{I}-\mathbf{T}_h \mathbf{M}_h \mathbf{S}_h})\phi_h
 *       = \mathbf{T}_hq_h \, ,
 *  @f]
 *  where \f$ h \f$ denotes a fine mesh quantity in space, angle, and energy.
 *  The multigroup TSA preconditioning process \$ \mathbf{\tilde{A}}^{-1} \$
 *  represents an approximation to \f$ \mathbf{A}^{-1} \f$.  Here,
 *  we represent that via the process
 *  @f[
 *   \mathbf{\tilde{A}} =
 *     \mathbf{P}
 *       (\mathbf{I} - \mathbf{T}_H \mathbf{M}_H \mathbf{S}_H)
 *     \mathbf{R} \, ,
 *  @f]
 *  where \f$ H \f$ represents the equivalent operator on a coarse
 *  mesh,
 *  \f$ \mathbf{R} \f$ is the restriction operator, and
 *  \f$ \mathbf{P} \f$ is the prolongation operator.
 *
 *  The restriction space is based on a user-specific level defining the
 *  number of fine cells per coarse cell and/or the number of fine groups
 *  per coarse group.  A coarser angular mesh is created simply by using
 *  a smaller quadrature set.  If only a new angle mesh is desired,
 *  similar to the original TSA methods, the original mesh and materials
 *  are used.
 *
 *  \todo The following should be common functions
 *
 *  The restriction in space is based either on a user-defined state vector,
 *  possibly from an initial guess, partial solution, or other approximation,
 *  or uniform volume weighting.
 *
 *  The restriction in energy is based either on the same user defined state
 *  vector or a set of material-dependent spectra based on infinite medium
 *  calculations for each material.
 */

class MGCMTSA: public MGPreconditioner
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef MGPreconditioner                  Base;
  typedef detran_utilities::SP<MGCMTSA>     SP_pc;
  typedef ScatterSource::SP_scattersource   SP_scattersource;
  typedef DiffusionLossOperator             Operator_T;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *
   *  Assuming the within-group transport problem is set up,
   *  a KSP object exists from which the PC is extracted.  This
   *  PC is passed here to be constructed and for its application
   *  operator to be assigned.
   *
   *  @param input      Input database
   *  @param material   Material database
   *  @param mesh       Cartesian mesh
   *  @param source     Scattering source
   *  @param cutoff     First group included in solve
   */
  MGCMTSA(SP_input input,
        SP_material material,
        SP_mesh mesh,
        SP_scattersource source,
        size_t cutoff,
        bool include_fission);

  /// virtual destructor
  virtual ~MGCMTSA(){}

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL PRECONDITIONERS MUST IMPLEMENT THIS
  //-------------------------------------------------------------------------//

  /// solve Px = b
  void apply(Vector &b, Vector &x);

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL SHELL MATRICES MUST IMPLEMENT THIS
  //-------------------------------------------------------------------------//

  // the client must implement the action y <-- A * x
  void multiply(const Vector &x,  Vector &y)
  {
    Vector b(x.size(), 0.0);
    b.copy(x);
    apply(b, y);
  }

  // the client must implement the action y <-- A' * x
  void multiply_transpose(const Vector &x, Vector &y)
  {
    THROW("NOT IMPLEMENTED");
  }

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Scatter source
  SP_scattersource d_scattersource;

};

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//


#endif /* MGCMTSA_HH_ */
