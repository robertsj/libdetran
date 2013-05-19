//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  MGDSA.hh
 *  @brief MGDSA class definition
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_MGDSA_HH_
#define detran_MGDSA_HH_

#include "MGPreconditioner.hh"
#include "DiffusionLossOperator.hh"
#include "transport/ScatterSource.hh"
#include "callow/solver/LinearSolver.hh"
#include "callow/preconditioner/Preconditioner.hh"
#include "callow/matrix/MatrixShell.hh"

namespace detran
{

/**
 *  @class MGDSA
 *  @brief Multigroup diffusion synthetic acceleration
 *
 *  The multigroup DSA preconditioning process \$ \mathbf{P}^{-1} \$
 *  is defined to be
 *  @f[
 *      (\mathbf{I} - \mathbf{C}^{-1} \mathbf{S}) \, ,
 *  @f]
 *  where \f$ \mathbf{C} \f$ is the multigroup diffusion operator.  This
 *  operator treats group-to-group scattering and, if requested,
 *  fission implicitly.
 *
 *  Because this performs a diffusion solve on the same mesh as
 *  the transport problem, the resulting system can be very
 *  large.  Hence, it is likely to perform best for relatively
 *  small systems.  A coarse mesh version is under development, but
 *  this implementation will serve as the upper bound for the
 *  efficacy of diffusion-based multigroup preconditioning.
 *
 *  @note This inherits from the shell matrix (for now) so that the
 *        action can be used to construct an explicit operator for
 *        detailed numerical studies
 */

class MGDSA: public MGPreconditioner
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef MGPreconditioner                  Base;
  typedef detran_utilities::SP<MGDSA>       SP_pc;
  typedef ScatterSource::SP_scattersource   SP_scattersource;
  typedef DiffusionLossOperator             Operator_T;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

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
  MGDSA(SP_input input,
        SP_material material,
        SP_mesh mesh,
        SP_scattersource source,
        size_t cutoff,
        bool include_fission);

  /// virtual destructor
  virtual ~MGDSA(){}

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL PRECONDITIONERS MUST IMPLEMENT THIS
  //--------------------------------------------------------------------------//

  /// solve Px = b
  void apply(Vector &b, Vector &x);

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL SHELL MATRICES MUST IMPLEMENT THIS
  //--------------------------------------------------------------------------//

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

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Scatter source
  SP_scattersource d_scattersource;

};

} // end namespace detran

#endif // detran_MGDSA_HH_

//----------------------------------------------------------------------------//
//              end of file MGDSA.hh
//----------------------------------------------------------------------------//
