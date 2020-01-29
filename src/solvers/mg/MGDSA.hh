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
 *      (\mathbf{I} + \mathbf{C}^{-1} \mathbf{S}) \, ,
 *  @f]
 *  where \f$ \mathbf{C} \f$ is the multigroup diffusion operator.  This
 *  operator treats group-to-group scattering and, if requested,
 *  fission implicitly.
 *
 *  Because this performs a diffusion solve on the same mesh as
 *  the transport problem, the resulting system can be very
 *  large.  Hence, it is likely to perform best for relatively
 *  small systems.  A coarse mesh version is available (@ref MGCMDSA),
 *  so MGDSA should provide an upper bound for the efficacy of
 *  a diffusion-based preconditioner with respect to spectral sharpening,
 *  though its net computational efficiency may be poor.
 */

class MGDSA: public MGPreconditioner
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef MGPreconditioner                  Base;
  typedef detran_utilities::SP<MGDSA>       SP_pc;
  typedef DiffusionLossOperator             Operator_T;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param input            Input database
   *  @param material         Material database
   *  @param mesh             Cartesian mesh
   *  @param ssource          Scattering source
   *  @param fsource          Fission source
   *  @param cutoff           First group included in solve
   *  @param include_fission  Flag to treat fission like scatter
   *  @param adjoint          Flag for adjoint problems
   */
  MGDSA(SP_input         input,
        SP_material      material,
        SP_mesh          mesh,
        SP_scattersource ssource,
        SP_fissionsource fsource,
        size_t           cutoff,
        bool             include_fission,
        bool             adjoint);

  /// virtual destructor
  virtual ~MGDSA(){}

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL PRECONDITIONERS MUST IMPLEMENT THIS
  //--------------------------------------------------------------------------//

  /// Solve Px = b
  void apply(Vector &b, Vector &x);

  /// build the preconditioner for a new keff
  void build(const double keff = 1.0, SP_state state = SP_state(0));

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

};

} // end namespace detran

#endif // detran_MGDSA_HH_

//----------------------------------------------------------------------------//
//              end of file MGDSA.hh
//----------------------------------------------------------------------------//
