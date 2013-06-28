//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Preconditioner.hh
 *  @brief Preconditioner class definition
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef callow_PRECONDITIONER_HH_
#define callow_PRECONDITIONER_HH_

#include "callow/callow_export.hh"
#include "callow/callow_config.hh"
#include "callow/vector/Vector.hh"
#include <string>

namespace callow
{

/**
 *  @class Preconditioner
 *  @brief Defines a preconditioner for linear solves
 *
 *  Consider the linear system
 *  \f[
 *      \mathbf{A}x = b \, .
 *  \f]
 *  When using iterative methods to solve this system,
 *  one can often make the system easier to solve.  Suppose
 *  we define an operator \f$ \mathbf{P} \f$ such that
 *  \f$ \mathbf{P}^{-1} \mathbf{A} \approx \mathbf{I} \f$.
 *  If the action of \f$ \mathbf{P}^{-1} \f$ is relatively
 *  easy to compute (as compared to inverting \mathbf{A}),
 *  \f$ \mathbf{P} \f$ is a good preconditioner.
 *
 *  We apply a precondition on the left to obtain the
 *  modified system
 *  \f[
 *      \mathbf{P^{-1} A}x = \mathbf{P}^{-1} b \, ,
 *  \f]
 *  or on the right to get
 *  \f[
 *      \mathbf{AP^{-1}} \overbrace{\mathbf{P}x}^{y} =  b \, ,
 *  \f]
 *  following which we solve
 *  \f[
 *      x = \mathbf{P}^{-1} y \, .
 *  \f]
 *
 *  Within callow, the Jacobi and ILU(0) preconditioners are
 *  available along with user-defined shell preconditioners.
 *  If built with PETSc, all preconditioners are available (to PETSc)
 *  as shells.  Otherwise, the user can set PETSc preconditioners
 *  with PetscSolver parameters.  If built with SLEPc, preconditioners
 *  are available for spectral transformations.
 */

class CALLOW_EXPORT Preconditioner
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<Preconditioner>  SP_preconditioner;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /// Constructor
  Preconditioner(std::string name)
    : d_name(name), d_petsc_pc(NULL), d_slepc_st(NULL) {}
  /// Virtual destructor
  virtual ~Preconditioner(){};
  /// set PETSc preconditioner and other setup (called by PetscSolver)
  void set_petsc_pc(PC pc);
  /// set SLEPc spectral transformer and other setup (called by SlepcSolver)
  void set_slepc_st(ST st);
  /// return PETSc preconditioner
  PC petsc_pc() {return d_petsc_pc;}
  /// return SLEPc spectral transformer
  ST slepc_st() {return d_slepc_st;}
  /// return the preconditioner name
  std::string name() const {return d_name;}

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL PRECONDITIONERS MUST IMPLEMENT THIS
  //--------------------------------------------------------------------------//

  /// solve Px = b
  virtual void apply(Vector &b, Vector &x) = 0;

protected:

  /// pc name
  std::string d_name;
  /// PETSc preconditioner
  PC d_petsc_pc;
  /// SLEPc spectral transformation
  ST d_slepc_st;

};

//@{
/// These are the methods actually called by PETSc or SLEPc.  The callow
/// preconditioners are viewed as shells.
PetscErrorCode pc_apply_wrapper(PC pc, Vec b, Vec x);
PetscErrorCode st_apply_wrapper(ST st, Vec b, Vec x);
//@}

CALLOW_TEMPLATE_EXPORT(detran_utilities::SP<Preconditioner>)

} // end namespace callow

#include "Preconditioner.i.hh"

#endif // callow_PRECONDITIONER_HH_

//----------------------------------------------------------------------------//
//              end of file Preconditioner.hh
//----------------------------------------------------------------------------//
