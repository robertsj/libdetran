//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Preconditioner.i.hh
 *  @brief Preconditioner inline member definitions
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef callow_PRECONDITIONER_I_HH_
#define callow_PRECONDITIONER_I_HH_

namespace callow
{

#ifdef DETRAN_ENABLE_PETSC

//----------------------------------------------------------------------------//
inline void Preconditioner::set_petsc_pc(PC pc)
{
  // save the pc object
  d_petsc_pc = pc;

  PetscErrorCode ierr;
  // set the shell preconditioner
  ierr = PCSetType(pc, PCSHELL);
  // set the PC context
  ierr = PCShellSetContext(d_petsc_pc, this);
  Insist(!ierr, "Error setting shell preconditioner context.");
  // set the PC operator
  ierr = PCShellSetApply(d_petsc_pc, pc_apply_wrapper);
  Insist(!ierr, "Error setting shell preconditioner operator.");
  // set the PC name for good measure
  ierr = PCShellSetName(d_petsc_pc, d_name.c_str());
  Insist(!ierr, "Error within-group preconditioner name.");

  Ensure(!ierr);
}

//----------------------------------------------------------------------------//
inline PetscErrorCode pc_apply_wrapper(PC pc, Vec b, Vec x)
{
  // get the context and cast
  PetscErrorCode ierr;
  void *context;
  ierr = PCShellGetContext(pc, &context); CHKERRQ(ierr);
  Preconditioner *foo = (Preconditioner *) context;
  // wrap the petsc vectors
  Vector B(b);
  Vector X(x);
  // call the actual apply operator.
  foo->apply(B, X);
  return ierr;
}

#else

inline void Preconditioner::set_petsc_pc(PC pc){}
inline PetscErrorCode pc_apply_wrapper(PC pc, Vec b, Vec x){return 0;}

#endif

#ifdef DETRAN_ENABLE_SLEPC

//----------------------------------------------------------------------------//
inline void Preconditioner::set_slepc_st(ST st)
{
  std::cout << "SETTING ST WRAPPER!!!" << std::endl;

  // save the pc object
  d_slepc_st = st;

  PetscErrorCode ierr;
  // set the shell preconditioner
  //ierr = STSetType(st, STSHELL);
  // set the PC context
  ierr = STShellSetContext(d_slepc_st, this);
  Insist(!ierr, "Error setting shell spectral transform context.");
  // set the PC operator
  ierr = STShellSetApply(d_slepc_st, st_apply_wrapper);
  Insist(!ierr, "Error setting shell spectral transform operator.");

  Ensure(!ierr);
}

//----------------------------------------------------------------------------//
inline PetscErrorCode st_apply_wrapper(ST st, Vec b, Vec x)
{
  std::cout << "APPLYING ST WRAPPER!!!" << std::endl;
  // get the context and cast
  PetscErrorCode ierr;
  void *context;
  ierr = STShellGetContext(st, &context); CHKERRQ(ierr);
  Preconditioner *foo = (Preconditioner *) context;
  // wrap the petsc vectors
  Vector B(b);
  Vector X(x);
  // call the actual apply operator.
  foo->apply(B, X);
  return ierr;
}

#else

inline void Preconditioner::set_slepc_st(ST st){}
inline PetscErrorCode st_apply_wrapper(ST st, Vec b, Vec x){return 0;}

#endif


} // end namespace detran

#endif // callow_PRECONDITIONER_I_HH_

//----------------------------------------------------------------------------//
//              end of file Preconditioner.i.hh
//----------------------------------------------------------------------------//
