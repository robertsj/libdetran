//-----------------------------------*-C-*-----------------------------------//
/**
 *  @file   callow_config.hh
 *  @author Jeremy Roberts
 *  @brief  Definitions for callow packages.
 */
//---------------------------------------------------------------------------//

/**
 *  @page callow callow: a naive linear algebra package
 *
 *  Originally, detran was built with an increasingly ugly
 *  layer of PETSc solvers for using Krylov solvers on various
 *  problems.  Now, the linear algebra, including the solvers
 *  and matrix operators, has been mostly abstracted away into
 *  the callow namespace.  In addition, a few
 *  solvers have been implemented (including GMRES(m)) that
 *  lets users have solver options without building with PETSc
 *  (though PETSc provides many more options along with much
 *  better preconditioners).
 *
 *  The callow library can be used as a stand alone utility,
 *  depending on only the utility functions of detran.  Hence,
 *  it offers a simple sparse linear algebra library for C++
 *  applications.
 *
 */

#ifndef CALLOW_CONFIG_HH_
#define CALLOW_CONFIG_HH_

// disable signed-to-unsigned
#pragma warning (disable : 4018)
// disable multiple copy ctor
#pragma warning (disable : 4521)
// safe fopen, etc.
#pragma warning (disable : 4996)
// int to bool
#pragma warning (disable : 4800)

// include detran configuration
#include "detran_config.hh"

// callow export stuff
#include "callow/callow_export.hh"

#ifdef DETRAN_ENABLE_PETSC
#define CALLOW_ENABLE_PETSC
#include "petsc.h"
#else
namespace callow
{
typedef int* KSP;
typedef int* PC;
typedef int* Mat;
typedef int* Vec;
typedef int PetscErrorCode;
}
#endif

#ifdef DETRAN_ENABLE_SLEPC
#define CALLOW_ENABLE_SLEPC
#include "slepc.h"
#else
namespace callow
{
typedef int* EPS;
typedef int* ST;
}
#endif

//#define CALLOW_ENABLE_PETSC_OPS

/// Provides linear algebra support for detran
namespace callow
{
  /* ... */
}

#endif // CALLOW_CONFIG_HH_ 

//---------------------------------------------------------------------------//
//              end of callow_config.hh
//---------------------------------------------------------------------------//
