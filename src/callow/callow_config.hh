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

// include detran configuration
#include "config/detran_config.hh"

#ifdef DETRAN_ENABLE_PETSC
#define CALLOW_ENABLE_PETSC
#include "petsc.h"
#endif

#ifdef DETRAN_ENABLE_SLEPC
#define CALLOW_ENABLE_SLEPC
#include "slepc.h"
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
