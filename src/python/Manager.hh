//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Manager.hh
 * \brief  Manager for external libraries.
 * \author Jeremy Roberts
 * \date   Jul 26, 2012
 */
//---------------------------------------------------------------------------//

#ifndef Manager_HH_
#define Manager_HH_

// Configuration
#include "detran_config.h"

#ifdef DETRAN_ENABLE_PETSC
#include "petsc.h"
#endif

#ifdef DETRAN_ENABLE_SLEPC
#include "slepc.h"
#endif

namespace detran
{

/*!
 *  \class Manager
 *  \brief Manager for initializing and finalizing external libraries.
 *
 *  This is useful if the user wants finer control on using
 *  detran classes than PyExecute affords but still requires
 *  libraries to be initialized and finalized.
 */
class Manager
{

public:

  /// Manager PETSc and SLEPc, if present.
  static void initialize(int argc, char *argv[])
  {
    #ifdef DETRAN_ENABLE_PETSC
    // Start PETSc.
    PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
    #endif

    #ifdef DETRAN_ENABLE_SLEPC
    // Start SLEPc.
    SlepcInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
    #endif
  }

  /// Manager PETSc and SLEPc, if present.
  static void finalize()
  {
    #ifdef DETRAN_ENABLE_SLEPC
    // Finish SLEPc.
    SlepcFinalize();
    #endif

    #ifdef DETRAN_ENABLE_PETSC
    // Finish PETSc.
    PetscFinalize();
    #endif
  }

};

} // end namespace detran

#endif // MANAGER_HH_

//---------------------------------------------------------------------------//
//              end of file Manager.hh
//---------------------------------------------------------------------------//
