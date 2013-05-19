//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  Manager.hh
 *  @brief Manager for external libraries.
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#ifndef detran_MANAGER_HH_
#define detran_MANAGER_HH_

#include "detran_config.hh"
#include "utilities/Profiler.hh"
#include "callow/utils/Initialization.hh"
#include <iostream>

namespace detran
{

/**
 *  @class Manager
 *  @brief Manager for initializing and finalizing external libraries.
 *
 *  This is useful if the user wants finer control on using
 *  detran classes than PyExecute affords but still requires
 *  libraries to be initialized and finalized.
 */
class Manager
{

public:

  /// Initialize libraries (PETSc/SLEPc through callow)
  static void initialize(int argc, char *argv[])
  {
    START_PROFILER();
    callow_initialize(argc, argv);
  }

  /// Finalize libraries
  static void finalize()
  {
    callow_finalize();
    STOP_PROFILER();
  }

};

} // end namespace detran

#endif // detran_MANAGER_HH_

//---------------------------------------------------------------------------//
//              end of file Manager.hh
//---------------------------------------------------------------------------//
