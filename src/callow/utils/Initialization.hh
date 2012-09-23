//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Initialization.hh
 * \author robertsj
 * \date   Sep 18, 2012
 * \brief  Initialization class definition.
 */
//---------------------------------------------------------------------------//

#ifndef callow_INITIALIZATION_HH_
#define callow_INITIALIZATION_HH_

#include "callow_config.hh"

#ifdef CALLOW_ENABLE_GPERFTOOLS
#include "/home/robertsj/opt/gperftools/include/google/profiler.h"
#define START_PROFILER() ProfilerStart("callow.prof")
#define STOP_PROFILER()  ProfilerStop()
#else
#define START_PROFILER() ((void) 0)
#define STOP_PROFILER()  ((void) 0)
#endif // CALLOW_ENABLE_GPERFTOOLS

/// Initialize external packages, if enabled
inline void callow_initialize(int argc, char** argv)
{
#ifdef CALLOW_ENABLE_PETSC
  PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
#endif
#ifdef CALLOW_ENABLE_SLEPC
  SlepcInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
#endif
  START_PROFILER();
}

/// Finalize external packages, if enabled
inline void callow_finalize()
{
  STOP_PROFILER();
#ifdef CALLOW_ENABLE_PETSC
  PetscFinalize();
#endif
#ifdef CALLOW_ENABLE_SLEPC
  SlepcFinalize();
#endif
}

#endif /* callow_INITIALIZATION_HH_ */
