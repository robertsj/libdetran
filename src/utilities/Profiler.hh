//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Profiler.hh
 * \author Jeremy Roberts
 * \date   Apr 19, 2012
 * \brief  Wrapper for gperftools profiler.
 * \todo   Fix this so it's not a fixed header location
 */
//---------------------------------------------------------------------------//

#ifndef PROFILER_HH_
#define PROFILER_HH_

#ifdef DETRAN_ENABLE_GPERFTOOLS

#include "/home/robertsj/opt/gperftools/include/google/profiler.h"
#define START_PROFILER() ProfilerStart("detran.prof")
#define STOP_PROFILER()  ProfilerStop()

#else

#define START_PROFILER() ((void) 0)
#define STOP_PROFILER()  ((void) 0)

#endif // DETRAN_ENABLE_GPERFTOOLS

#endif /* PROFILER_HH_ */
