//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Profiler.hh
 * \author Jeremy Roberts
 * \date   Apr 19, 2012
 * \brief  Wrapper for gperftools profiler.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

#ifndef PROFILER_HH_
#define PROFILER_HH_

#ifdef DETRAN_ENABLE_GPERFTOOLS
#include "/home/robertsj/opt/gperftools/include/google/profiler.h"
#endif


#ifdef DETRAN_ENABLE_GPERFTOOLS

#define START_PROFILER() ProfilerStart("detran.prof")
#define STOP_PROFILER()  ProfilerStop()

#else

#define START_PROFILER() ((void) 0)
#define STOP_PROFILER()  ((void) 0)

#endif


#endif /* PROFILER_HH_ */
