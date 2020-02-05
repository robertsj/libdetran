//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  CallowDefinitions.hh
 *  @brief Various definitions to support callow solvers
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef callow_CALLOWDEFINITIONS_HH_
#define callow_CALLOWDEFINITIONS_HH_

namespace callow
{

//----------------------------------------------------------------------------//
// ENUMERATIONS
//----------------------------------------------------------------------------//

/// return codes for solvers
enum solver_status
{
  RUNNING, SUCCESS, MAXIT, DIVERGE, END_SOLVER_STATUS
};

/**
 *  @brief various vector norms available
 *
 *  For a vector @f$ v \in (m, 1) @f$:
 *    - @f$ ||v||_1         \equiv \sum_i |v_i| @f$
 *    - @f$ ||v||_2         \equiv \sqrt{\sum_i v_i^2} @f$
 *    - @f$ ||v||_{\infty}  \equiv \max_{v_i} |v_i| @f$
 *    - @f$ ||v||_{1g}      \equiv \sum_i |v_i| / m @f$
 *    - @f$ ||v||_{2g}      \equiv \sqrt{ m^{-1} \sum_i v_i^2} @f$
 *
 *  Additional options specific to residual norms include:
 *    - @f$ ||x-y||_{1r}         \equiv \sum_i | (x_i-y_i)/y_i | @f$
 *    - @f$ ||x-y||_{2r}         \equiv \sqrt{\sum_i ((x_i-y_i)/y_i)^2} @f$
 *    - @f$ ||x-y||_{\infty r}   \equiv \max |(x_i-y_i)/y_i| @f$
 */
enum vector_norm_types
{
  L1, L2, LINF, L1GRID, L2GRID, L1REL, L2REL, LINFREL, END_VECTOR_NORM_TYPES
};


} // end namespace callow

#endif // callow_CALLOWDEFINITIONS_HH_

//----------------------------------------------------------------------------//
//              end of file CallowDefinitions.hh
//----------------------------------------------------------------------------//
