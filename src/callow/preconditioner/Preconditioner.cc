//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Preconditioner.cc
 *  @brief Preconditioner member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "Preconditioner.hh"

namespace callow
{

//----------------------------------------------------------------------------//
Preconditioner::Preconditioner(const std::string &name)
  : d_name(name)
  , d_petsc_pc(NULL)
  , d_slepc_st(NULL)
  , d_size(0)
{

}

//----------------------------------------------------------------------------//
void Preconditioner::display(const std::string &name)
{
  Assert(d_size > 0);
  MatrixShellPC M(this, d_size);
  M.compute_explicit(name);
}

} // end namespace callow

//----------------------------------------------------------------------------//
//              end of file Preconditioner.hh
//----------------------------------------------------------------------------//
