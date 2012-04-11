//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   detran.cc
 * \author robertsj
 * \date   Apr 11, 2012
 * \brief  detran executable driver.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// Detran
#include "StupidParser.hh"
#include "Execute.hh"
#include "Traits.hh"

int main(int argc, char **argv)
{

  // Create parser.
  detran::StupidParser parser(argc, argv);

  // Create problem manager.
  detran::Execute boss(parser);

  // Solve a 1D, 2D. or 3D problem.
  if (boss.dimension() == 1)
    boss.solve<detran::_1D>();
  else if (boss.dimension() == 2)
    boss.solve<detran::_2D>();
  else
    boss.solve<detran::_3D>();

  return 0;
}


