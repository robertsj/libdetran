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
#include "Timer.hh"
#include "Traits.hh"
#include "Profiler.hh"

// System
#include <iostream>
#include <ctime>

int main(int argc, char **argv)
{

  START_PROFILER();

  std::time_t t;
  std::time(&t);
  std::cout << std::endl << std::endl;
  std::cout << "       _|              _|                                    " << std::endl;
  std::cout << "   _|_|_|    _|_|    _|_|_|_|  _|  _|_|    _|_|_|  _|_|_|    " << std::endl;
  std::cout << " _|    _|  _|_|_|_|    _|      _|_|      _|    _|  _|    _|  " << std::endl;
  std::cout << " _|    _|  _|          _|      _|        _|    _|  _|    _|  " << std::endl;
  std::cout << "   _|_|_|    _|_|_|      _|_|  _|          _|_|_|  _|    _|  " << std::endl;
  std::cout << " a DETerministic TRANsport tool" << std::endl;
  std::cout << " Run on: " << std::ctime(&t);
  std::cout << std::endl << std::endl;
  // Timer
  detran::Timer timer;

  // Start timer.
  timer.tic();

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

  timer.toc(true);

  STOP_PROFILER();

  return 0;
}


