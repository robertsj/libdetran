//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  detran.cc
 *  @brief Detran executable driver
 *  @note  Copyright (C) Jeremy Roberts 2012-2013
 */
//----------------------------------------------------------------------------//

#include "detran_config.hh"
#include "StupidParser.hh"
#include "Execute.hh"
#include "solvers/Manager.hh"
#include "utilities/Timer.hh"
#include "utilities/Profiler.hh"
#include "utilities/Definitions.hh"
#include <iostream>
#include <ctime>
#include <cstdio>

void print_welcome();

//----------------------------------------------------------------------------//
int main(int argc, char **argv)
{
 
  START_PROFILER();

  detran::Manager::initialize(argc, argv);

  print_welcome();

  detran_utilities::Timer timer;
  timer.tic();

  try
  {
	  detran::StupidParser parser(argc, argv);
    detran::Execute boss(parser);
    if (boss.dimension() == 1)
      boss.solve<detran::_1D>();
    else if (boss.dimension() == 2)
      boss.solve<detran::_2D>();
    else
      boss.solve<detran::_3D>();
  }
  catch (detran_utilities::GenException &e)
  {
	std::cout << e.what() << std::endl;
	return 0;
  }
  catch (...)
  {
    std::cout << "Unknown exception!" << std::endl;
	return 0;
  }

  timer.toc(true);

  detran::Manager::finalize();

  STOP_PROFILER();

  return 0;
}

//----------------------------------------------------------------------------//
void print_welcome()
{
  using std::cout;
  using std::endl;
  std::time_t t;
  std::time(&t);
  cout << endl << endl;
  cout << "       _|              _|                                  " << endl;
  cout << "   _|_|_|    _|_|    _|_|_|_|  _|  _|_|    _|_|_|  _|_|_|  " << endl;
  cout << " _|    _|  _|_|_|_|    _|      _|_|      _|    _|  _|    _|" << endl;
  cout << " _|    _|  _|          _|      _|        _|    _|  _|    _|" << endl;
  cout << "   _|_|_|    _|_|_|      _|_|  _|          _|_|_|  _|    _|" << endl;
  cout << " a DETerministic TRANsport tool"    << endl;
  cout << " Copyright (C) Jeremy Roberts 2012-2013" << endl;
  cout << " Built on: " << DETRAN_COMPILED_M << "/"
                        << DETRAN_COMPILED_D << "/"
                        << DETRAN_COMPILED_Y << endl;
  cout << " Git SHA1: " << DETRAN_GIT_SHA1   << endl;
  cout << "   Run on: " << std::ctime(&t) << endl << endl;
}


