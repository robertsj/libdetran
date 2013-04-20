//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   detran.cc
 *  @author robertsj
 *  @date   Apr 11, 2012
 *  @brief  detran executable driver.
 */
//---------------------------------------------------------------------------//

#include "detran_config.hh"
#include "StupidParser.hh"
#include "Execute.hh"
#include "solvers/Manager.hh"
#include "utilities/Timer.hh"
#include "utilities/Profiler.hh"
#include <iostream>
#include <ctime>

void print_welcome();

int main(int argc, char **argv)
{

  START_PROFILER();

  detran::Manager::initialize(argc, argv);

  // Print the welcome header.
  print_welcome();

  // Timer
  detran_utilities::Timer timer;

  // Start timer.
  timer.tic();

  // Create parser.
  detran::StupidParser parser(argc, argv);

  // Create problem manager.
  detran::Execute boss(parser);

  // Solve a 1D, 2D, or 3D problem.
  if (boss.dimension() == 1)
    boss.solve<detran::_1D>();
  else if (boss.dimension() == 2)
    boss.solve<detran::_2D>();
  else
    boss.solve<detran::_3D>();

  timer.toc(true);

  detran::Manager::finalize();

  STOP_PROFILER();

  return 0;
}

void print_welcome()
{
  std::time_t t;
  std::time(&t);
  std::cout << std::endl << std::endl;
  std::cout << "       _|              _|                                    " << std::endl;
  std::cout << "   _|_|_|    _|_|    _|_|_|_|  _|  _|_|    _|_|_|  _|_|_|    " << std::endl;
  std::cout << " _|    _|  _|_|_|_|    _|      _|_|      _|    _|  _|    _|  " << std::endl;
  std::cout << " _|    _|  _|          _|      _|        _|    _|  _|    _|  " << std::endl;
  std::cout << "   _|_|_|    _|_|_|      _|_|  _|          _|_|_|  _|    _|  " << std::endl;
  std::cout << " a DETerministic TRANsport tool" << std::endl;
  std::cout << " Built on: " << DETRAN_COMPILED_M << "/"
                             << DETRAN_COMPILED_D << "/"
                             << DETRAN_COMPILED_Y << std::endl;
  std::cout << "   Run on: " << std::ctime(&t);
  std::cout << std::endl << std::endl;
}


