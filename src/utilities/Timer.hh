/*
 * Timer.hh
 *
 *  Created on: Apr 11, 2012
 *      Author: robertsj
 */

#ifndef TIMER_HH_
#define TIMER_HH_

#include "DBC.hh"

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <map>
#include <string>

namespace detran
{

class Timer
{

public:

  Timer() : d_started(false), d_value(0.0) {}

  /// Begin the timer around a single code block
  void tic()
  {
    d_started = true;
    d_value = (double) std::clock() / (double)CLOCKS_PER_SEC;
  }

  /// Return the time elapsed around a single code block
  double toc(bool flag = false)
  {
    Require(d_started);
    double etime = wtime() - d_value;
    if (flag) std::cout << " Elapsed time: " << etime << std::endl;
    return etime;
  }

  /// Return the current wall time (only differences are meaningful)
  double wtime()
  {
    return (double) std::clock() / (double)CLOCKS_PER_SEC;
  }

  /// Begin the timer at the beginning of a function call.
  void function_tic()
  {
    d_started = true;
    d_value = (double) std::clock() / (double)CLOCKS_PER_SEC;
  }

  /// Log the function time.
  void function_toc(std::string function)
  {
    Require(d_started);
    double etime = wtime() - d_value;
    d_function_map[function] += etime;
  }

  /// Display function elapsed times.
//  void display()
//  {
//    for (int i = 0; i < d_function_map.size(); i++)
//    {
//      std::cout << " Function " << (d_function_map.begin()+i)->first
//                << " Elapsed time: " << (d_function_map.begin()+i)->second
//                << std::endl;
//    }
//  }

private:

  bool d_started;
  double d_value;
  std::map<std::string, double> d_function_map;

};

}



#endif /* TIMER_HH_ */
