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

namespace detran
{

class Timer
{

public:

  Timer() : d_started(false), d_value(0.0) {}

  void tic()
  {
    d_started = true;
    d_value = (double) std::clock() / (double)CLOCKS_PER_SEC;
  }

  double toc(bool flag = false)
  {
    Require(d_started);
    double etime = wtime() - d_value;
    if (flag) std::cout << " Elapsed time: " << etime << std::endl;
    return etime;
  }

  double wtime()
  {
    return (double) std::clock() / (double)CLOCKS_PER_SEC;
  }

private:

  bool d_started;
  double d_value;

};

}



#endif /* TIMER_HH_ */
