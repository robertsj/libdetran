//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Random.hh
 *  @brief Random class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_utilities_RANDOM_HH_
#define detran_utilities_RANDOM_HH_

namespace detran_utilities
{

/**
 *  @class Random
 *  @brief Linear congruential generator for pseudo-random numbers
 *
 *  Initial implementation uses a fixed set of parameters.  Easy to extend
 *  when needed.
 */
class Random
{

public:

  typedef long        I8;
  typedef double      R8;

  Random(const int S0 = 123)
    : d_g(19073486328125)
    , d_c(0)
    , d_stride(152917)
    , d_S0(S0)
    , d_two_to_m_minus_1(281474976710655)
    , d_one_over_two_to_m(1.0 / 281474976710656.0)
    , d_S(19073486328125)
  {
    /* ... */
  }

  virtual ~Random(){}

  virtual double rnd()
  {
    d_S = (((d_g * d_S) & d_two_to_m_minus_1) + d_c) & d_two_to_m_minus_1;
    return d_S * d_one_over_two_to_m;
  }

protected:

  I8 d_g;               // 5^19
  I8 d_c;               // 0
  I8 d_stride;          // default stride
  I8 d_S0;              // 5^19
  I8 d_two_to_m_minus_1;  // 2^48 - 1
  R8 d_one_over_two_to_m;  // 1/2^48

  /// Current seed, needs to be private to a thread
  I8 d_S;
#pragma omp threadprivate (d_S);

};

} // end namespace detran_utilities

#endif /* detran_utilities_RANDOM_HH_ */

//----------------------------------------------------------------------------//
//              end of file Random.hh
//----------------------------------------------------------------------------//
