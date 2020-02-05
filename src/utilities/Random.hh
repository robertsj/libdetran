//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Random.hh
 *  @brief Random class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_utilities_RANDOM_HH_
#define detran_utilities_RANDOM_HH_

#include <cmath>

namespace detran_utilities
{

/**
 *  @class Random
 *  @brief Linear congruential generator for pseudo-random numbers
 *
 *  LCG's are one of the most common generators in practice and are
 *  defined by the sequence
 *
 *  @f[
 *      S_{i+1} = \Big ( S_{i} \cdot g + c \Big ) \: \textrm{mod} \: p
 *  @f]
 *
 *  where
 *  - @f$ S_0 @f$ is the <em>seed</em>
 *  - @f$ g   @f$ is <em>generator</em> (or multiplier)
 *  - @f$ c   @f$ is <em>adder</em> (or increment)
 *  - @f$ p   @f$ is <em>modulus</em>
 *
 *  Given the seed @f$ S_0 @f$, any subsequent value @f$ S_k @f$ can
 *  be defined from which we obtain the random number @f$ r_k = S_k / p @f$.
 *
 *  In this implementation, we limit the modulus to be a power of 2, i.e.
 *  @f$ p = 2^m @f$ for a positive integer @f$ m @f$.  Doing so allows us
 *  to employ fast bit operations.
 *
 *  The default parameters selected are those used (at least at one
 *  point) in the MCNP code.
 *
 */
class Random
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef long        I8;
  typedef double      R8;

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param   S    Seed
   *  @param   g    generator
   *  @param   c    adder
   *  @param   m    log-2 of the modulus
   */
  Random(const I8  S = 123,
         const I8  g = 19073486328125,
         const I8  c = 0,
         const I8  m = 48)
    : d_S(S)
    , d_g(19073486328125)
    , d_c(0)
    , d_two_to_m_minus_1(std::pow(2, m)-1)
    , d_one_over_two_to_m(1.0/std::pow(2.0, m))
  {
    /* ... */
  }

  virtual ~Random(){}

  /// Return a value from a uniform distribution over [0, 1]
  virtual double rnd()
  {
    d_S = (((d_g * d_S) & d_two_to_m_minus_1) + d_c) & d_two_to_m_minus_1;
    return d_S * d_one_over_two_to_m;
  }

protected:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Seed
  I8 d_S;
//#pragma omp threadprivate (d_S);
  /// Generator (or multiplier)
  I8 d_g;
  /// Adder (or increment)
  I8 d_c;
  /// Modulus minus one
  I8 d_two_to_m_minus_1;
  /// One over the modulus
  R8 d_one_over_two_to_m;

};

} // end namespace detran_utilities

#endif /* detran_utilities_RANDOM_HH_ */

//----------------------------------------------------------------------------//
//              end of file Random.hh
//----------------------------------------------------------------------------//
