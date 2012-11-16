//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   BDFCoefficients.hh
 *  @brief  BDFCoefficients
 *  @author Jeremy Roberts
 *  @date   Nov 15, 2012
 */
//---------------------------------------------------------------------------//

#ifndef detran_BDFCOEFFICIENTS_HH_
#define detran_BDFCOEFFICIENTS_HH_

namespace detran
{

/**
 * The BDF formulas approximate
 *   a(n+1)*y(n+1) - a(n)*y(n) - a(n-1)*y(n) ... = dt*f(y(n+1))
 * This table gives those parameters
 */
const double bdf_coefs[6][7] =
{ // n+1          n-0         n-1         n-2          n-3        n-4        n-5
  {  1.0     ,   1.0     ,    0.0     ,   0.0     ,    0.0     ,  0.0     ,   0.0     }, // order 1
  {  3.0/ 2.0,   4.0/ 2.0,   -1.0/ 2.0,   0.0     ,    0.0     ,  0.0     ,   0.0     }, // order 2
  { 11.0/ 6.0,  18.0/ 6.0,   -9.0/ 6.0,   2.0/ 6.0,    0.0     ,  0.0     ,   0.0     }, // order 3
  { 25.0/12.0,  48.0/12.0,  -36.0/12.0,  16.0/12.0,   -3.0/12.0,  0.0     ,   0.0     }, // order 4
  {137.0/60.0, 300.0/60.0, -300.0/60.0, 200.0/60.0,  -75.0/60.0, 12.0/60.0,   0.0     }, // order 5
  {147.0/60.0, 360.0/60.0, -450.0/60.0, 400.0/60.0, -225.0/60.0, 72.0/60.0, -10.0/60.0}, // order 6
};

} // end namespace detran

#endif // detran_BDFCOEFFICIENTS_HH_

//---------------------------------------------------------------------------//
//              end of file BDFCoefficients.hh
//---------------------------------------------------------------------------//
