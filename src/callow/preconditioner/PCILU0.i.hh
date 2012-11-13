//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   PCILU0.i.hh
 *  @brief  PCILU0 inline member definitions
 *  @author Jeremy Roberts
 *  @date   Sep 18, 2012
 */
//---------------------------------------------------------------------------//

#ifndef callow_PCILU0_I_HH_
#define callow_PCILU0_I_HH_

#include <iostream>

namespace callow
{

//---------------------------------------------------------------------------/
inline void PCILU0::apply(Vector &b, Vector &x)
{
  // solve LUx = x --> x = inv(U)*inv(L)*x
  THROW("DONE");
  // forward substitution
  //   for i = 0:m-1
  //     x[i] = 1/L[i,i] * ( b[i] - sum(k=0:i-1, L[i,k]*y[k]) )
  // but note that in our ILU(0) scheme, L is *unit* lower triangle,
  // meaning L has ones on the diagonal (whereas U does not)
  for (int i = 0; i < d_P->number_rows(); ++i)
  {
    // start index
    int s = d_P->start(i);
    // diagonal index
    int d = d_P->diagonal(i);
    // invert row
    d_y[i] = b[i];
    for (int p = s; p < d; ++p)
    {
      // column index
      int c = d_P->column(p);
      d_y[i] -= d_P->values()[p] * d_y[c];
    }
  }

  // backward substitution
  //   for i = m-1:0
  //     y[i] = 1/U[i,i] * ( b[i] - sum(k=i+1:m-1, U[i,k]*y[k]) )
  for (int i = d_P->number_rows() - 1; i >= 0; --i)
  {
    // diagonal index
    int d = d_P->diagonal(i);
    // end index
    int e = d_P->end(i);
    // invert row
    x[i] = d_y[i];
    for (int p = d; p < e; ++p)
    {
      // column index
      int c = d_P->column(p);
      x[i] -= d_P->values()[p] * x[c];
    }
    x[i] /= d_P->values()[d];
  }

}

} // end namespace detran

#endif // callow_PCILU0_I_HH_

//---------------------------------------------------------------------------//
//              end of file PCILU0.i.hh
//---------------------------------------------------------------------------//
