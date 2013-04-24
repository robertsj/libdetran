//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Precursors.cc
 *  @brief  Precursors
 *  @author Jeremy Roberts
 *  @date   Nov 15, 2012
 */
//---------------------------------------------------------------------------//

#include "Precursors.hh"
#include <cstdio>

namespace detran
{

//---------------------------------------------------------------------------//
Precursors::Precursors(const size_t number_precursor_groups,
                       const size_t number_cells)
  : d_number_precursor_groups(number_precursor_groups)
  , d_number_cells(number_cells)
  , d_C(d_number_precursor_groups, vec_dbl(d_number_cells, 0.0))
{
  /* ... */
}

//---------------------------------------------------------------------------//
void Precursors::display() const
{

  for (size_t i = 0; i < d_number_cells; i++)
  {
    printf("%10i", i);
    for (size_t g = 0; g < d_number_precursor_groups; g++)
    {
      printf(" %12.5e ", d_C[g][i]);
    }
    printf("\n");
  }

}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file Precursors.cc
//---------------------------------------------------------------------------//
