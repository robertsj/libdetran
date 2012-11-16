//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Precursors.cc
 *  @brief  Precursors
 *  @author Jeremy Roberts
 *  @date   Nov 15, 2012
 */
//---------------------------------------------------------------------------//

#include "Precursors.hh"

namespace detran
{

Precursors::Precursors(const size_t number_precursor_groups,
                       const size_t number_cells)
  : d_number_precursor_groups(number_precursor_groups)
  , d_number_cells(number_cells)
  , d_C(d_number_precursor_groups, vec_dbl(d_number_cells, 0.0))
{

}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file Precursors.cc
//---------------------------------------------------------------------------//
