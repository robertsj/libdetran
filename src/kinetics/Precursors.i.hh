//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Precursors.i.hh
 *  @brief  Precursors.i
 *  @author Jeremy Roberts
 *  @date   Nov 15, 2012
 */
//---------------------------------------------------------------------------//

#ifndef detran_PRECURSORS_I_HH_
#define detran_PRECURSORS_I_HH_

namespace detran
{

//---------------------------------------------------------------------------//
inline const Precursors::vec_dbl&
Precursors::C(const size_t i) const
{
  Require(i < d_number_precursor_groups);
  return d_C[i];
}

//---------------------------------------------------------------------------//
inline Precursors::vec_dbl&
Precursors::C(const size_t i)
{
  // Cast away return type
  return const_cast<vec_dbl&>
  (
    // Add const to *this's type and call const version
    static_cast<const Precursors*>(this)->C(i)
  );
}



} // end namespace detran

#endif // detran_PRECURSORS_I_HH_

//---------------------------------------------------------------------------//
//              end of file Precursors.i.hh
//---------------------------------------------------------------------------//
