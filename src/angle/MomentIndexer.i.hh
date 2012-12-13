//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   MomentIndexer.i.hh
 *  @author robertsj
 *  @date   Dec 13, 2012
 *  @brief  MomentIndexer.i class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_angle_MOMENTINDEXER_I_HH_
#define detran_angle_MOMENTINDEXER_I_HH_

namespace detran_angle
{

//---------------------------------------------------------------------------//
inline int MomentIndexer::l(const size_t i) const
{
   Require(i < d_number_moments);
   return d_l[i];
}

//---------------------------------------------------------------------------//
inline int MomentIndexer::m(const size_t i) const
{
  Require(i < d_number_moments);
  return d_m[i];
}

//---------------------------------------------------------------------------//
inline MomentIndexer::size_t
MomentIndexer::index(const size_t l, const int m) const
{
  Require(l <= d_legendre_order);
  Require(m <= l);
  Require(-m <= l);
  THROW("lalala");
  return 0;
}


//---------------------------------------------------------------------------//
inline const MomentIndexer::vec_int&
MomentIndexer::m_index(const size_t l) const
{
  Require(l <= d_legendre_order);
  return d_m_index[l];
}

} // end detran_angle

#endif /* detran_angle_MOMENTINDEXER_I_HH_ */
