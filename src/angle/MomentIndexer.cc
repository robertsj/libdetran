//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  MomentIndexer.cc
 *  @brief MomentIndexer class definition
 *  @note  Copyright (C) Jeremy Roberts 2012-2013
 */
//---------------------------------------------------------------------------//

#include "MomentIndexer.hh"

namespace detran_angle
{

//---------------------------------------------------------------------------//
MomentIndexer::MomentIndexer(const size_t dimension,
                             const size_t legendre_order)
  : d_legendre_order(legendre_order)
{
  // Preconditions
  Require(dimension >= 1);
  Require(dimension <= 3);

  if (dimension == 1)
    construct_1D();
  else if (dimension == 2)
    construct_2D();
  else
    construct_3D();

  // Fill the l and m vectors
  d_l.resize(d_number_moments, 0);
  d_m.resize(d_number_moments, 0);
  size_t j = 0;
  for (size_t l = 0; l <= d_legendre_order; ++l)
  {
    for (size_t i = 0; i < d_m_index[l].size(); ++i, ++j)
    {
      d_l[j] = l;
      d_m[j] = d_m_index[l][i];
    }
  }
}

//---------------------------------------------------------------------------//
MomentIndexer::SP_momentindexer
MomentIndexer::Create(const size_t dimension,
                      const size_t legendre_order)
{
  SP_momentindexer p(new MomentIndexer(dimension, legendre_order));
  return p;
}

//---------------------------------------------------------------------------//
void MomentIndexer::display() const
{
  std::cout << " MOMENTS INDEXER: " << std::endl;
  for (size_t i = 0; i < d_number_moments; ++i)
  {
    std::cout << " " << i << " " << l(i) << " " << m(i) << std::endl;
  }
}

//---------------------------------------------------------------------------//
void MomentIndexer::construct_1D()
{
  d_number_moments = d_legendre_order + 1;
  d_m_index.resize(d_legendre_order + 1);
  for (size_t l = 0; l <= d_legendre_order; ++l)
    d_m_index[l].resize(1, 0); // always zero
}

//---------------------------------------------------------------------------//
void MomentIndexer::construct_2D()
{
  d_number_moments = (d_legendre_order + 1) * (d_legendre_order + 2) / 2;
  d_m_index.resize(d_legendre_order + 1);
  for (int l = 0; l <= (int)d_legendre_order; ++l)
    for (int m = -l; m <= l; m++)
      if ((m + l) % 2 == 0) d_m_index[l].push_back(m);
}

//---------------------------------------------------------------------------//
void MomentIndexer::construct_3D()
{
  d_number_moments = (d_legendre_order + 1) * (d_legendre_order + 1);
  d_m_index.resize(d_legendre_order + 1);
  for (int l = 0; l <= (int)d_legendre_order; ++l)
  {
    for (int m = -l; m <= l; ++m)
      d_m_index[l].push_back(m);
  }
}

} // end namespace detran_angle

