//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   OrthogonalBasis.i.hh
 *  @brief  OrthogonalBasis inline member definitions
 *  @author Jeremy Roberts
 *  @date   Jan 8, 2013
 */
//---------------------------------------------------------------------------//

#ifndef detran_orthog_ORTHOGONALBASIS_I_HH_
#define detran_orthog_ORTHOGONALBASIS_I_HH_

namespace detran_orthog
{

//---------------------------------------------------------------------------//
inline double OrthogonalBasis::fold(const size_t row, const Vector &V)
{
  // Preconditions
  Require(V.size() == d_basis->number_columns());

  // Copy and vector and multiply point-wise by the weights
  Vector Y(V);
  if (d_w) Y.multiply(d_w);

  // Extract the basis row and return the dot product
  Vector P_l(d_size, &(*d_basis)(row, 0));
  return P_l.dot(Y);
}
//---------------------------------------------------------------------------//
inline double OrthogonalBasis::fold(const size_t row, SP_vector V)
{
  return fold(row, *V);
}
//---------------------------------------------------------------------------//
inline double OrthogonalBasis::fold(const size_t row, const vec_dbl &V)
{
  Vector v(V);
  return fold(row, v);
}

//---------------------------------------------------------------------------//
inline double OrthogonalBasis::unfold(const size_t column, const Vector &V)
{
  // Preconditions
  Require(V.size() == d_basis->number_rows());

  // Copy and vector and multiply by the inverse normalization coefficients
  Vector Y(V);
  if (d_a) Y.multiply(d_a);

  // Extract the basis column and return the dot product
  Vector P_l(d_basis->number_rows(), 0.0);
  for (size_t i = 0; i < P_l.size(); ++i)
    P_l = (*d_basis)(i, column);
  return P_l.dot(Y);
}
//---------------------------------------------------------------------------//
inline double OrthogonalBasis::unfold(const size_t column, SP_vector V)
{
  return unfold(column, *V);
}
//---------------------------------------------------------------------------//
inline double OrthogonalBasis::unfold(const size_t column, const vec_dbl &V)
{
  Vector v(V);
  return unfold(column, v);
}

//---------------------------------------------------------------------------//
inline void OrthogonalBasis::transform(const Vector &f, Vector &ftilde)
{
  // Copy and vector and multiply point-wise by the weights
  Vector wf(f);
  if (d_w) wf.multiply(d_w);
  d_basis->multiply(wf, ftilde);
}
//---------------------------------------------------------------------------//
inline void OrthogonalBasis::transform(SP_vector f, SP_vector ftilde)
{
  transform(*f, *ftilde);
}
//---------------------------------------------------------------------------//
inline void OrthogonalBasis::transform(const vec_dbl &f, vec_dbl &ftilde)
{
  Vector F(f);
  Vector Ftilde(ftilde);
  transform(F, Ftilde);
}

//---------------------------------------------------------------------------//
inline void OrthogonalBasis::inverse(const Vector &ftilde, Vector &f)
{
  // Copy and vector and multiply point-wise by the weights
  Vector aft(ftilde);
  if (d_a) aft.multiply(d_a);
  d_basis->multiply_transpose(aft, f);
}
//---------------------------------------------------------------------------//
inline void OrthogonalBasis::inverse(SP_vector ftilde, SP_vector f)
{
  inverse(*f, *ftilde);
}
//---------------------------------------------------------------------------//
inline void OrthogonalBasis::inverse(const vec_dbl &ftilde, vec_dbl &f)
{
  Vector F(f);
  Vector Ftilde(ftilde);
  inverse(Ftilde, F);
}

} // end namespace detran_orthog

#endif // detran_orthog_ORTHOGONALBASIS_I_HH_

//---------------------------------------------------------------------------//
//              end of file OrthogonalBasis.i.hh
//---------------------------------------------------------------------------//
