//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Operator.i.hh
 * \brief  Operator inline member definitions
 * \author Jeremy Roberts
 * \date   Sep 8, 2012
 */
//---------------------------------------------------------------------------//

#ifndef detran_OPERATOR_I_HH_
#define detran_OPERATOR_I_HH_

namespace detran
{

inline void Operator::multiply(Vector &x, Vector &y)
{
  // Preconditions
  Require(is_assembled());
  Require(x.is_assembled());
  Require(x.size() == number_columns());
  PetscErrorCode ierr = MatMult(d_A, x.V(), y.V());
  // Postconditions
  Ensure(!ierr);
}

inline void Operator::multiply_transpose(Vector &x, Vector &y)
{
  // Preconditions
  Require(is_assembled());
  Require(x.is_assembled());
  Require(x.size() == number_rows());
  PetscErrorCode ierr = MatMultTranspose(d_A, x.V(), y.V());
  // Postconditions
  Ensure(!ierr);
}

inline Mat Operator::A()
{
  return d_A;
}

inline Operator::size_t Operator::number_rows() const
{
  return d_number_rows;
}

inline Operator::size_t Operator::number_columns() const
{
  return d_number_columns;
}

inline bool Operator::is_assembled() const
{
  return d_is_assembled;
}

} // end namespace detran

#endif // detran_OPERATOR_I_HH_

//---------------------------------------------------------------------------//
//              end of file Operator.i.hh
//---------------------------------------------------------------------------//
