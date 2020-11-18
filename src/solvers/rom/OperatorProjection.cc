//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  OperatorProjection.cc
 *  @brief OperatorProjection class definition
 *  @note  Copyright(C) 2020 Jeremy Roberts
 */
//----------------------------------------------------------------------------//



#include "OperatorProjection.hh"

namespace detran
{

OperatorProjection::OperatorProjection(int a)
:d_r(0)
,d_n(0)
{

}

void OperatorProjection::SetOperators(SP_matrix A, SP_matrix U)
{
  Require(A);
  Require(U);

  d_A = A;
  d_U = U;

  d_n = d_U->number_rows();
  d_r = d_U->number_columns();

  Ensure(d_A->number_rows() == d_A->number_columns());
  Ensure(d_A->number_rows() == d_U->number_rows());
  // make sure that the rank is less than problem size.
  Ensure(d_n > d_r);
}

void OperatorProjection::Project(SP_matrixDense Ar)
{
  // compute AU
  Vector y(d_n, 0.0);
  Vector y2(d_r, 0.0);

  Vector v(d_n);
  for (int i=0; i<d_r; i++)
  {
   Vector v(d_n, &(*d_U)(0, i));
   // compute AU
  d_A->multiply(v, y);
  // compute UTAU
  d_U->multiply_transpose(y, y2);
  double *y_ = &y2[0];
  Ar->insert_col(i, y_, 0);
 }
} // end namespace detran

} // end namespace detran
