/*
 * ProjectedOperator.cc
 *
 *  Created on: Jun 24, 2020
 *      Author: rabab
 */


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
    // maybe need to make sure that the rank is less than problem size.
    Ensure(d_n > d_r);
}

void OperatorProjection::Project(SP_matrixDense Ar)
{
  // compute AU
  callow::Vector y(d_n, 0.0);
  callow::Vector y2(d_r, 0.0);
  callow::MatrixDense AU(d_n, d_r);

  callow::Vector v(d_n);
  for (int i=0; i<d_r; i++)
  {
   for (int j=0; j<d_n; j++)
   {
    v[j] = (*d_U)(j, i);
   }

  d_A->multiply(v, y);
  d_U->multiply_transpose(y, y2);
  double *y_ = &y2[0];
  Ar->insert_col(i, y_, 0);
 }
}

}

