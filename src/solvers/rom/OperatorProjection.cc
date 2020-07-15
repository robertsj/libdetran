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
	callow::MatrixDense AU = OperatorProjection::ComputeAU();
	// compute UTAU
	for (int i=0; i<d_r; i++)
	 {
	  for (int k=0; k<d_r ; k++)
	  {
	    for (int j=0; j<d_n; j++)
		{
		 double v = (*d_U)(j, i)*AU(j, k);
		 Ar->insert(i, k, v, 1);
		}
	  }
   }
}


callow::MatrixDense OperatorProjection::ComputeAU()
{
  callow::Vector y(d_n, 0.0);
  callow::MatrixDense AU(d_n, d_r);

  for (int i=0; i<d_r; i++)
  {
	callow::Vector v(d_n);
	for (int j=0; j<d_n; j++)
	{
	  v[j] = (*d_U)(j, i);
	}

    d_A->multiply(v, y);

    double *y_ = &y[0];

    AU.insert_col(i, y_, 0);
  }
  //AU.print_matlab("AU.txt");

  return AU;
}
}


