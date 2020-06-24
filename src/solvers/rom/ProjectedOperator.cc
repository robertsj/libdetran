/*
 * ProjectedOperator.cc
 *
 *  Created on: Jun 24, 2020
 *      Author: rabab
 */


#include "ProjectedOperator.hh"
#include "callow/vector/Vector.hh"

ProjectedOPerator::ProjectedOPerator(callow::Matrix A, callow::MatrixDense U)
:d_U(U),
 d_A(A)
{
 d_r = U.number_columns;
 d_n = U.number_rows();
}

/*
vector<vector<double>> basis_vecs;

 for (int i=0; i<r; i++)
 {
   vector<double> v;
   for (int j=0; j< n; j++)
   {
     v.push_back(U[j][i]);
   }
     basis_vecs.push_back(v);
 }
 */

callow::MatrixDense ProjectedOPerator::ComputeAU()
{

  callow::Vector y(d_n, 0.0);
  callow::MatrixDense A_(d_n, d_r);
  for (int i=0; i<d_r; i++)
  {
	callow::Vector v(d_n);
	for (int j=0; j<d_n; j++)
	{
	  v[j] = U[i][j];
	}
    d_A.multiply(v, y);
    // need to cast from callow vector to double
    A_.insert_col(i, y, 0);
  }
  return A_;
}


callow::MatrixDense ProjectedOPerator::ComputeUTAU()
{
 callow::MatrixDense::SP_matrix A_r;
 A_r = new callow::MatrixDense(d_r, d_r);
 for (int i=0; i<d_r; i++)
 {
  for (int k=0; k<d_r ; k++)
  {
    for (int j=0; j<d_n; j++)
	{
	  double v = U(j, i)*A_r1(j, k);
	  //std::cout << U[j][i] << " & ";
	  A_r->insert(i, k, v, 1);
	}
  }

 }
 return A_r;
}


