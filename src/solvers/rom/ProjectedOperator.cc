/*
 * ProjectedOperator.cc
 *
 *  Created on: Jun 24, 2020
 *      Author: rabab
 */


#include "ProjectedOperator.hh"

namespace detran
{

template <class T>
ProjectedOperator<T>::ProjectedOperator(int a)
:d_r(0)
,d_n(0)
{
 std::cout << "here we go" << "\n";
 //build();
}

template <class T>
void ProjectedOperator<T>::SetOperators(SP_matrix A, SP_matrix U)
{
	Require(A);
	Require(U);

    d_A = A;
    d_U = U;

    Ensure(d_A->number_rows() == d_A->number_columns());
    Ensure(d_A->number_rows() == d_U->number_columns());
    d_n = d_U->number_rows();
    d_r = d_U->number_columns();

}



template <class T>
void ProjectedOperator<T>::Project(SP_matrix Ar)
{
	callow::MatrixDense A_ = ProjectedOperator<T>::ComputeAU();
	// compute UTAU
	for (int i=0; i<d_r; i++)
	 {
	  for (int k=0; k<d_r ; k++)
	  {
	    for (int j=0; j<d_n; j++)
		{
		 double v = (*d_U)(j, i)*A_(j, k);
		  //std::cout << U[j][i] << " & ";
		 Ar->insert(i, k, v, 1);
		}
	  }

	 }
}


template <class T>
callow::MatrixDense ProjectedOperator<T>::ComputeAU()
{
  callow::Vector y(d_n, 0.0);
  callow::MatrixDense A_(d_n, d_r);

  for (int i=0; i<d_r; i++)
  {
	callow::Vector v(d_n);
	for (int j=0; j<d_n; j++)
	{
	  v[j] = (*d_U)(i, j);
	}
    d_A->multiply(v, y);

    // need to cast from callow vector to double
    double *a = &v[0];
    A_.insert_col(i, a, 0);
  }
  return A_;
}



SOLVERS_INSTANTIATE_EXPORT(ProjectedOperator<callow::MatrixDense>)
}




