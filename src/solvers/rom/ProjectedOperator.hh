/*
 * ProjectedOperator.hh
 *
 *  Created on: Jun 24, 2020
 *      Author: rabab
 */


#ifndef SOLVERS_ROM_PROJECTEDOPERATOR_HH_
#define SOLVERS_ROM_PROJECTEDOPERATOR_HH_

#include "solvers/solvers_export.hh"
#include "callow/matrix/MatrixBase.hh"
#include "callow/matrix/MatrixDense.hh"
#include "callow/matrix/MatrixShell.hh"
#include "callow/matrix/Matrix.hh"
#include "callow/vector/Vector.hh"

namespace detran
{
template <class T>
class ProjectedOperator
{

  public:
	typedef detran_utilities::SP<T>  SP_matrix;

	// constructor
	ProjectedOperator();

	void SetOperators(SP_matrix A, SP_matrix U);

	void Project(SP_matrix Ar);


  private:
	int d_r;
	int d_n;
	SP_matrix d_A;
	SP_matrix d_U;
	callow::MatrixDense ComputeAU();

};
}

#endif /* SOLVERS_ROM_PROJECTEDOPERATOR_HH_ */
