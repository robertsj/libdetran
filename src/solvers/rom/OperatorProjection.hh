/*
 * ProjectedOperator.hh
 *
 *  Created on: Jun 24, 2020
 *      Author: rabab
 */


#ifndef SOLVERS_ROM_OPERATORPROJECTION_HH_
#define SOLVERS_ROM_OPERATORPROJECTION_HH_

#include "solvers/solvers_export.hh"
#include "callow/matrix/MatrixBase.hh"
#include "callow/matrix/MatrixDense.hh"
#include "callow/matrix/MatrixShell.hh"
#include "callow/matrix/Matrix.hh"
#include "callow/vector/Vector.hh"

namespace detran
{
class OperatorProjection
{

  public:
	typedef callow::MatrixBase::SP_matrix  SP_matrix;
	typedef callow::MatrixDense::SP_matrix  SP_matrixDense;

	// constructor
	OperatorProjection(int a);

	void SetOperators(SP_matrix A, SP_matrix U);

	void Project(SP_matrixDense Ar);


  private:
	int d_r;
	int d_n;
	SP_matrix d_A;
	SP_matrixDense d_U;
	callow::MatrixDense ComputeAU();

};
}

#endif /* SOLVERS_ROM_OPERATORPROJECTION_HH_ */
