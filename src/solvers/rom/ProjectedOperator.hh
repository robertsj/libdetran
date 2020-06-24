/*
 * ProjectedOperator.hh
 *
 *  Created on: Jun 24, 2020
 *      Author: rabab
 */

#ifndef SOLVERS_ROM_PROJECTEDOPERATOR_HH_
#define SOLVERS_ROM_PROJECTEDOPERATOR_HH_
#include "callow/matrix/MatrixBase.hh"
#include "callow/matrix/MatrixDense.hh"
#include "callow/matrix/MatrixShell.hh"
#include "callow/matrix/Matrix.hh"

class ProjectedOPerator
{

  public:
	ProjectedOPerator(callow::Matrix A, callow::MatrixDense U);
	//ProjectedOPerator(callow::MatrixShell A, ROMBasis U);
	//ProjectedOPerator(matrixshell A, ROMBasis U);

	callow::MatrixDense ComputeAU();

	callow::MatrixDense ComputeUTAU();


  private:
	callow::MatrixDense d_U;
	callow::Matrix d_A;
	int d_r;
	int d_n;


};



#endif /* SOLVERS_ROM_PROJECTEDOPERATOR_HH_ */
