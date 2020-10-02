//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  OperatorProjection.hh
 *  @brief OperatorProjection class definition
 *  @note  Copyright(C) 2020 Jeremy Roberts
 */
//----------------------------------------------------------------------------//


#ifndef SOLVERS_ROM_OPERATORPROJECTION_HH_
#define SOLVERS_ROM_OPERATORPROJECTION_HH_

#include "solvers/solvers_export.hh"
#include "callow/matrix/MatrixBase.hh"
#include "callow/matrix/MatrixDense.hh"
#include "callow/matrix/MatrixShell.hh"
#include "callow/matrix/Matrix.hh"
#include "callow/vector/Vector.hh"

using namespace callow;

namespace detran
{
 /**
 *  @class OperatorProjection
 *  @brief project matrix onto the space of a given basis
 */
class OperatorProjection
{

  public:
	typedef MatrixBase::SP_matrix  SP_matrix;
	typedef MatrixDense::SP_matrix  SP_matrixDense;

	/// constructor
	OperatorProjection(int a);
    /// set the operator pointer and the basis
	void SetOperators(SP_matrix A, SP_matrix U);
    /// generate the reduced matrix UTAU
	void Project(SP_matrixDense Ar);


  private:
	/// basis rank
	int d_r;
	/// operator size
	int d_n;
	/// operator to be projected
	SP_matrix d_A;
	/// basis onto which the operator is projected
	SP_matrixDense d_U;

};

} // end namespace detran

#endif /* SOLVERS_ROM_OPERATORPROJECTION_HH_ */
