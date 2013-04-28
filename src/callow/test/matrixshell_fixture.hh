//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   matrixshell_fixture.hh
 * \author robertsj
 * \date   Sep 20, 2012
 * \brief  matrixshell_fixture class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef MATRIXSHELL_FIXTURE_HH_
#define MATRIXSHELL_FIXTURE_HH_

#include "matrix/Matrix.hh"
#include "matrix/MatrixShell.hh"
#include "test/matrix_fixture.hh"

namespace callow
{

//---------------------------------------------------------------------------//
class TestMatrixShell: public MatrixShell
{
public:
  TestMatrixShell(const size_t m)
    : MatrixShell(this, m, m),
      d_B(test_matrix_1(m))
  {
    d_B->display();
    d_B->print_matlab("test.out");
  }
  // THESE MUST BE IMPLEMENTED BY SHELL OPERATORS
  void multiply(const Vector &x,  Vector &y)
  {
    d_B->multiply(x, y);
  }
  void multiply_transpose(const Vector &x,  Vector &y)
  {
    d_B->multiply_transpose(x, y);
  }
private:
  // Regular matrix.
  Matrix::SP_matrix d_B;
};

} // end namespace callow

#endif /* MATRIXSHELL_FIXTURE_HH_ */
