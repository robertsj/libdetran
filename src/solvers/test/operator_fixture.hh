//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   operator_fixture.hh
 * \brief  Test implementations of matrix and shell operators
 * \author Jeremy Roberts
 * \date   Sep 9, 2012
 */
//---------------------------------------------------------------------------//

#ifndef OPERATOR_FIXTURE_HH_
#define OPERATOR_FIXTURE_HH_

#include "OperatorMatrix.hh"
#include "OperatorShell.hh"

namespace detran
{

//---------------------------------------------------------------------------//
class TestOperatorMatrix: public OperatorMatrix
{

public:

  TestOperatorMatrix(size_t n)
    : OperatorMatrix(n, n)
  { /* ... */ }

private:

  // THIS MUST BE IMPLEMENTED BY MATRIX OPERATORS
  void build()
  {
    // Doing a small 1D second order difference
    for (int row = 0; row < number_rows(); row++)
    {
      if (row == 0)
      {
        int c[]    = {0, 1};
        double v[] = {-2, 1};
        insert_values(1, &row, 2, c, v);
      }
      else if (row == number_rows() - 1)
      {
        int    c[] = {number_rows() - 2, number_rows() - 1};
        double v[] = {1, -2};
        insert_values(1, &row, 2, c, v);
      }
      else
      {
        int    c[] = {row - 1, row, row + 1};
        double v[] = {1, -2, 1};
        insert_values(1, &row, 3, c, v);
      }
    }
    assemble();
  }
};

//---------------------------------------------------------------------------//
class TestOperatorShell: public OperatorShell
{

public:

  TestOperatorShell(const size_t m)
      : OperatorShell(m, m, this), d_B(m, m)
  {
    /* ... */
  }

  // THESE MUST BE IMPLEMENTED BY SHELL OPERATORS

  PetscErrorCode
  shell_multiply(Vec x, Vec y)
  {
    return MatMult(d_B.A(), x, y);
  }

  PetscErrorCode
  shell_multiply_transpose(Vec x, Vec y)
  {
    return MatMultTranspose(d_B.A(), x, y);
  }

private:

  // Regular matrix.
  TestOperatorMatrix d_B;

};

} // end namespace detran

#endif // OPERATOR_FIXTURE_HH_ 

//---------------------------------------------------------------------------//
//              end of file operator_fixture.hh
//---------------------------------------------------------------------------//
