//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   OperatorMatrix.i.hh
 * \brief  OperatorMatrix inline member definitions
 * \author Jeremy Roberts
 * \date   Sep 8, 2012
 */
//---------------------------------------------------------------------------//

#ifndef detran_OPERATORMATRIX_I_HH_
#define detran_OPERATORMATRIX_I_HH_

namespace detran
{

inline void OperatorMatrix::
insert_values(const size_t number_rows,
              const int *rows,
              const size_t number_columns,
              const int *columns,
              const double *values,
              InsertMode insert_t)
{
  // Preconditions
  Require(number_rows > 0);
  Require(number_columns > 0);

  // Flush if needed
  flush(insert_t);

  d_is_assembled = false; // changing values = no longer ready

  PetscErrorCode ierr;
  ierr = MatSetValues(d_A, number_rows, rows, number_columns, columns,
                      values, insert_t);

  // Postconditions
  Ensure(!ierr);
}


} // end namespace detran


#endif // detran_OPERATORMATRIX_I_HH_

//---------------------------------------------------------------------------//
//              end of file OperatorMatrix.i.hh
//---------------------------------------------------------------------------//
