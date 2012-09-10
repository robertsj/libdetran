//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   OperatorMatrix.hh
 * \brief  OperatorMatrix class definition
 * \author Jeremy Roberts
 * \date   Sep 8, 2012
 */
//---------------------------------------------------------------------------//

#ifndef detran_OPERATORMATRIX_HH_
#define detran_OPERATORMATRIX_HH_

#include "Operator.hh"

namespace detran
{

/*!
 *  \class OperatorMatrix
 *  \brief An operator with an explicit matrix representation
 *
 *
 */
class OperatorMatrix: public Operator
{

public:

  //---------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //---------------------------------------------------------------------------//

  /*!
   *  \brief Constructor
   *  \param m    number of rows
   *  \param n    number of columns
   */
  OperatorMatrix(const size_t m, const size_t n);

  // Virtual destructor
  virtual ~OperatorMatrix(){}

  //---------------------------------------------------------------------------//
  // ABSTRACT INTERFACE
  //---------------------------------------------------------------------------//

  /// Assemble the matrix.
  void assemble();

  /*!
   *  \brief Display the matrix to screen (or to output)
   *  \param  output    Flag indicating (stdout=0, ascii=1, binary=2)
   *  \param  name      File name for ascii or binary file
   */
  virtual void display(const int output = 0,
                       const std::string name = "matrix.out") const;


protected:

  //---------------------------------------------------------------------------//
  // IMPLEMENTATION
  //---------------------------------------------------------------------------//

  /*!
   *  \brief Preallocate the matrix
   *  \param number_nonzeros    Number of nonzeros per row
   */
  void preallocate(vec_int &number_nonzeros);

  /*!
   *  \brief Insert values (=)
   *
   *  \param number_rows    Number of column indices
   *  \param rows           Indices of global rows
   *  \param number_columns Number of column indices
   *  \param columns        Indices of global columns
   *  \param values         Logically 2D array of values to insert
   *  \param insert_t       Insertion type
   */
   void insert_values(const size_t number_rows,
                      const int *rows,
                      const size_t number_columns,
                      const int *columns,
                      const double *values,
                      InsertMode = INSERT_VALUES);

    /// Build the underlying matrix.
    virtual void build() = 0;

private:

    /// INSERT or ADD as last insert type
    InsertMode d_insert_mode;

    /// Flush the matrix insertion type
    void flush(InsertMode insert_t);
};

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE MEMBERS
//---------------------------------------------------------------------------//

#include "OperatorMatrix.i.hh"

#endif // detran_OPERATORMATRIX_HH_

//---------------------------------------------------------------------------//
//              end of file OperatorMatrix.hh
//---------------------------------------------------------------------------//
