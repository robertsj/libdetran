//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Operator.hh
 * \brief  Operator class definition
 * \author Jeremy Roberts
 * \date   Sep 8, 2012
 */
//---------------------------------------------------------------------------//

#ifndef detran_OPERATOR_HH_
#define detran_OPERATOR_HH_

#include "Vector.hh"
#include "petsc.h"

namespace detran
{

/*!
 *  \class Operator
 *  \brief Base operator class based on PETSc matrix
 */
class Operator
{

public:

  //---------------------------------------------------------------------------//
  // ENUMERATIONS
  //---------------------------------------------------------------------------//

  /// Output flag options for matrix viewer
  enum MATOUT
  {
    STDOUT, ASCII, BINARY
  };

  //---------------------------------------------------------------------------//
  // TYPEDEFS
  //---------------------------------------------------------------------------//

  typedef detran_utilities::SP<Operator>  SP_operator;
  typedef detran_utilities::size_t        size_t;
  typedef detran_utilities::vec_int       vec_int;

  //---------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //---------------------------------------------------------------------------//

  /*!
   *  \brief Constructor
   *  \param m local number of rows
   *  \param n local number of columns
   */
  Operator(const size_t m,
           const size_t n);

  /// Virtual destructor
  virtual ~Operator(){};

  //---------------------------------------------------------------------------//
  // ABSTRACT INTERFACE
  //---------------------------------------------------------------------------//

  /// Assemble the matrix.
  virtual void assemble() = 0;

  /*!
   *  \brief Display the matrix to screen (or to output)
   *  \param  output    Flag indicating (stdout=0, ascii=1, binary=2)
   *  \param  name      File name for ascii or binary file
   */
  virtual void display(const int output = 0,
                       const std::string name = "matrix.out") const = 0;

  //---------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //---------------------------------------------------------------------------//

  /*!
   *  \brief Matrix-vector multiplication
   *  \param x  Input vector
   *  \param y  Output vector
   */
  void multiply(Vector &x, Vector &y);

  /*!
   *  \brief Matrix-vector multiplication using matrix transpose.
   *  \param x  Input vector
   *  \param y  Output vector
   */
  void multiply_transpose(Vector &x, Vector &y);

  /// Get PETSc Mat object.
  Mat A();

  /// Get number of rows
  size_t number_rows() const;

  /// Get number of columns
  size_t number_columns() const;

  /// Am I assembled?
  bool is_assembled() const;

protected:

  //---------------------------------------------------------------------------//
  // DATA
  //---------------------------------------------------------------------------//

  /// PETSc matrix
  Mat d_A;

  /// Number of matrix rows on this process
  const size_t d_number_rows;

  /// Number of matrix columns on this process
  const size_t d_number_columns;

  /// Flag for whether the matrix is assembled
  bool d_is_assembled;

};

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE MEMBERS
//---------------------------------------------------------------------------//

#include "Operator.i.hh"

#endif // detran_OPERATOR_HH_

//---------------------------------------------------------------------------//
//              end of file Operator.hh
//---------------------------------------------------------------------------//
