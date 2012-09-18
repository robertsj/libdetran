//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Vector.hh
 * \brief  Vector class definition
 * \author Jeremy Roberts
 * \date   Sep 8, 2012
 */
//---------------------------------------------------------------------------//

#ifndef detran_VECTOR_HH_
#define detran_VECTOR_HH_

#include "utilities/Definitions.hh"
#include "utilities/SP.hh"
#include "petsc.h"

namespace detran
{

/*!
 *  \class Vector
 *  \brief Lightweight wrapper for PETSc Vec.
 *
 *  Use of Vector and the corresponding \ref Operator class should, in
 *  theory, eliminate a lot of PETSc code from the rest of Detran.
 *
 */
/*!
 *  \example solvers/test/test_Vector.cc
 *
 *  Test of Vector class.
 */

class Vector
{

public:

  //---------------------------------------------------------------------------//
  // ENUMERATIONS
  //---------------------------------------------------------------------------//

  enum vec_norm_types
  {
    L1, L2, LINF, END_VEC_NORM_TYPES
  };

  //---------------------------------------------------------------------------//
  // TYPEDEFS
  //---------------------------------------------------------------------------//

  typedef detran_utilities::SP<Vector>      SP_vector;
  typedef detran_utilities::size_t          size_t;

  //---------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //---------------------------------------------------------------------------//

  /*!
   *  \brief Constructor
   *  \param m      Local number of rows
   *  \param val    Optional initial value
   */
  Vector(const size_t m, const double val = 0.0);

  /// Destructor
  ~Vector();

  //---------------------------------------------------------------------------//
  // SETTERS
  //---------------------------------------------------------------------------//

  /*!
   *  \brief Insert values
   *  \param values   Array of values to insert
   *  \param number   Number of values to insert
   *  \param rows     Indices of rows where values are inserted
   */
  void insert_values(const unsigned int number,
                     const int *rows,
                     const double *values);


  /// Assemble the vector.
  void assemble();

  //---------------------------------------------------------------------------//
  // VECTOR OPERATIONS
  //---------------------------------------------------------------------------//

  /// Dot product of another Vector with me
  double dot(Vector &x);

  /// Scale the Vector
  void scale(const double factor);

  /// Add a vector to me
  void add(Vector &x);

  /// Subtract a vector from me
  void subtract(Vector &x);

  /// Residual norm
  double residual_norm(Vector &x, const int type);

  //---------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //---------------------------------------------------------------------------//

  /*
   *  \brief Const access to local array
   *  \param i  Local index
   */
  const double& operator[](const size_t i) const;

  /*
   *  \brief Mutable access to local array
   *  \param i  Local index
   */
  double& operator[](const size_t i);

  /// Return the PETSc vector
  Vec V() const
  {
    return d_V;
  }

  /// Return the local size.
  size_t size() const
  {
    return d_size;
  }

  /// Return assembled flag.
  bool is_assembled() const
  {
    return d_is_assembled;
  }

  /// View via standard output.
  void display() const;

private:

  //---------------------------------------------------------------------------//
  // PRIVATE DATA
  //---------------------------------------------------------------------------//

  /// PETSc vector
  Vec d_V;

  /// Pointer to underlying local array.  Be *careful* when using.
  double* d_array;

  /// Local size
  const size_t d_size;

  /// Am I assembled?
  bool d_is_assembled;

};

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "Vector.i.hh"

#endif // detran_VECTOR_HH_

//---------------------------------------------------------------------------//
//              end of file Vector.hh
//---------------------------------------------------------------------------//
