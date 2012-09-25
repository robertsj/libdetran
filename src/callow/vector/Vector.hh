//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Vector.hh
 * \author robertsj
 * \date   Sep 13, 2012
 * \brief  Vector class definition.
 */
//---------------------------------------------------------------------------//

#ifndef callow_VECTOR_HH_
#define callow_VECTOR_HH_

#include "callow/callow_config.hh"
#include "callow/utils/CallowDefinitions.hh"
#include "utilities/SP.hh"
#include <vector>

namespace callow
{

/*!
 *  \class Vector
 *  \brief Dense vector object
 *
 *
 */
template <class T>
class Vector
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<Vector<T> >    SP_vector;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /// Various constructor interfaces
  /// @{
  Vector();
  Vector(const int n, T v = 0);
  Vector(const Vector &x);
  Vector(Vector &x);
  Vector(std::vector<T> &x);
  //Vector(SP_vector x);
#ifdef CALLOW_ENABLE_PETSC
  Vector(Vec pv);
#endif
  /// @}

  /// SP constructor
  static SP_vector
  Create(const int n, const T v = 0.0)
  {
    SP_vector p(new Vector<T>(n, v));
    return p;
  }

  /// Virtual destructor
  virtual ~Vector();

  /// Wipe out the contents and resize
  void resize(const int n, const T v = 0.0);

  //-------------------------------------------------------------------------//
  // ACCESS
  //-------------------------------------------------------------------------//

  const T& operator[](const int i) const;
  T& operator[](const int i);
  const T& operator()(const int i) const;
  T& operator()(const int i);
  const T& value(const int i) const;
  T& value(const int i);
#ifdef CALLOW_ENABLE_PETSC
  Vec petsc_vector() {return d_petsc_vector;}
#endif

  //-------------------------------------------------------------------------//
  // VECTOR OPERATIONS
  //-------------------------------------------------------------------------//

  /// Inner product of this vector with vector x
  T dot(const Vector<T>& x);
  T dot(SP_vector x);
  /// Norm of this vector
  T norm(const int type = L2);
  /**
   *  @brief Norm of the difference of this and another
   *
   *  Relative norms are with respect to this vector, and zeros
   *  are not checked.
   */
  T norm_residual(const Vector<T>& x, const int type = L2);
  T norm_residual(SP_vector x, const int type = L2);

  /// Set all elements of this vector to a value v
  void set(const T v);
  /// Multiply all elements of this vector by a value v
  void scale(const T v);
  /// Add a vector x to this vector
  void add(const Vector<T>& x);
  void add(SP_vector x);
  /// Subtract a vector x from this vector
  void subtract(const Vector<T>& x);
  void subtract(SP_vector x);
  /// Multiply this vector pointwise with a vector x
  void multiply(const Vector<T>& x);
  void multiply(SP_vector x);
  /// Multiply this vector pointwise with a vector x
  void divide(const Vector<T>& x);
  void divide(SP_vector x);
  /// Copy a vector x to this vector
  void copy(const Vector<T>& x);
  void copy(SP_vector x);
  /// Add a vector x times a scalar a to this vector
  void add_a_times_x(const T a, const Vector<T>& x);
  void add_a_times_x(const T a, SP_vector x);

  //-------------------------------------------------------------------------//
  // QUERY
  //-------------------------------------------------------------------------//

  int size() const { return d_size; }

  //-------------------------------------------------------------------------//
  // IO
  //-------------------------------------------------------------------------//

  /// Pretty print to stdout
  void display() const;
  /// Formatted write to ascii for matlab
  void print_matlab(std::string filename="vector.out") const;

protected:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Size of the vector
  int d_size;
  /// Values of the vector
  T* d_value;
#ifdef CALLOW_ENABLE_PETSC
  /// PETSc vector object
  Vec d_petsc_vector;
#endif
  // Is this a temporary wrapper around a pointer? i.e. no delete at dtor?
  bool d_temporary;

};

} // end namespace callow

#include "Vector.i.hh"

#endif /* callow_VECTOR_HH_ */
