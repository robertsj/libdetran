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

  //---------------------------------------------------------------------------//
  // ENUMERATIONS
  //---------------------------------------------------------------------------//

  enum vec_norm_types
  {
    L1, L2, LINF, L1GRID, L2GRID, L2REL, END_VEC_NORM_TYPES
  };

  //---------------------------------------------------------------------------//
  // TYPEDEFS
  //---------------------------------------------------------------------------//

  typedef detran_utilities::SP<Vector<T> >    SP_vector;

  //---------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //---------------------------------------------------------------------------//

  Vector();
  Vector(const int n, T v = 0);
  Vector(const Vector &x);
  Vector(Vector &x);
  Vector(std::vector<T> &x);
  //Vector(SP_vector x);
#ifdef CALLOW_ENABLE_PETSC
  Vector(Vec pv);
#endif

  /// SP constructor
  static SP_vector
  Create(const int n, const T v = 0.0)
  {
    SP_vector p(new Vector<T>(n, v));
    return p;
  }


  virtual ~Vector();
  void resize(const int n, const T v = 0.0);

  //---------------------------------------------------------------------------//
  // ACCESS
  //---------------------------------------------------------------------------//

  const T& operator[](const int i) const;
  T& operator[](const int i);
  const T& operator()(const int i) const;
  T& operator()(const int i);
  const T& value(const int i) const;
  T& value(const int i);
#ifdef CALLOW_ENABLE_PETSC
  Vec petsc_vector() {return d_petsc_vector;}
#endif

  //---------------------------------------------------------------------------//
  // VECTOR OPERATIONS
  //---------------------------------------------------------------------------//

  // scalar
  T dot(const Vector<T>& x);
  T dot(SP_vector x);
  T norm(const int type = L2);
  T norm_residual(const Vector<T>& x, const int type = L2);
  T norm_residual(SP_vector x, const int type = L2);

  // element-wise operations
  void set(const T v);
  void scale(const T v);
  //
  void add(const Vector<T>& x);
  void subtract(const Vector<T>& x);
  void multiply(const Vector<T>& x);
  void divide(const Vector<T>& x);
  void copy(const Vector<T>& x);
  void add_a_times_x(const T a, const Vector<T>& x);
  //
  void add(SP_vector x);
  void subtract(SP_vector x);
  void multiply(SP_vector x);
  void divide(SP_vector x);
  void copy(SP_vector x);
  void add_a_times_x(const T a, SP_vector x);

  //---------------------------------------------------------------------------//
  // QUERY
  //---------------------------------------------------------------------------//

  int size() const { return d_size; }

  //---------------------------------------------------------------------------//
  // IO
  //---------------------------------------------------------------------------//

  /// pretty print to stdout
  void display() const;
  /// formatted write to ascii for matlab
  void print_matlab(std::string filename="vector.out") const;

protected:

  //---------------------------------------------------------------------------//
  // DATA
  //---------------------------------------------------------------------------//

  int d_size;
  T* d_value;
#ifdef CALLOW_ENABLE_PETSC
  Vec d_petsc_vector;
#endif
  // is this a temporary wrapper around a pointer? i.e. no delete?
  bool d_temporary;

};

} // end namespace callow

#include "Vector.i.hh"

#endif /* callow_VECTOR_HH_ */
