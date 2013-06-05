//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  Vector.hh
 *  @brief Vector class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
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

/**
 *  @class Vector
 *  @brief Dense vector object
 */
class CALLOW_EXPORT Vector
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<Vector>    SP_vector;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /// Various constructor interfaces
  /// @{

  /// Default constructor
  Vector();
  Vector(const int n, double v = 0);
  Vector(const Vector &x);
  //Vector(Vector &x);
  /// Temporarily wrap a std::vector
  Vector(const std::vector<double> &x);
  /// Temporarily wrap a dumb point
  Vector(const int n, double* v);
  //Vector(SP_vector x);
#ifdef CALLOW_ENABLE_PETSC
  Vector(Vec pv);
#endif

  /// @}

  /// SP constructor
  static SP_vector
  Create(const int n, const double v = 0.0)
  {
    SP_vector p(new Vector(n, v));
    return p;
  }

  /// Virtual destructor
  virtual ~Vector();

  /// Wipe out the contents and resize
  void resize(const int n, const double v = 0.0);

  //-------------------------------------------------------------------------//
  // ACCESS
  //-------------------------------------------------------------------------//

  const double& operator[](const int i) const;
  double& operator[](const int i);
  const double& operator()(const int i) const;
  double& operator()(const int i);
  const double& value(const int i) const;
  double& value(const int i);
#ifdef CALLOW_ENABLE_PETSC
  Vec petsc_vector() {return d_petsc_vector;}
#endif

  //-------------------------------------------------------------------------//
  // VECTOR OPERATIONS
  //-------------------------------------------------------------------------//

  /// Inner product of this vector with vector x
  double dot(const Vector& x);
  double dot(SP_vector x);
  /// Norm of this vector
  double norm(const int type = L2);
  /**
   *  @brief Norm of the difference of this and another
   *
   *  Relative norms are with respect to this vector, and zeros
   *  are not checked.
   */
  double norm_residual(const Vector& x, const int type = L2);
  double norm_residual(SP_vector x, const int type = L2);

  /// Set all elements of this vector to a value v
  void set(const double v);
  /// Multiply all elements of this vector by a value v
  void scale(const double v);
  /// Add a vector x to this vector
  void add(const Vector& x);
  void add(SP_vector x);
  /// Subtract a vector x from this vector
  void subtract(const Vector& x);
  void subtract(SP_vector x);
  /// Multiply this vector pointwise with a vector x
  void multiply(const Vector& x);
  void multiply(SP_vector x);
  /// Multiply this vector pointwise with a vector x
  void divide(const Vector& x);
  void divide(SP_vector x);
  /// Copy a vector x to this vector
  void copy(const Vector& x);
  void copy(SP_vector x);
  /// Add a vector x times a scalar a to this vector
  void add_a_times_x(const double a, const Vector& x);
  void add_a_times_x(const double a, SP_vector x);

  //-------------------------------------------------------------------------//
  // QUERY
  //-------------------------------------------------------------------------//

  int size() const { return d_size; }

  //-------------------------------------------------------------------------//
  // IO
  //-------------------------------------------------------------------------//

  /// Pretty print to stdout
  void display(const std::string &name = "") const;
  /// Formatted write to ascii for matlab
  void print_matlab(std::string filename="vector.out") const;

protected:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Size of the vector
  int d_size;
  /// Values of the vector
  double* d_value;
#ifdef CALLOW_ENABLE_PETSC
  /// PETSc vector object
  Vec d_petsc_vector;
#endif
  // Is this a temporary wrapper around a pointer? i.e. no delete at dtor?
  bool d_temporary;
  // Is this also temporary around an extant PETSC vector?
  bool d_temporary_petsc;

};

CALLOW_TEMPLATE_EXPORT(detran_utilities::SP<Vector>);

} // end namespace callow

#include "Vector.i.hh"

#endif /* callow_VECTOR_HH_ */
