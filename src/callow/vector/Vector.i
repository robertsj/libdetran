//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Vector.i
 * \author Jeremy Roberts
 * \brief  Python interface for callow Vector
 */
//---------------------------------------------------------------------------//

%include "detran_utilities.i"

%include "Vector.hh"

namespace callow
{

//class Vector
//{
//public:
//  typedef detran_utilities::SP<Vector<T> >    SP_vector;
//  Vector();
//  Vector(const int n, T v = 0);
//  Vector(const Vector &x);
//  Vector(Vector &x);
//  Vector(std::vector<T> &x);
//  static SP_vector Create(const int n, const T v = 0.0);
//  virtual ~Vector();
//  void resize(const int n, const T v = 0.0);
//  const T& operator[](const int i) const;
//  T& operator[](const int i);
//  const T& operator()(const int i) const;
//  T& operator()(const int i);
//  const T& value(const int i) const;
//  T& value(const int i);
//  /// Inner product of this vector with vector x
//  T dot(const Vector<T>& x);
//  T dot(SP_vector x);
//  T norm(const int type = L2);
//  T norm_residual(const Vector<T>& x, const int type = L2);
//  T norm_residual(SP_vector x, const int type = L2);
//  void set(const T v);
//  void scale(const T v);
//  void add(const Vector<T>& x);
//  void add(SP_vector x);
//  void subtract(const Vector<T>& x);
//  void subtract(SP_vector x);
//  void multiply(const Vector<T>& x);
//  void multiply(SP_vector x);
//  void divide(const Vector<T>& x);
//  void divide(SP_vector x);
//  void copy(const Vector<T>& x);
//  void copy(SP_vector x);
//  void add_a_times_x(const T a, const Vector<T>& x);
//  void add_a_times_x(const T a, SP_vector x);
//  int size() const { return d_size; }
//  void display() const;
//  void print_matlab(std::string filename="vector.out") const;
//};

%extend Vector
{
   double  __getitem__(int i) 
   { 
     return (*self)[i]; 
   }
   void __setitem__(int i, double v) 
   { 
     (*self)[i] = v; 
   }
}

}

%template(VectorSP) detran_utilities::SP<callow::Vector>;

