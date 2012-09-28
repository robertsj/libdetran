//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Matrix.i
 * \author Jeremy Roberts
 * \brief  Python interface for callow Matrix
 */
//---------------------------------------------------------------------------//

%include "detran_utilities.i"

%include "Matrix.hh"

namespace callow
{
//template <class T>
//class Matrix: public MatrixBase<T>
//{
//public:
//  enum insert_type
//  {
//    INSERT, ADD, END_INSERT_TYPE
//  };
//  typedef detran_utilities::SP<Matrix<T> >  SP_matrix;
//  typedef detran_utilities::SP<Vector<T> >  SP_vector;
//  Matrix();
//  Matrix(const int m, const int n);
//  Matrix(const int m, const int n, const int nnz);
//  Matrix(Matrix<T> &A);
//  virtual ~Matrix();
//  static SP_matrix Create(const int m, const int n);
//  static SP_matrix Create(const int m, const int n, const int nnz);
//  void preallocate(const int nnz_row);
//  void preallocate(int *nnz_rows);
//  bool insert(int  i, int  j, T  v, const int type = INSERT);
//  bool insert(int  i, int *j, T *v, int n, const int type = INSERT);
//  bool insert(int *i, int  j, T *v, int n, const int type = INSERT);
//  bool insert(int *i, int *j, T *v, int n, const int type = INSERT);
//  int start(const int i) const;
//  int diagonal(const int i) const;
//  int end(const int i) const ;
//  int column(const int p) const;
//  T operator[](const int p) const;
//  T operator()(const int i, const int j) const;
//  int number_nonzeros() const;
//  bool allocated() const {return d_allocated;}
//  void print_matlab(std::string filename = "matrix.out") const;
//  void assemble();
//  void multiply(const Vector<T> &x,  Vector<T> &y);
//  void multiply_transpose(const Vector<T> &x, Vector<T> &y);
//  void display() const;
//  // base
//  void set_size(const int m, const int n);
//  int number_rows() const { return d_m; }
//  int number_columns() const { return d_n; }
//  bool is_ready() const { return d_is_ready; }
////  void multiply(SP_vector x,  SP_vector y);
////  void multiply_transpose(SP_vector x, SP_vector y);
//};

%extend Matrix
{
   double  __getitem__(int i, int j) 
   { 
     return (*self)(i, j);
   }
}

}

//%template(MatrixDouble)   callow::Matrix<double>; 
%template(MatrixSP) detran_utilities::SP<callow::Matrix>;

