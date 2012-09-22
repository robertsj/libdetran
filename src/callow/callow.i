//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   callow.i
 * \author Jeremy Roberts
 * \brief  Python interface for callow library.
 */
//---------------------------------------------------------------------------//

%include "detran_utilities.i"
%include "callow_config.hh"

//---------------------------------------------------------------------------//
// setup for numerical arrays
//---------------------------------------------------------------------------//

%{
#define SWIG_FILE_WITH_INIT
%}
%include "numpy.i"
%init %{
  import_array();
%}
%numpy_typemaps(double, NPY_DOUBLE, int)
%numpy_typemaps(in,     NPY_INT, int)
%apply (int* IN_ARRAY1, double* IN_ARRAY2, int DIM1) 
       {(int *j, double *v, int n)}
%apply (int* IN_ARRAY1, int* IN_ARRAY2, double* IN_ARRAY3, int DIM1) 
       {(int *i, int *j, double* v, int n)}

//---------------------------------------------------------------------------//
// vector
//---------------------------------------------------------------------------//

%include "vector/Vector.hh"
%template(VectorDouble) callow::Vector<double>;
namespace callow
{
%extend Vector<double>{
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
%template(VectorDoubleSP) detran_utilities::SP<callow::Vector<double> >;

//---------------------------------------------------------------------------//
// matrix
//---------------------------------------------------------------------------//
%include "matrix/Matrix.hh"
%template(MatrixDouble)    callow::Matrix<double>;
%template(MatrixDoubleSP)  detran_utilities::SP<callow::Matrix<double> >;



