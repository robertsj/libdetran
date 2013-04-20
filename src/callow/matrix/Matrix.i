//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Matrix.i
 *  @author Jeremy Roberts
 *  @brief  Python interface for callow Matrix
 */
//---------------------------------------------------------------------------//

%include "detran_utilities.i"

%include "MatrixBase.hh"
%include "Matrix.hh"
%include "MatrixDense.hh"

%extend Matrix
{
   double  __getitem__(int i, int j) 
   { 
     return (*self)(i, j);
   }
}

%extend MatrixDense
{
   double  __getitem__(int i, int j) 
   { 
     return (*self)(i, j);
   }
}

%template(MatrixBaseSP)  detran_utilities::SP<callow::MatrixBase>;
%template(MatrixSP)      detran_utilities::SP<callow::Matrix>;
%template(MatrixDenseSP) detran_utilities::SP<callow::MatrixDense>;
