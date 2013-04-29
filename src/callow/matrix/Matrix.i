//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  Matrix.i
 *  @brief Python interface for callow Matrix
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

%include "detran_utilities.i"

%ignore *::operator[];

%include "MatrixBase.hh"
%include "Matrix.hh"
%include "MatrixDense.hh"


%extend callow::Matrix
{
   double  __getitem__(int i, int j) 
   { 
     return (*self)(i, j);
   }
}

%extend callow::MatrixDense
{
   double  __getitem__(int i, int j) 
   { 
     return (*self)(i, j);
   }
}

%template(MatrixBaseSP)  detran_utilities::SP<callow::MatrixBase>;
%template(MatrixSP)      detran_utilities::SP<callow::Matrix>;
%template(MatrixDenseSP) detran_utilities::SP<callow::MatrixDense>;