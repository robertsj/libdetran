//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  Vector.i
 *  @brief Python interface for callow Vector
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

%include "detran_utilities.i"

%ignore *::operator[];

%include "Vector.hh"

%extend callow::Vector
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

%template(VectorSP) detran_utilities::SP<callow::Vector>;

