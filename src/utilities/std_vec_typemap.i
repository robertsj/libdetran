//-----------------------------------------------------------------------------
// Macro for defining an in typemap for a std::vector of primitives passed 
// by value.
//
// TYPE       : The primitive C/C++ type.
// TYPE_CODE  : 
// 
// This was largely inspired by one of the typemaps in Dolfin, a C++/Python
// interface for FEniCS.
//-----------------------------------------------------------------------------

%define %pyseq_to_stlvec_typemaps(TYPE, TYPE_NAME)

%typecheck(SWIG_TYPECHECK_ ## TYPE_NAME ## _ARRAY) std::vector<TYPE> INPUTVECTOR
{
  $1 = PySequence_Check($input) ? 1 : 0;
}

%typemap (in) 
         (std::vector<TYPE> INPUTVECTOR)
         (std::vector<TYPE> tmp_vec, PyObject* item, TYPE value, int i)
{

  // First ensure the argument is actually a Python sequence.
  if (!PySequence_Check($input))
    SWIG_exception(SWIG_TypeError, "expected a sequence for argument $argnum");

  // Get the sequence length.
  Py_ssize_t pyseq_length = PySequence_Size($input);

  // Make and fill the copy.
  tmp_vec.reserve(pyseq_length);
  for (i = 0; i < pyseq_length; i++)
  {
    item = PySequence_GetItem($input, i);

    if(!SWIG_IsOK(SWIG_AsVal(TYPE)(item, &value)))
    {
      SWIG_exception(SWIG_TypeError, "expected items of sequence to be of type "\
		     "\"TYPE\" in argument $argnum");
    }
    tmp_vec.push_back(value);
    Py_DECREF(item);
  }
  $1 = tmp_vec;
}
%enddef // end macro

// Instantiations
%pyseq_to_stlvec_typemaps(int,    INT32)
%pyseq_to_stlvec_typemaps(double, DOUBLE)

