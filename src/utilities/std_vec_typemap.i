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

/*//-----------------------------------------------------------------------------*/
/*// Macro for defining an in typemap for a const std::vector& of primitives*/
/*// The typemaps takes a NumPy array of that primitive*/
/*//*/
/*// TYPE       : The primitive type*/
/*// TYPE_UPPER : The SWIG specific name of the type used in the array type checks*/
/*//              values SWIG use: INT32 for integer, DOUBLE for double aso.*/
/*// ARG_NAME   : The name of the argument that will be maped as an 'argout' argument*/
/*// NUMPY_TYPE : The type of the NumPy array that will be returned*/
/*// TYPE_NAME  : The name of the pointer type, 'double' for 'double', 'uint' for*/
/*//              'dolfin::uint'*/
/*// DESCR      : The char descriptor of the NumPy type*/
/*//-----------------------------------------------------------------------------*/
/*%define IN_TYPEMAP_STD_VECTOR_OF_PRIMITIVES(TYPE, TYPE_UPPER, ARG_NAME, NUMPY_TYPE, TYPE_NAME)*/

/*// The typecheck*/
/*%typecheck(SWIG_TYPECHECK_ ## TYPE_UPPER ## _ARRAY)  \*/
/*const std::vector<TYPE>&  ARG_NAME*/
/*{*/
/*  $1 = PyArray_Check($input) ? 1 : 0;*/
/*}*/

/*// The typemap*/
/*%typemap(in) const std::vector<TYPE>& ARG_NAME (std::vector<TYPE> temp)*/
/*{*/
/*  // IN_TYPEMAP_STD_VECTOR_OF_PRIMITIVES(TYPE, TYPE_UPPER, ARG_NAME,*/
/*  //                                     NUMPY_TYPE, TYPE_NAME, DESCR)*/
/*  {*/
/*    if (PyArray_Check($input))*/
/*    {*/
/*      PyArrayObject *xa = reinterpret_cast<PyArrayObject*>($input);*/
/*      if ( PyArray_TYPE(xa) == NUMPY_TYPE )*/
/*      {*/
/*        const unsigned int size = PyArray_DIM(xa, 0);*/
/*        temp.resize(size);*/
/*        TYPE* array = static_cast<TYPE*>(PyArray_DATA(xa));*/
/*        if (PyArray_ISCONTIGUOUS(xa))*/
/*          std::copy(array, array + size, temp.begin());*/
/*        else*/
/*        {*/
/*          const npy_intp strides = PyArray_STRIDE(xa, 0)/sizeof(TYPE);*/
/*          for (int i=0; i<size; i++)*/
/*            temp[i] = array[i*strides];*/
/*        }*/
/*        $1 = &temp;*/
/*      }*/
/*      else*/
/*      {*/
/*        SWIG_exception(SWIG_TypeError, "numpy array of 'TYPE_NAME' expected."\*/
/*          " Make sure that the numpy array use dtype.");*/
/*      }*/
/*    }*/
/*    else*/
/*    {*/
/*      SWIG_exception(SWIG_TypeError, "numpy array of 'TYPE_NAME' expected. "\*/
/*		     "Make sure that the numpy array use dtype.");*/
/*    }*/
/*  }*/
/*}*/

/*%enddef*/

/*IN_TYPEMAP_STD_VECTOR_OF_PRIMITIVES(double, DOUBLE, NPY_DOUBLE, double)*/
/*IN_TYPEMAP_STD_VECTOR_OF_PRIMITIVES(int, INT32, NPY_UINT, int)*/

