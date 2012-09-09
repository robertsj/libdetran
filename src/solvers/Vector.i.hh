//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Vector.i.hh
 * \brief  Vector inline member definitions
 * \author Jeremy Roberts
 * \date   Sep 8, 2012
 */
//---------------------------------------------------------------------------//

#ifndef detran_VECTOR_I_HH_
#define detran_VECTOR_I_HH_

namespace detran
{

// Value Setting

inline void Vector::insert_values(const size_t number,
                                  const int *rows,
                                  const double *values)
{
  VecSetValues(d_V, number, rows, values, INSERT_VALUES);
  d_is_assembled = false;
}

//---------------------------------------------------------------------------//
// Vector Operations
//---------------------------------------------------------------------------//

inline double Vector::dot(Vector &x)
{
  // Preconditions
  Require(x.is_assembled());
  Require(is_assembled());
  Require(x.size() == d_size);

  PetscErrorCode ierr;
  double val;
  ierr = VecDot(d_V, x.V(), &val);

  // Postconditions
  Ensure(!ierr);
  return val;
}

inline void Vector::scale(const double factor)
{
  // Preconditions
  Require(is_assembled());

  PetscErrorCode ierr;
  ierr = VecScale(d_V, factor);

  // Postconditions
  Ensure(!ierr);
}

//---------------------------------------------------------------------------//
// Accessors
//---------------------------------------------------------------------------//

inline const double&
Vector::operator[](const size_t i) const
{
  Require(i < d_size);
  return d_array[i];
}

inline double&
Vector::operator[](const size_t i)
{
  // Cast away return type
  return const_cast<double&>
  (
    // Add const to *this's type and call const version
    static_cast<const Vector&>(*this)[i]
  );
}

} // end namespace detran

#endif // detran_VECTOR_I_HH_

//---------------------------------------------------------------------------//
//              end of file Vector.i.hh
//---------------------------------------------------------------------------//
