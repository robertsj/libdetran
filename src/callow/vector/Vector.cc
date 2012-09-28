//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Vector.cc
 * \author robertsj
 * \date   Sep 13, 2012
 * \brief  Vector member definitions
 */
//---------------------------------------------------------------------------//

#include "Vector.hh"

namespace callow
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR & DESTRUCTOR
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
Vector::Vector()
  : d_size(0)
  , d_temporary(false)
{
  /* ... */
}

//---------------------------------------------------------------------------//
Vector::Vector(const int n, double v)
  : d_size(n)
  , d_temporary(false)
{
  // Preconditions
  Require(d_size > 0);

  // Size the vector
#ifdef CALLOW_ENABLE_PETSC
  // Create the vector and initialize values
  PetscErrorCode ierr;
  ierr = VecCreateSeq(PETSC_COMM_SELF, n, &d_petsc_vector);
  ierr = VecSet(d_petsc_vector, v);
  // Grab the underlying storage
  ierr = VecGetArray(d_petsc_vector, &d_value);
  ierr = VecRestoreArray(d_petsc_vector, PETSC_NULL);
#else
  d_value = new double[n];
  set(v);
#endif
  // Postconditions
  Ensure(d_value);
}

//---------------------------------------------------------------------------//
Vector::Vector(const Vector &x)
  : d_size(x.size())
  , d_temporary(false)
{
  // nothing to do if x has no elements
  if (!d_size) return;
#ifdef CALLOW_ENABLE_PETSC
  // create the vector and initialize values
  PetscErrorCode ierr;
  Vec x_v = const_cast<Vector* >(&x)->petsc_vector();
  ierr = VecDuplicate(x_v, &d_petsc_vector);
  ierr = VecCopy(x_v, d_petsc_vector);
  // Grab the underlying storage
  ierr = VecGetArray(d_petsc_vector, &d_value);
  ierr = VecRestoreArray(d_petsc_vector, PETSC_NULL);
#else
  d_value = new double[d_size];
  set(0.0);
  add(x);
#endif
}

//---------------------------------------------------------------------------//
Vector::Vector(Vector &x)
  : d_size(x.size())
  , d_temporary(false)
{
  // nothing to do if x has no elements
  if (!d_size) return;
#ifdef CALLOW_ENABLE_PETSC
  // create the vector and initialize values
  PetscErrorCode ierr;
  Vec x_v = x.petsc_vector();
  ierr = VecDuplicate(x_v, &d_petsc_vector);
  ierr = VecCopy(x_v, d_petsc_vector);
  // Grab the underlying storage
  ierr = VecGetArray(d_petsc_vector, &d_value);
  ierr = VecRestoreArray(d_petsc_vector, PETSC_NULL);
#else
  d_value = new double[d_size];
  set(0.0);
  add(x);
#endif
}

//---------------------------------------------------------------------------//
Vector::Vector(std::vector<double> &x)
  : d_size(0)
  , d_temporary(true)
{
  d_value = &x[0];
  d_size  = x.size();
}

//---------------------------------------------------------------------------//
#ifdef CALLOW_ENABLE_PETSC
Vector::Vector(Vec pv)
  : d_size(0)
  , d_temporary(true)
{
  PetscErrorCode ierr;
  ierr = VecGetArray(pv, &d_value);
  ierr = VecRestoreArray(pv, PETSC_NULL);
  ierr = VecGetSize(pv, &d_size);
  d_petsc_vector = pv;
  Ensure(!ierr);
}
#endif

//---------------------------------------------------------------------------//
Vector::~Vector()
{
  if (!d_size or d_temporary) return;
#ifdef CALLOW_ENABLE_PETSC
  VecDestroy(&d_petsc_vector);
#else
  delete [] d_value;
#endif
}

//---------------------------------------------------------------------------//
void Vector::resize(const int n, const double v)
{
  Insist(!d_temporary, "Cannot resize a temporary Vector!");
#ifdef CALLOW_ENABLE_PETSC
  // destroy the vector if it's built
  if (d_size)
  {
    d_value = 0;
    VecDestroy(&d_petsc_vector);
  }
  d_size = n;
  // create a new vector
  VecCreateSeq(PETSC_COMM_SELF, n, &d_petsc_vector);
  // grab the underlying storage
  PetscErrorCode ierr;
  ierr = VecGetArray(d_petsc_vector, &d_value);
  ierr = VecRestoreArray(d_petsc_vector, PETSC_NULL);
#else
  if (d_size) delete [] d_value;
  d_size = n;
  if (d_size) d_value = new double[d_size];
  set(v);
#endif
}

//---------------------------------------------------------------------------//
// IO
//---------------------------------------------------------------------------//

void Vector::display() const
{
  printf(" Vector \n");
  printf(" ---------------------------\n");
  printf("      number rows = %5i \n\n", d_size);
  if (d_size > 100)
  {
    printf("  *** vector not printed for size > 20 *** ");
    return;
  }
  for (int i = 0; i < d_size; i++)
  {
    printf(" row  %3i | %13.6e \n", i, d_value[i]);
  }
  printf("\n");
}

void Vector::print_matlab(std::string filename) const
{
  FILE * f;
  f = fopen (filename.c_str(), "w");
  for (int i = 0; i < d_size; i++)
  {
    fprintf(f, "%23.16e \n", d_value[i]);
  }
  fprintf(f, "\n");
  fclose (f);
}


} // end namespace callow


