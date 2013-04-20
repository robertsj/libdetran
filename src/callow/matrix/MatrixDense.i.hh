//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   MatrixDense.i.hh
 *  @brief  MatrixDense.i
 *  @author Jeremy Roberts
 *  @date   Jan 7, 2013
 */
//---------------------------------------------------------------------------//

#ifndef callow_MATRIXDENSE_I_HH_
#define callow_MATRIXDENSE_I_HH_

namespace callow
{

//---------------------------------------------------------------------------//
// ACCESS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
inline const double& MatrixDense::operator()(const int i, const int j) const
{
  // Preconditions
  Require(d_is_ready);
  Require(i >= 0 and i < d_m);
  Require(j >= 0 and j < d_n);

  return d_values[j + i * d_n];
}
//---------------------------------------------------------------------------//
inline double& MatrixDense::operator()(const int i, const int j)
{
  // Preconditions
  Require(d_is_ready);
  Require(i >= 0 and i < d_m);
  Require(j >= 0 and j < d_n);

  return d_values[j + i * d_n];
}

//---------------------------------------------------------------------------//
inline const double& MatrixDense::operator[](const int p) const
{
  // Preconditions
  Require(d_is_ready);
  Require(p >= 0 and p < d_n * d_m);

  return d_values[p];
}
//---------------------------------------------------------------------------//
inline double& MatrixDense::operator[](const int p)
{
  // Preconditions
  Require(d_is_ready);
  Require(p >= 0 and p < d_n * d_m);

  return d_values[p];
}

//---------------------------------------------------------------------------//
// MULTIPLY
//---------------------------------------------------------------------------//

// These and other routines could be substituted via BLAS

//---------------------------------------------------------------------------//
inline void MatrixDense::multiply(const Vector &x, Vector &y)
{
  // Preconditions
  Require(d_is_ready);
  Require(x.size() == d_n);
  Require(y.size() == d_m);
#ifdef CALLOW_ENABLE_PETSC_OPS
  MatMultTranspose(d_petsc_matrix,
                   const_cast<Vector* >(&x)->petsc_vector(),
                   y.petsc_vector());
#else
  // clear the output vector
  y.scale(0);
  // for all rows
  int p = 0;
  for (int i = 0; i < d_m; ++i)
  {
    double temp = y[i];
    // for all columns
    for (int j = 0; j < d_n; ++j, ++p)
      temp += x[j] * d_values[p];
    y[i] = temp;
  }
#endif
}

//---------------------------------------------------------------------------//
inline void MatrixDense::multiply_transpose(const Vector &x, Vector &y)
{
  // Preconditions [m, n] -> [n, m] * x m
  Require(d_is_ready);
  Require(x.size() == d_m);
  Require(y.size() == d_n);

#ifdef CALLOW_ENABLE_PETSC_OPS
  MatMult(d_petsc_matrix,
          const_cast<Vector* >(&x)->petsc_vector(),
          y.petsc_vector());
#else
  // clear the output vector
  y.scale(0);
  // for all rows (now columns)
  int p = 0;
  for (int i = 0; i < d_m; ++i)
  {
    // for all columns (now rows)
    for (int j = 0; j < d_n; ++j, ++p)
      y[j] += x[i] * d_values[p];
  }
#endif
}

//---------------------------------------------------------------------------//
// INSERTING VALUES
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
inline bool MatrixDense::insert(int i, int j, double v, const int type)
{
  // Preconditions
  Require(d_is_ready);
  Require(i >= 0 and i < d_m);
  Require(j >= 0 and j < d_n);

  if (type == ADD)
    d_values[j + i * d_n] += v;
  else
    d_values[j + i * d_n] = v;

  return true;
}

//---------------------------------------------------------------------------//
inline bool MatrixDense::insert_row(int i, double *v, const int type)
{
  // Preconditions
  Require(d_is_ready);
  Require(i >= 0 and i < d_m);

  for (int j = 0; j < d_n; ++j)
    insert(i, j, v[j], type);

  return true;
}

//---------------------------------------------------------------------------//
inline bool MatrixDense::insert_col(int j, double *v, const int type)
{
  // Preconditions
  Require(d_is_ready);
  Require(j >= 0 and j < d_n);

  for (int i = 0; i < d_n; ++i)
    insert(i, j, v[i], type);

  return true;
}

} // end namespace callow

#endif // callow_MATRIXDENSE_I_HH_

//---------------------------------------------------------------------------//
//              end of file MatrixDense.i.hh
//---------------------------------------------------------------------------//
