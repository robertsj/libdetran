//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  MatrixDense.i.hh
 *  @brief MatrixDense inline member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef callow_MATRIXDENSE_I_HH_
#define callow_MATRIXDENSE_I_HH_

namespace callow
{

#ifdef MATRIXDENSE_COLMAJ
#define MDIDX(i, j) i + j * d_m
#else
#define MDIDX(i, j) j + i * d_n
#endif

//----------------------------------------------------------------------------//
// ACCESS
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
inline const double& MatrixDense::operator()(const int i, const int j) const
{
  Require(d_is_ready);
  Requirev(i >= 0 && i < d_m,  "0 <= " + AsString(i) + " < " + AsString(d_m));
  Requirev(j >= 0 && j < d_n,  "0 <= " + AsString(j) + " < " + AsString(d_n));
  return d_values[MDIDX(i, j)];
}
//----------------------------------------------------------------------------//
inline double& MatrixDense::operator()(const int i, const int j)
{
  Require(d_is_ready);
  Requirev(i >= 0 && i < d_m,  "0 <= " + AsString(i) + " < " + AsString(d_m));
  Requirev(j >= 0 && j < d_n,  "0 <= " + AsString(j) + " < " + AsString(d_n));
  return d_values[MDIDX(i, j)];
}

//----------------------------------------------------------------------------//
inline const double& MatrixDense::operator[](const int p) const
{
  Require(d_is_ready);
  Require(p >= 0 && p < d_n * d_m);
  return d_values[p];
}
//----------------------------------------------------------------------------//
inline double& MatrixDense::operator[](const int p)
{
  Require(d_is_ready);
  Require(p >= 0 && p < d_n * d_m);
  return d_values[p];
}

//----------------------------------------------------------------------------//
// MULTIPLY
//----------------------------------------------------------------------------//

// These and other routines could be substituted via BLAS

//----------------------------------------------------------------------------//
inline void MatrixDense::multiply(const Vector &v_in, Vector &v_out)
{
  Require(d_is_ready);
  Require(v_in.size() == d_n);
  Require(v_out.size() == d_m);
#ifdef CALLOW_ENABLE_PETSC_OPS
  MatMult(d_petsc_matrix,
         const_cast<Vector* >(&v_in)->petsc_vector(),
         y.petsc_vector());
#else
  // clear the output vector
  v_out.scale(0);
  // for all columns
  int p = 0;
  for (int j = 0; j < d_n; ++j)
  {
    double temp = v_in[j];
    // for all rows
    for (int i = 0; i < d_m; ++i, ++p)
      v_out[i] += temp * d_values[p];
  }
#endif
}

//----------------------------------------------------------------------------//
inline void MatrixDense::multiply_transpose(const Vector &v_in, Vector &v_out)
{
  Require(d_is_ready);
  Require(v_in.size() == d_m);
  Require(v_out.size() == d_n);

#ifdef CALLOW_ENABLE_PETSC_OPS
  MatMultTranspose(d_petsc_matrix,
                   const_cast<Vector*>(&v_in)->petsc_vector(),
                   y.petsc_vector());
#else
  // clear the output vector
  v_out.scale(0);
  // for all rows (now columns)
  int p = 0;
  for (int j = 0; j < d_n; ++j)
  {
    // for all columns (now rows)
    double temp = 0.0;
    for (int i = 0; i < d_m; ++i, ++p)
      temp += v_in[i] * d_values[p];
    v_out[j] = temp;
  }
#endif
}

//----------------------------------------------------------------------------//
// INSERTING VALUES
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
inline bool MatrixDense::insert(int i, int j, double v, const int type)
{
  // Preconditions
  Require(d_is_ready);
  Require(i >= 0 && i < d_m);
  Require(j >= 0 && j < d_n);

  if (type == ADD)
    d_values[MDIDX(i, j)] += v;
  else
    d_values[MDIDX(i, j)] = v;

  return true;
}

//----------------------------------------------------------------------------//
inline bool MatrixDense::insert_row(int i, double *v, const int type)
{
  Require(d_is_ready);
  Require(i >= 0 && i < d_m);

  for (int j = 0; j < d_n; ++j)
    insert(i, j, v[j], type);

  return true;
}

//----------------------------------------------------------------------------//
inline bool MatrixDense::insert_col(int j, double *v, const int type)
{
  Require(d_is_ready);
  Require(j >= 0 && j < d_n);

  for (int i = 0; i < d_n; ++i)
    insert(i, j, v[i], type);

  return true;
}

} // end namespace callow

#endif // callow_MATRIXDENSE_I_HH_

//----------------------------------------------------------------------------//
//              end of file MatrixDense.i.hh
//----------------------------------------------------------------------------//
