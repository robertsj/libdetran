//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  matrix_fixture.hh
 *  @brief Matrices for testing
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
#ifndef MATRIX_FIXTURE_HH_
#define MATRIX_FIXTURE_HH_

#include "matrix/Matrix.hh"
#include <iostream>
namespace callow
{

// test 1D finite difference stencil
Matrix::SP_matrix test_matrix_1(int n = 5)
{

  Matrix::SP_matrix A;
  A = new Matrix(n, n);
  A->preallocate(3);

  double l = -0.20;
  double d =  0.50;
  double u = -0.21;

  for (int row = 0; row < n; row++)
  {
    if (row == 0)
    {
      int c[]    = { 0, 1 };
      double v[] = { d, u };
      A->insert(row, c, v, 2);
    }
    else if (row == n - 1)
    {
      int c[]    = { n - 2, n - 1 };
      double v[] = {     l,     d };
      A->insert(row, c, v, 2);
    }
    else
    {
      int c[]    = { row - 1, row, row + 1 };
      double v[] = {       l,   d,       u };
      A->insert(row, c, v, 3);
    }
  }
  A->assemble();
  return A;
}

/*
 *  sample 2D neutron diffusion matrix
 *
 *  100 cm by 100 cm with reflecting conditions on left and right
 *  group:    0    1
 *  D         1.5  0.4
 *  SigmaR    0.1  0.1
 *  Sigma21   0.1  n/a
 */
Matrix::SP_matrix test_matrix_2(int n = 10)
{
  using std::cout;
  using std::endl;

  cout << " preallocating... " << endl;
  // total cells = 2 * n * n
  // mesh size   = n * n
  // num group   = 2
  Matrix::SP_matrix A;
  int size = 2 * n * n;
  A = new Matrix(size, size);
  A->preallocate(2*2+2);

  // cell width
  double h = 100.0 / n;
  // diffusion coefficient
  double D[] = {1.0 / ( 3*0.2263), 1.0 / ( 3*1.0119)};
  // removal cross section
  double sigma_r[] = {0.2263 - 0.2006,  1.0119-0.9355};
  double sigma_s21 = 0.0161;
  // boundary condition
  double albedo[] = {1, 0, 1, 0};

  cout << " constructing..." << endl;
  for (int g = 0; g < 2; g++)
  {
    // Loop over all cells.
    for (int cell = 0; cell < n * n; cell++)
    {
      // Compute row index.
      int row = cell + g * n * n;
      // Get the directional indices.
      int i = cell % n;
      int j = cell % (n * n);
      double tmp = std::floor(double(j)/double(n));
      j = int(tmp);
      //cout << " row = " << row << " i = " << i << " j = " << j << endl;
      // Direction-specific leakage coefficients.
      double jo[4] = {0.0, 0.0, 0.0, 0.0};
      // Index arrays to help determine if a cell surface is on the boundary.
      int bound[4] = {i, i, j, j};
      int nxyz[2][2] = {{0, n-1}, {0, n-1}};
      // leak --> 0=-x, 1=+x, 2=-y, 3=+y
      for (int leak = 0; leak < 4; leak++)
      {
        // Determine whether this is a left/right, bottom/top,
        // or south/north boundary.
        int xyz_idx = std::floor(leak / 2);
        // Determine the direction, e.g. left (-) vs right (+).
        int dir_idx = 1 - ((leak + 1) % 2);
        // Put i, j, k into an array.  The neighbor indices are a perturbation of this.
        int neig_idx[2] = {i, j};
        // Determine whether the neighbor is positive (+1) or negative (-1)
        // relative to the surface under consideration, and then decrement
        // the appropriate x, y, or z index.
        int shift_idx      = -2 * ((leak + 1) % 2) + 1;
        neig_idx[xyz_idx] += shift_idx;
        // Compute coupling coefficient
        double dtilde = 0.0;
        if (bound[leak] == nxyz[xyz_idx][dir_idx])
        {
          dtilde = ( 2.0 * D[g] * (1.0 - albedo[leak]) ) /
                   ( 4.0 * D[g] * (1.0 + albedo[leak]) +
                    (1.0 - albedo[leak]) * h);
        }
        else // not a boundary
        {
          // Get the neighbor data.
          int ii = neig_idx[0];
          int jj = neig_idx[1];
          int neig_cell = ii + jj * n;
          // Compute dtilde.
          dtilde = ( 2.0 * D[g] * D[g] ) / ( h * D[g] + h * D[g] );
          // Compute and set the off-diagonal matrix value.
          double val = - dtilde / h;
          int neig_row = neig_cell + g * n * n;
          A->insert(row, neig_row, val);
        }
        // Compute leakage coefficient for this cell and surface.
        jo[leak] = (double)shift_idx * dtilde;
      } // leak loop
      // Net leakage coefficient.
      double jnet = (jo[1] - jo[0]) / h + (jo[3] - jo[2]) / h ;

     // Compute and set the diagonal matrix value.
     double val = jnet + sigma_r[g];
     A->insert(row, row, val);
     // Add downscatter component.
     if (g == 1)
     {
       int col = cell;
       double val = -sigma_s21;
       A->insert(row, col, val);
     }

    } // row loop

  } // group loop

  cout << " assembling..." << endl;
  A->assemble();
  cout << " done." << endl;
  return A;
}

/*
 *  sample 2D neutron diffusion matrix
 *
 *  100 cm by 100 cm with reflecting conditions on left and right
 *  group:    0    1
 *  D         1.5  0.4
 *  SigmaR    0.1  0.1
 *  Sigma21   0.1  n/a
 */
Matrix::SP_matrix test_matrix_3(int n = 10)
{
  using std::cout;
  using std::endl;

  cout << " preallocating... " << endl;
  // total cells = 2 * n * n
  // mesh size   = n * n
  // num group   = 2
  Matrix::SP_matrix A;
  int size = 2 * n * n;
  A = new Matrix(size, size, 2);
  double h = 100.0 / n;
  double sigma_f[] = {0.0067,  0.1241};
  double chi[] = {1.0, 0.0};
  for (int g = 0; g < 2; g++)
  {
    // Loop over all cells.
    for (int cell = 0; cell < n*n; cell++)
    {
      // Compute row index.
      int row = cell + g * n * n;
      // Loop through source group.
      for (int gp = 0; gp < 2; gp++)
      {
        // Compute column index.
        int col = cell + gp *  n * n;
        // Fold the fission density with the spectrum.
        double val = sigma_f[g] * chi[gp];
        // Set the value.
        bool flag = A->insert(col, row, val, A->INSERT);
        Assert(flag);
      }
    } // row loop
  } // group loop
  A->assemble();
  return A;
}


} // end namespace detran

#endif // MATRIX_FIXTURE_HH_

//---------------------------------------------------------------------------//
//              end of file matrix_fixture.hh
//---------------------------------------------------------------------------//
