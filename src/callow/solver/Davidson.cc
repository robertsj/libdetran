//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Davidson.cc
 *  @brief Davidson member definitions
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "Davidson.hh"
#include "callow/matrix/Matrix.hh"
#include "callow/matrix/MatrixDense.hh"
#include "callow/solver/Eispack.hh"
#include "EigenSolverCreator.hh"

namespace callow
{

//----------------------------------------------------------------------------//
Davidson::Davidson(const double    tol,
                   const int       maxit,
                   const int       subspace_size)
  : Base(tol, maxit, "davidson")
  , d_subspace_size(subspace_size)
{
  // create projected system solver
  SP_db db = detran_utilities::InputDB::Create();
  db->put<std::string>("eigen_solver_type",   "eispack");
  db->put<double>("eigen_solver_tol",         d_tolerance);
  db->put<int>("eigen_solver_maxit",          10000);
  db->put<double>("linear_solver_atol",       d_tolerance);
  db->put<double>("linear_solver_rtol",       d_tolerance);
  db->put<int>("linear_solver_monitor_level", 0);
  db->put<int>("eigen_solver_monitor_level",  0);
  d_projected_solver = EigenSolverCreator::Create(db);
}

//----------------------------------------------------------------------------//
void Davidson::set_operators(SP_matrix A,
                             SP_matrix B,
                             SP_db     db)
{
  Insist(A, "The operator A cannot be null");
  Require(A->number_rows() == A->number_columns());
  d_A = A;
  if (B)
  {
    Assert(B->number_rows() == B->number_columns());
    Assert(A->number_rows() == B->number_columns());
    d_B = B;
  }
  else
  {
    // Just use the identity matrix for the right hand operator so
    // one set of code works for standard and generalized eigenvalue problems
    int m = d_A->number_rows();
    Matrix::SP_matrix BB(new Matrix(m, m, 1));
    for (int i = 0; i < m; ++i)
      BB->insert(i, i, 1.0, BB->INSERT);
    d_B = BB;
    d_B->assemble();
  }

  // Create the residual function
  d_A_minus_ritz_times_B = new DavidsonResidual(d_A, d_B, this);

  // Create the default preconditioner
  d_P = new DavidsonDefaultP(d_A_minus_ritz_times_B, NULL, db);
}

//----------------------------------------------------------------------------//
void Davidson::set_preconditioner(SP_pc P, const int side)
{
  Require(P);
  d_P = P;
}

//----------------------------------------------------------------------------//
void Davidson::solve_impl(Vector &u, Vector &x0)
{
  // problem size
  const size_t m = u.size();

  // easy way to pass the eigenvalue guess
  d_lambda = u[0];

  // working residual
  Vector r(m, 0.0);
  Vector t(m, 0.0);

  // initialize the orthogonal basis
  std::vector<Vector> V(d_subspace_size,   Vector(m, 0.0));
  std::vector<Vector> V_a(d_subspace_size, Vector(m, 0.0));
  std::vector<Vector> V_b(d_subspace_size, Vector(m, 0.0));
  V[0].copy(x0);
  V[0].scale(1.0/V[0].norm());
  d_A->multiply(V[0], V_a[0]); // aka  Dinv(L)MF
  d_B->multiply(V[0], V_b[0]); // aka  I - Dinv(L)MS

  // projected operators
  MatrixDense::SP_matrix A;
  MatrixDense::SP_matrix B;

  // perform outer iterations
  size_t it = 1;
  size_t outers = d_maximum_iterations / d_subspace_size + 1;
  for (size_t o = 0; o < outers; ++o)
  {
    for (size_t i = 0; i < d_subspace_size; ++i, ++it)
    {
      std::cout << "it=" << it << " o=" << o << " i=" << i << std::endl;
      // construct the projected problem, Ax=eBx
      A = new MatrixDense(i+1, i+1, 0.0);
      B = new MatrixDense(i+1, i+1, 0.0);
      for (size_t row = 0; row < i+1; ++row)
      {
        for (size_t col = 0; col < i+1; ++col)
        {
          (*A)(row, col) = V[row].dot(V_a[col]);
          (*B)(row, col) = V[row].dot(V_b[col]);
        }
      }
      d_projected_solver->set_operators(A, B);
//      A->display();
//      B->display();
      Vector y1(i+1, 0.0);
      Vector y0(i+1, 1.0);
      y0.scale(1.0/y0.norm());

      // solve the projected problem
      d_projected_solver->solve(y1, y0);
      d_lambda = d_projected_solver->eigenvalue();
      //y1.display("Y");
      //V[0].display("V0");
      // compute ritz vector, V*y1
      u.scale(0.0);
      for (size_t j = 0; j < i + 1; ++j)
        u.add_a_times_x(y1[j], V[j]);
//     u.display("U");
      // update residual  u = Vy,  r = Au = AVy =
      r.scale(0.0);
      for (size_t j = 0; j < i + 1; ++j)
      {
        r.add_a_times_x( y1[j],          V_a[j]);
        r.add_a_times_x(-d_lambda*y1[j], V_b[j]);
      }
//      r.display("R");
//      d_A_minus_ritz_times_B->multiply(u, t);
//      t.display("R part 2");

      // check for convergence
      if (monitor(it, d_lambda, r.norm()))
      {
        it = 0;
        break;
      }
      // restart if necessary (saving u as is)
      t.copy(u);
      if (i == d_subspace_size - 1) break;
      // update the subspace
      d_P->apply(r, u);
      //u.display(" V = P\R");
      for (size_t j = 0; j < i+1; ++j)
        y1[j] = V[j].dot(u); // V'*u
      for (size_t j = 0; j < i+1; ++j)
        u.add_a_times_x(-y1[j], V[j]);
      u.scale(1.0/u.norm());
      //u.display(" V ");

      V[i+1].copy(u);
      d_A->multiply(V[i+1], V_a[i+1]);
      d_B->multiply(V[i+1], V_b[i+1]);
    } // end inners
    if (it == 0)
    {
      break;
    }
    // otherwise, restart with u
    V[0].copy(t);
    d_A->multiply(V[0], V_a[0]); // aka  Dinv(L)MF
    d_B->multiply(V[0], V_b[0]); // aka  I - Dinv(L)MS
    for (size_t i = 1; i < V.size(); ++i)
    {
      V[i].scale(0.0);
      V_a[i].scale(0.0);
      V_b[i].scale(0.0);
    }

  } // end outers
}

} // end namespace callow

//----------------------------------------------------------------------------//
//              end of file Davidson.cc
//----------------------------------------------------------------------------//
