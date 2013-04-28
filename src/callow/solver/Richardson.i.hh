//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Richardson.i.hh
 *  @author robertsj
 *  @date   Sep 13, 2012
 *  @brief  Richardson inline member definitions
 */
//---------------------------------------------------------------------------//

#ifndef callow_RICHARDSON_I_HH_
#define callow_RICHARDSON_I_HH_

#include <cmath>
#include <cstdio>

namespace callow
{

//---------------------------------------------------------------------------//
// SOLVE
//---------------------------------------------------------------------------//

inline void Richardson::solve_impl(const Vector &b, Vector &x)
{

  // precondition the right hand side or scale via the relaxation factor
  Vector B(b);
  if (d_P && d_pc_side == Base::LEFT)
  {
	//THROW("lalal");
    d_P->apply(*const_cast<Vector*>(&b), B);
  }
  else
  {
    B.scale(d_omega);
  }

  // temporary storage and pointers for swapping
  Vector temp(x.size(), 0.0);
  Vector* x0 = &x;
  Vector* x1 = &temp;
  Vector* swap;
  Vector res(x.size(), 0.0);

  // compute initial residual w(Ax - b) and its norm
  d_A->multiply((*x0), (*x1));
  x1->scale(d_omega);
  double r = x1->norm_residual(B, L2);
  if (monitor_init(r)) return;

  // perform iterations
  for (int iteration = 1; iteration <= d_maximum_iterations; ++iteration)
  {

    //---------------------------------------------------//
    // compute X1 <-- (I - w*A) * X0 + w*b
    //---------------------------------------------------//

    // X1 <-- A * X0
    d_A->multiply((*x0), (*x1));
	// Apply preconditioning or relaxation
	if (d_P && d_pc_side == Base::LEFT)
	{
		// X1 <-- P * X1 = P * A * X1
		Vector t(*x1);
		d_P->apply(t, *x1);
	}
	else
	{
      // X1 <-- w * X1 = w * A * X0
      x1->scale(d_omega);
	}
    // X1 <-- X1 - X0 =  (A - I) * X0
    x1->subtract(*x0);
    // X1 <-- X1 - b = (A - I) * X0 - b
    x1->subtract(B);
    // X1 <-- -X1 = (I - A) * X0 + b
    x1->scale(-1);

    //---------------------------------------------------//
    // compute residual norm
    //---------------------------------------------------//

    r = x1->norm_residual(*x0, d_norm_type);

    //---------------------------------------------------//
    // swap pointers
    //---------------------------------------------------//
    swap = x0;
    x0   = x1;
    x1   = swap;

    //---------------------------------------------------//
    // check convergence
    //---------------------------------------------------//

    if (monitor(iteration, r)) break;

  }

  // copy into the solution vector if needed
  if (x0 != &x) x.copy(*x0);

  return;
}

} // end namespace callow

#endif /* callow_RICHARDSON_I_HH_ */
