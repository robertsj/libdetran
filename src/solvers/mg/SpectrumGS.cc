//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  SpectrumGS.cc
 *  @brief SpectrumGS member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "SpectrumGS.hh"
#include "callow/matrix/MatrixDense.hh"
#include "callow/solver/Eispack.hh"

namespace detran
{

//----------------------------------------------------------------------------//
SpectrumGS::SpectrumGS(SP_input input, SP_material material, SP_mesh mesh)
  : SpectrumBase(input, material, mesh)
{

}

//----------------------------------------------------------------------------//
SpectrumGS::vec2_dbl SpectrumGS::spectrum(double keff)
{
  using namespace callow;

  size_t ng = d_material->number_groups();
  size_t nm = d_material->number_materials();

  vec2_dbl xi(ng, vec_dbl(nm, 0.0));

  MatrixBase::SP_matrix U_p(new MatrixDense(ng, ng, 0.0));
  MatrixBase::SP_matrix L_p(new MatrixDense(ng, ng, 0.0));

  MatrixDense &U = *(static_cast<MatrixDense*>(U_p.bp()));
  MatrixDense &L = *(static_cast<MatrixDense*>(L_p.bp()));

  Eispack solver(1e-9, 1e3);
  solver.set_operators(U_p, L_p);

  Vector::SP_vector xi_v(new Vector(ng, 0.0));
  Vector::SP_vector xi_v_0(new Vector(ng, 0.0));

  for (size_t m = 0; m < 1; ++m)
  {
    for (size_t g = 0; g < ng; ++g)
    {
      // Upper diagonal scattering + fission
      for (size_t gp = g+1; gp < ng; ++gp)
      {
        U(g, gp) = d_material->sigma_s(m, g, gp);
        if (d_include_fission)
        {
          U(g, gp) += d_material->chi(m, g) * d_material->sigma_f(m, gp);
        }
      }
      // Lower + central diagonal with total - scatter - fission
      L(g, g) = d_material->sigma_t(m, g);
      for (size_t gp = 0; gp <= g; ++gp)
      {
        L(g, gp) -= d_material->sigma_s(m, g, gp);
        if (d_include_fission)
        {
          L(g, gp) -= d_material->chi(m, g) * d_material->sigma_f(m, gp);
        }
      }
    }
    solver.solve(xi_v, xi_v_0);
    for (size_t g = 0; g < ng; ++g)
    {
      xi[g][m] = (*xi_v)[g];
    }
  }

  return xi;
}

} // end namespace detran

//----------------------------------------------------------------------------//
//              end of file SpectrumGS.cc
//----------------------------------------------------------------------------//
