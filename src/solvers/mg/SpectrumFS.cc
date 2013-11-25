//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  SpectrumFS.cc
 *  @brief SpectrumFS member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "SpectrumFS.hh"
#include "callow/matrix/MatrixDense.hh"
#include "callow/solver/GMRES.hh"

namespace detran
{

//----------------------------------------------------------------------------//
SpectrumFS::SpectrumFS(SP_input input, SP_material material, SP_mesh mesh)
  : SpectrumBase(input, material, mesh)
{

}

//----------------------------------------------------------------------------//
SpectrumFS::vec2_dbl SpectrumFS::spectrum(double keff)
{
  using namespace callow;

  size_t ng = d_material->number_groups();
  size_t nm = d_material->number_materials();

  vec2_dbl xi(ng, vec_dbl(nm, 0.0));


  MatrixBase::SP_matrix A_p(new MatrixDense(ng, ng, 0.0));
  MatrixDense &A = *(static_cast<MatrixDense*>(A_p.bp()));

  GMRES solver(1e-9, 1e-9, 100, 20);
  solver.set_operators(A_p);

  Vector::SP_vector x(new Vector(ng, 0.0));
  vec_dbl b_v(ng, 1.0);
  if (d_input->check("mgpc_spectrum_fixed_source"))
  {
    b_v = d_input->get<vec_dbl>("mgpc_spectrum_fixed_source");
    Assert(b_v.size() == ng);
  }
  else
  {
    THROW("FUCK");
  }
  Vector::SP_vector b(new Vector(b_v));
  b->display("B");

  for (size_t m = 0; m < nm; ++m)
  {
    for (size_t g = 0; g < ng; ++g)
    {
      for (size_t gp = 0; gp < ng; ++gp)
      {
        A(g, gp) = -d_material->sigma_s(m, g, gp);
        if (d_include_fission)
        {
          A(g, gp) -= d_material->chi(m, g) * d_material->sigma_f(m, gp) / keff;
        }
        if (g == gp)
        {
          A(g, gp) += d_material->sigma_t(m, g);
        }
      }
    }
    A.print_matlab("sfs.out");

    solver.solve(*b, *x);
    x->scale(1.0 / x->norm());
    x->print_matlab("SPECTRUM");
    for (size_t g = 0; g < ng; ++g)
    {
      xi[g][m] = (*x)[g];
    }
  }

  return xi;
}

} // end namespace detran

//----------------------------------------------------------------------------//
//              end of file SpectrumFS.cc
//----------------------------------------------------------------------------//
