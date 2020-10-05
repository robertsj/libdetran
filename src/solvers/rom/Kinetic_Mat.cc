//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Kinetic_Mat.cc
 *  @brief Kinetic_Mat class definition.
 *  @note  Copyright(C) 2020 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "callow/utils/Initialization.hh"
#include "callow/vector/Vector.hh"
#include "solvers/EigenvalueManager.hh"
#include "solvers/mg/DiffusionLossOperator.hh"
#include "solvers/mg/DiffusionGainOperator.hh"
#include "kinetics/KineticsMaterial.hh"
#include "callow/matrix/MatrixShell.hh"
#include "callow/vector/Vector.hh"
#include "callow/matrix/Matrix.hh"

#include "Kinetic_Mat.hh"
using namespace detran;
using namespace detran_material;
using namespace detran_geometry;
using namespace detran_utilities;
typedef callow::MatrixDense::SP_matrix             SP_matrix;


typedef callow::MatrixDense::SP_matrix		           SP_matrix;


Kinetic_Mat::Kinetic_Mat(SP_input inp, SP_mesh mesh, SP_material mat, SP_matrix basis_f, SP_matrix basis_p)
:d_input(inp),
 d_mesh(mesh),
 d_mat(mat),
 d_number_precursor_groups(mat->number_precursor_groups()),
 d_rp(basis_p->number_columns()),
 d_rf(basis_f->number_columns()),
 d_basis_f(basis_f),
 d_basis_p(basis_p)
{
 num_cells = mesh->number_cells();
 num_groups = mat->number_groups();
 // need to assert the number of rows in the basis is related to the number of cells
}

//----------------------------------------------------------------------------//

SP_matrix Kinetic_Mat::precursors_decay()
{
// precursors decay
 SP_matrix P;

 P = new callow::MatrixDense(d_rp, d_rp);
 double PP[d_rp][d_rp];
 // loop over rows
 for (int r=0; r< d_rp; r++)
 {
  // loop over columns
  for (int c=0; c<d_rp; c++)
  {
   double value = 0.0;
   for (int i=0; i < d_number_precursor_groups; i++)
   {
     double l =  d_mat->lambda(i);
     for (int j=0; j< num_cells ; j++)
     {
      value += -l*(*d_basis_p)((i*num_cells)+j, r)*(*d_basis_p)((i*num_cells+j), c);
      //PP[r][c] +=  -l*(*d_basis_p)(i*num_cells+j, r)*(*d_basis_p)(i*num_cells+j, c);
      // std::cout << (*d_basis_p)[r, i*num_cells+j] << "  " << value << "\n";
     }
     P->insert(r, c, value);
   }
   }
 }
return P;
}

//----------------------------------------------------------------------------//
SP_matrix Kinetic_Mat::delayed_production()
{
  // The reduced
  SP_matrix D;
  D = new callow::MatrixDense(d_rf, d_rp);
  for (int r=0; r<d_rf; r++)
  {
    for (int c=0; c<d_rp; c++)
    {
	  double value = 0;
      for (int n=0; n<num_cells; n++)
      {
        for (int i=0; i<d_number_precursor_groups; i++)
        {
          value += (*d_basis_f)(n, r)*d_mat->lambda(i)*(*d_basis_p)((i*num_cells+n), c);
          D->insert(r, c, value);
       }
      }
    }
  }
 return D;
}

//----------------------------------------------------------------------------//

SP_matrix Kinetic_Mat::precursors_production()
{
  int m;
  const vec_int &mat_map = d_mesh->mesh_map("MATERIAL");
  SP_matrix F;
  F = new callow::MatrixDense(d_rp, d_rf);
  for (int r=0; r<d_rp; r++)
  {
    for (int c=0; c<d_rf; c++)
    {
      double value = 0.0;
      for (int g=0; g<num_groups; g++)
      {
        for (int p=0; p<d_number_precursor_groups; p++)
        {
          for (int n=0; n<num_cells; n++)
          {
            m = mat_map[n];
            value = (*d_basis_p)(p*num_cells + n, r)*d_mat->beta(m, p)*d_mat->nu_sigma_f(m, g)*(*d_basis_f)((n + g*num_cells), c);
            F->insert(r, c, value, 1);
          }
        }
      }
    }
  }
  F->print_matlab("mat3.txt");
  return F;
}
