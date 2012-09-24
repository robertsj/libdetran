//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   GainOperator.cc
 * \brief  GainOperator 
 * \author Jeremy Roberts
 * \date   Jul 25, 2012
 */
//---------------------------------------------------------------------------//

// Configuration
#include "config/detran_config.hh"

#ifdef DETRAN_ENABLE_PETSC

// Detran
#include "GainOperator.hh"

// System
#include <string>
#include <cmath>
#include <iostream>

namespace detran_diffusion
{

GainOperator::GainOperator(SP_input    input,
                           SP_material material,
                           SP_mesh     mesh)
  : BaseOperator(input, material, mesh)
  , d_number_groups(d_material->number_groups())
  , d_group_size(d_mesh->number_cells())
{
  // Preconditions

  using std::string;

  // The matrix dimension
  d_size = d_group_size * d_number_groups;

  // Construct the matrix.
  construct();

}

//---------------------------------------------------------------------------//
// IMPLEMENTATION
//---------------------------------------------------------------------------//

void GainOperator::construct()
{

  using std::cout;
  using std::endl;

  // Create the matrix.
  int number_nz = d_number_groups;

  //MatCreateSeqAIJ(MPI_Comm comm,PetscInt m,PetscInt n,PetscInt nz,const PetscInt nnz[],Mat *A)
  MatCreateSeqAIJ(PETSC_COMM_SELF, d_size, d_size,
                  number_nz, PETSC_NULL, &d_operator);
//  MatCreate(PETSC_COMM_SELF, &d_operator);
//  MatSetType(d_operator, MATAIJ);
//  MatSetSizes(d_operator, PETSC_DECIDE, PETSC_DECIDE, d_size, d_size);

  // Get the material map.
  vec_int mat_map = d_mesh->mesh_map("MATERIAL");

  for (int g = 0; g < d_number_groups; g++)
  {

    // Loop over all cells.
    for (int cell = 0; cell < d_group_size; cell++)
    {

      // Compute row index.
      int row = cell + g * d_group_size;

      // Define the data for this cell.
      int m = mat_map[cell];

      // Get the directional indices.
      int i = d_mesh->cell_to_i(cell);
      int j = d_mesh->cell_to_j(cell);
      int k = d_mesh->cell_to_k(cell);

      // Loop through source group.
      for (int gp = 0; gp < d_number_groups; gp++)
      {
        // Compute column index.
        int col = cell + gp * d_group_size;

        // Fold the fission density with the spectrum.
        double val = d_material->nu_sigma_f(m, gp) *
                     d_material->chi(m, g);

        // Set the value.
        MatSetValue(d_operator, row, col, val, INSERT_VALUES);
      }

    } // row loop

  } // group loop

  // Assemble.
  MatAssemblyBegin(d_operator, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(d_operator, MAT_FINAL_ASSEMBLY);

}

} // end namespace detran

#endif

//---------------------------------------------------------------------------//
//              end of file GainOperator.cc
//---------------------------------------------------------------------------//
