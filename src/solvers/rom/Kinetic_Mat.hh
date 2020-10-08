//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Kinetic_Mat.hh
 *  @brief Kinetic_Mat class definition.
 *  @note  Copyright(C) 2020 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "kinetics/TimeDependentMaterial.hh"
#include "utilities/InputDB.hh"
#include "callow/matrix/MatrixBase.hh"
#include "callow/matrix/MatrixDense.hh"
#include "material/Material.hh"
#include "geometry/Mesh.hh"

using namespace detran;


class Kinetic_Mat
{
  public:
  typedef callow::MatrixDense::SP_matrix		 SP_matrix; 
  typedef callow::Vector::SP_vector                      SP_vector;
  typedef detran_utilities::InputDB::SP_input            SP_input; 
  typedef detran_geometry::Mesh::SP_mesh                 SP_mesh;
  typedef TimeDependentMaterial::SP_material             SP_material;

  Kinetic_Mat(SP_input inp, SP_mesh mesh, SP_material mat, SP_matrix basis_f, SP_matrix basis_p);
  SP_matrix precursors_decay();
  SP_matrix delayed_production();
  SP_matrix precursors_production();

 private:
   /// Input
   SP_input d_input;
   /// Mesh
   SP_mesh d_mesh;
   /// Material
   SP_material d_mat;
   /// Flux basis
   SP_matrix d_basis_f;
   /// Precursors basis
   SP_matrix d_basis_p; 
   /// Number of precursors group
   int d_number_precursor_groups;
   /// Number of cells
   int num_cells;
   /// Number of energy groups
   int num_groups;
   /// Precursors rank
   int d_rp;
   /// Flux rank
   int d_rf;
};
