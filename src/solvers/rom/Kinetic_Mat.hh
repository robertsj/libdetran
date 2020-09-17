/*
 * kinetic_mats.hh
 *
 *  Created on: Sep 1, 2020
 *      Author: rabab
 */
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
	typedef callow::MatrixDense::SP_matrix                 SP_matrix;
	typedef callow::Vector::SP_vector                      SP_vector;
	typedef detran_utilities::InputDB::SP_input            SP_input;
	typedef detran_geometry::Mesh::SP_mesh                 SP_mesh;
	typedef TimeDependentMaterial::SP_material             SP_material;

	Kinetic_Mat(SP_input inp, SP_mesh mesh, SP_material mat, SP_matrix basis_f, SP_matrix basis_p);
	SP_matrix Mat1();
	SP_matrix Mat2();
	SP_matrix Mat3();

  private:
	 size_t d_number_precursor_groups;
	 int num_cells;
	 int num_groups;
	 int d_rp;
	 int d_rf;
	 SP_input d_input;
	 SP_mesh d_mesh;
	 SP_material d_mat;
	 SP_matrix d_basis_f;
	 SP_matrix d_basis_p;

};



