//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   coarsemesh_fixture.hh
 * \brief  coarsemesh_fixture 
 * \author Jeremy Roberts
 * \date   Aug 9, 2012
 */
//---------------------------------------------------------------------------//

#ifndef COARSEMESH_FIXTURE_HH_
#define COARSEMESH_FIXTURE_HH_

// Detran
#include "CoarseMesh.hh"
#include "Mesh1D.hh"
#include "Mesh2D.hh"
#include "Mesh3D.hh"

namespace detran_test
{

detran::CoarseMesh::SP_coarsemesh coarsemesh_1d()
{
  using detran::CoarseMesh;
  using detran::Mesh1D;
  using detran::vec_dbl;
  using detran::vec_int;

  // Coarse mesh edges
  detran::vec_dbl cm(3, 0.0);
  cm[1] = 5.0;
  cm[2] = 10.0;

  // Fine mesh counts
  detran::vec_int fm(2, 5);
  fm[1] = 10;

  // Material (homogeneous)
  vec_int mt(2, 0);

  // Create the fine mesh
  Mesh1D::SP_mesh mesh(new Mesh1D(fm, cm, mt));

  // Create the coarse mesh and return
  CoarseMesh::SP_coarsemesh coarse(new CoarseMesh(mesh, 2));
  return coarse;

}

detran::CoarseMesh::SP_coarsemesh coarsemesh_2d()
{
  using detran::CoarseMesh;
  using detran::Mesh2D;
  using detran::vec_dbl;
  using detran::vec_int;

  // Coarse mesh edges
  detran::vec_dbl cm(3, 0.0);
  cm[1] = 5.0;
  cm[2] = 10.0;

  // Fine mesh counts
  detran::vec_int fm(2, 5);
  fm[1] = 10;

  // Material (homogeneous)
  vec_int mt(4, 0);

  // Create the fine mesh
  Mesh2D::SP_mesh mesh(new Mesh2D(fm, fm, cm, cm, mt));

  // Create the coarse mesh and return
  CoarseMesh::SP_coarsemesh coarse(new CoarseMesh(mesh, 2));
  return coarse;

}


detran::CoarseMesh::SP_coarsemesh coarsemesh_2d_b()
{
  using detran::CoarseMesh;
  using detran::Mesh2D;
  using detran::vec_dbl;
  using detran::vec_int;

  // Coarse mesh edges
  detran::vec_dbl cm(3, 0.0);
  cm[1] = 5.0;
  cm[2] = 5.0;

  // Fine mesh counts
  detran::vec_int fm(2, 2);

  // Material (homogeneous)
  vec_int mt(4, 0);

  // Create the fine mesh
  Mesh2D::SP_mesh mesh(new Mesh2D(fm, fm, cm, cm, mt));

  // Create the coarse mesh and return
  CoarseMesh::SP_coarsemesh coarse(new CoarseMesh(mesh, 2));
  return coarse;

}

detran::CoarseMesh::SP_coarsemesh coarsemesh_3d()
{
  using detran::CoarseMesh;
  using detran::Mesh3D;
  using detran::vec_dbl;
  using detran::vec_int;

  // Coarse mesh edges
  detran::vec_dbl cm(3, 0.0);
  cm[1] = 5.0;
  cm[2] = 10.0;

  // Fine mesh counts
  detran::vec_int fm(2, 5);
  fm[1] = 10;

  // Material (homogeneous)
  vec_int mt(8, 0);

  // Create the fine mesh
  Mesh3D::SP_mesh mesh(new Mesh3D(fm, fm, fm, cm, cm, cm, mt));

  // Create the coarse mesh and return
  CoarseMesh::SP_coarsemesh coarse(new CoarseMesh(mesh, 2));
  return coarse;

}

detran::CoarseMesh::SP_coarsemesh coarsemesh_3d_b()
{
  using detran::CoarseMesh;
  using detran::Mesh3D;
  using detran::vec_dbl;
  using detran::vec_int;

  // Coarse mesh edges
  detran::vec_dbl cm(3, 0.0);
  cm[1] = 5.0;
  cm[2] = 10.0;

  // Fine mesh counts
  detran::vec_int fm(2, 2);

  // Material (homogeneous)
  vec_int mt(8, 0);

  // Create the fine mesh
  Mesh3D::SP_mesh mesh(new Mesh3D(fm, fm, fm, cm, cm, cm, mt));

  // Create the coarse mesh and return
  CoarseMesh::SP_coarsemesh coarse(new CoarseMesh(mesh, 2));
  return coarse;

}

} // end namespace detran

#endif // COARSEMESH_FIXTURE_HH_ 

//---------------------------------------------------------------------------//
//              end of file coarsemesh_fixture.hh
//---------------------------------------------------------------------------//
