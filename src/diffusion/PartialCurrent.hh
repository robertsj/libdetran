//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PartialCurrent.hh
 * \author robertsj
 * \date   Sep 5, 2012
 * \brief  PartialCurrent class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef PARTIALCURRENT_HH_
#define PARTIALCURRENT_HH_

namespace detran_diffusion
{

/*!
 *  \class PartialCurrent
 *  \brief Given the flux, compute the partial current at global boundaries
 *
 */
class PartialCurrent
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran::InputDB::SP_input     SP_input;
  typedef detran::Material::SP_material SP_material;
  typedef detran::Mesh::SP_mesh         SP_mesh;

  //--------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //--------------------------------------------------------------------------//


  /*!
   *  \brief Constructor
   *  \param input
   *  \param material
   *  \param mesh
   */
  PartialCurrent(SP_input    input,
                 SP_material material,
                 SP_mesh     mesh);

  /*!
   *  \brief Compute partial current
   *  \param x  Flux vector
   *  \param b  Boundary source vector
   */
  void compute(Vec x, Vec b);

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Input
  SP_input d_input;

  /// Materials
  SP_material d_material;

  /// Mesh
  SP_mesh d_mesh;

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

};

} // end namespace detran_diffusion


#endif /* PARTIALCURRENT_HH_ */
