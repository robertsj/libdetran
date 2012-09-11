//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   BoundaryTraits.hh
 * \author robertsj
 * \date   Sep 11, 2012
 * \brief  BoundaryTraits class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef BOUNDARYTRAITS_HH_
#define BOUNDARYTRAITS_HH_

namespace detran
{


/*!
 *  \brief Boundary traits to simplify type access.
 *
 *  For fine mesh discretizations, the boundary data is stored for
 *  each surface in the form of a value (1-D), 1-D vector (2-D), or
 *  2-D vector (3-D).  This would be used for SN and diffusion.
 *
 */
template <class D>
class BoundaryTraits
{
public:
  typedef detran_utilities::vec2_dbl value_type;
};
template <>
class BoundaryTraits<_2D>
{
public:
  typedef detran_utilities::vec_dbl value_type;
};
template <>
class BoundaryTraits<_1D>
{
public:
  typedef double value_type;
};

} // end namespace detran

#endif /* BOUNDARYTRAITS_HH_ */
