//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   BoundaryMOC.i.hh
 * \brief  BoundaryMOC inline member definitions.
 * \author Jeremy Roberts
 * \date   Jun 26, 2012
 */
//---------------------------------------------------------------------------//

#ifndef BOUNDARYMOC_I_HH_
#define BOUNDARYMOC_I_HH_


namespace detran
{

//---------------------------------------------------------------------------//
// Inherited Interface
//---------------------------------------------------------------------------//

template <class D>
inline void BoundaryMOC<D>::set(int g)
{
  for(int side = 0; side < 2*D::dimension; side++)
  {
    d_bc[side]->set(g);
  }
}

template <class D>
inline void BoundaryMOC<D>::update(int g)
{
  for(int side = 0; side < 2*D::dimension; side++)
  {
    d_bc[side]->set(g);
  }
}

template <class D>
inline void BoundaryMOC<D>::update(int g, int o, int a)
{
  for(int side = 0; side < 2*D::dimension; side++)
  {
    d_bc[side]->update(g, o, a);
  }
}

template <class D>
inline void BoundaryMOC<D>::clear(int g)
{
  /* ... */
}

//---------------------------------------------------------------------------//
// Explicit Instantiations
//---------------------------------------------------------------------------//

// only 2d for now
template class BoundaryMOC<_2D>;


} // end namespace detran

#endif // BOUNDARYMOC_I_HH_ 

//---------------------------------------------------------------------------//
//              end of file BoundaryMOC.i.hh
//---------------------------------------------------------------------------//
