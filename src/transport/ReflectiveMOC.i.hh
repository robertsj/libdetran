//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ReflectiveMOC.i.hh
 * \brief  ReflectiveMOC inline member definitions
 * \author Jeremy Roberts
 * \date   Jun 27, 2012
 */
//---------------------------------------------------------------------------//

#ifndef REFLECTIVEMOC_I_HH_
#define REFLECTIVEMOC_I_HH_

#include <iostream>

namespace detran
{

template <class D>
void ReflectiveMOC<D>::update(int g)
{
  using std::cout;
  using std::endl;

  u_int oo = 0;
  u_int aa = 0;
  u_int tt = 0;

  // List of (o, az, t) indices incident on my side
  const vec2_int &side_index = d_boundary.side_indices(d_side);

  int nt = side_index.size();

  for (int i = 0; i < nt; i++)
  {

    // Incident indices
    u_int o  = side_index[i][0];
    u_int az = side_index[i][1];
    u_int t  = side_index[i][2];

    for (int p = 0; p < d_quadrature->number_polar_octant(); p++)
    {
      u_int a = d_quadrature->angle(az, p);
      d_boundary.feed_from(o, a, t, oo, aa, tt);

//      cout << "*** i = " << i
//           << " o = " << o
//           << " a = " << a
//           << " t = " << t << endl
//           << " feeds from "
//           << " oo = " << oo
//           << " aa = " << aa
//           << " tt = " << tt << endl;
//
//      cout << " out = " << d_boundary(g, oo, aa, BoundaryMOC<D>::OUT, tt) << endl;
//      cout << "  in = " << d_boundary(g, oo, aa, BoundaryMOC<D>::IN, tt) << endl;
      d_boundary(g, o, a, BoundaryMOC<D>::IN, t) =
        d_boundary(g, oo, aa, BoundaryMOC<D>::OUT, tt);

    }

  }

}

template <class D>
void ReflectiveMOC<D>::update(int g, int o, int a)
{

  // Cheating for now.
  update(g);

}

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//---------------------------------------------------------------------------//

template class ReflectiveMOC<_2D>;

} // end namespace detran

#endif // REFLECTIVEMOC_I_HH_ 

//---------------------------------------------------------------------------//
//              end of file ReflectiveMOC.i.hh
//---------------------------------------------------------------------------//
