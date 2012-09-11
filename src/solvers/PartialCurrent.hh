//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PartialCurrent.hh
 * \brief  PartialCurrent class definition
 * \author Jeremy Roberts
 * \date   Sep 10, 2012
 */
//---------------------------------------------------------------------------//

#ifndef PARTIALCURRENT_HH_
#define PARTIALCURRENT_HH_

namespace detran
{

void compute_partial_current()
{


  for (int g = 0; g < d_number_groups; g++)
  {
  // side --> 0=-x, 1=+x, 2=-y, 3=+y, 4=-z, 5=+z
    for (int side = 0; side < d_dimension * 2; side++)
    {

      // Loop over all cells of this side
      for (int dim0 = dim0_start; dim0 < dim0_end; dim0 += dim0_delta)
      {
        for (int dim1 = dim1_start; dim1 < dim1_end; dim1 += dim1_delta)
        {

          // compute cell and row
          // get mat index
          // get i, j, k
          // compute dtilde
        }
      }

    } // end side loop
  } // end group loop

}

} // end namespace detran

#endif // PARTIALCURRENT_HH_ 

//---------------------------------------------------------------------------//
//              end of file PartialCurrent.hh
//---------------------------------------------------------------------------//
