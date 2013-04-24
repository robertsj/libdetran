#include "BoundaryBase.hh"

namespace detran
{

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//---------------------------------------------------------------------------//

BOUNDARY_TEMPLATE class BOUNDARY_EXPORT BoundaryBase<_1D>;
BOUNDARY_TEMPLATE class BOUNDARY_EXPORT BoundaryBase<_2D>;
BOUNDARY_TEMPLATE class BOUNDARY_EXPORT BoundaryBase<_3D>;

} // end namespace detran