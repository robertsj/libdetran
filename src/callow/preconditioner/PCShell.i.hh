//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PCShell.i.hh
 * \author robertsj
 * \date   Sep 20, 2012
 * \brief  PCShell.i class definition.
 */
//---------------------------------------------------------------------------//

#ifndef PCSHELL_I_HH_
#define PCSHELL_I_HH_

namespace callow
{

template <class T>
PCShell<T>::PCShell(std::string name)
  : Preconditioner<T>(name)
{
  /* ... */
}

} // end namespace callow

#endif /* PCSHELL_I_HH_ */
