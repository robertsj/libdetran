/*
 * WithinGroupAcceleration.hh
 *
 *  Created on: May 18, 2012
 *      Author: robertsj
 */

#ifndef WITHINGROUPACCELERATION_HH_
#define WITHINGROUPACCELERATION_HH_

// Detran
#include "CoarseMesh.hh"
#include "CurrentTally.hh"
#include "Material.hh"
#include "SweepSource.hh"

// Utilities
#include "DBC.hh"
#include "InputDB.hh"
#include "SP.hh"

// System
#include "petsc.h"

namespace detran
{

/*!
 *  \class WithinGroupAcceleration
 *  \brief Base class for within-group coarse mesh acceleration schemes.
 *
 *
 */
template <class D>
class WithinGroupAcceleration
{

public:

  /// \name Useful Typedefs
  // \{

  typedef SP<WithinGroupAcceleration>     SP_acceleration;
  typedef InputDB::SP_input               SP_input;
  typedef Material::SP_material           SP_material;
  typedef CoarseMesh::SP_coarsemesh       SP_coarsemesh;
  typedef CurrentTally::SP_currenttally   SP_currenttally;
  typedef typename
          SweepSource<D>::SP_sweepsource  SP_sweepsource;

  // \}

  /*!
   *  \brief Constructor
   *
   *  @param input          Input database
   *  @param material       Material database
   *  @param coarsemesh     Coarse mesh
   *  @param currenttally   Current tally
   */
  WithinGroupAcceleration(SP_input input,
                          SP_material material,
                          SP_coarsemesh coarsemesh,
                          SP_currenttally currenttally)
  :  d_input(input)
  ,  d_material(material)
  ,  d_coarsemesh(coarsemesh)
  ,  d_currenttally(currenttally)
  {
    /* ... */
  }

  /// Virtual destructor.
  virtual ~WithinGroupAcceleration()
  {
    VecDestroy(&d_x);
    VecDestroy(&d_b);
    KSPDestroy(&d_solver);
  };

  /// \name Public Interface
  /// \{

  /*!
   *  \brief Solve the low order equation to update the scalar flux moments
   *
   *  @param  group   Energy group for this solve
   *  @param  phi     Reference to the group flux to be updated.
   *  @param  source  Smart pointer to up-to-date sweep source.
   *
   */
  virtual void accelerate(u_int group,
                          State::moments_type &phi,
                          SP_sweepsource source) = 0;

  /// \}

protected:

  /// \name Protected Data
  /// \{

  /// Input database
  SP_input d_input;

  /// Material database
  SP_material d_material;

  /// Coarse mesh
  SP_coarsemesh d_coarsemesh;

  /// Current tally
  SP_currenttally d_currenttally;

  /// Unknowns
  Vec d_x;

  /// Right hand side
  Vec d_b;

  /// Linear Solver
  KSP d_solver;

  /// \}


};

} // end namespace detran


#endif /* WITHINGROUPACCELERATION_HH_ */
