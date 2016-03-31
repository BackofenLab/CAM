/*
 * Constraint Cclass to guarantee the condition that the ITS structure should be an alternating cycle.
 * An simple cyclic alternating ITS has a ring form in which each even ring pairs corresponds to
 * bond formulation (value = 1) i.e. edges in products the were not present in educts. While odd
 * ring pairs correspond to bond breaking (value = -1) i.e. edges in educts that are removed in
 * products.
 *
 * Created on: 2012
 *  	@author Feras Nahar <naharf@informatik.uni-freiburg.de>
 */

#ifndef ALTERNATINGCYCLE_H_
#define ALTERNATINGCYCLE_H_

#include "ValenceMatrix.h"
#include <gecode/int.hh>

class AlternatingCycle: public Gecode::Propagator
{

protected:

	  //! eduFirst, eduSecond represent educts ring pair in the sequence of bonds
	  //! proFirst, proSecond represent products ring pair in the sequence of bonds
	Gecode::Int::IntView eduFirst,
						 eduSecond,
						 proFirst,
						 proSecond;

	  //! Sparse matrices of educts and products graphs
	const ValenceMatrix & matE;
	const ValenceMatrix & matP;

	  //! the targeted bond valence difference to ensure
	const int valenceDifference;

public:
	  //! Construction
	  //! @param home space
	  //! @param _eduFirst first educt variable of the educts pair
	  //! @param _eduSecond second educt variable of the educts pair
	  //! @param _proFirst first product variable of the products pair
	  //! @param _proSecond second product variable of the products pair
	  //! @param eduMatrixS sparse matrix of educts graph
	  //! @param proMatrixS sparse matrix of product graph
	  //! @param valenceDifference the targeted bond valence difference to ensure
	AlternatingCycle(Gecode::Space & home
			    , Gecode::Int::IntView _eduFirst
		 	    , Gecode::Int::IntView _eduSecond
		 	    , Gecode::Int::IntView _proFirst
		 	    , Gecode::Int::IntView _proSecond
				, const ValenceMatrix& eduMatrixS
				, const ValenceMatrix& proMatrixS
				, const int valenceDifference);

	  //! Propagator disposal, canceling subscriptions made by this propagator i.e. the propagator
	  //! will not perform propagation any more.
	  //! @param home space of alternating cycle
	  //! @return the size of the disposed propagator
	virtual
	size_t
	dispose(Gecode::Space& home);

	  //! Propagator post function of the aternating cycle constraint
	  //! @param home space
	  //! @param _eduFirst first educt variable of the educts pair
	  //! @param _eduSecond second educt variable of the educts pair
	  //! @param _proFirst first product variable of the products pair
	  //! @param _proSecond second product variable of the products pair
	  //! @param eduMatrixS sparse matrix of educts graph
	  //! @param proMatrixS sparse matrix of product graph
	  //! @param valenceDifference the targeted bond valence difference to ensure
	static
	Gecode::ExecStatus
	post(Gecode::Space & home
		    , Gecode::Int::IntView _eduFirst
	 	    , Gecode::Int::IntView _eduSecond
	 	    , Gecode::Int::IntView _proFirst
	 	    , Gecode::Int::IntView _proSecond
			, const ValenceMatrix& eduMatrixS
			, const ValenceMatrix& proMatrixS
			, const int valenceDifference);

	  //! Copy constructor: returns a copy of the current space and updates CSP result during search
	  //! i.e. space must be capable to return a copy of itself
	  //! æparam space home of the alternating cycle propagator
	  //! @param bool share states whether a shared copy is constructed
	  //! @param its is a reference to an object of the this propagator
	AlternatingCycle(Gecode::Space& home, bool share, AlternatingCycle& cycle);

	  //! The function is responsible for propagator cloning. It calls the previously mentioned copy ctor
	  //! @param home space of the current propagator
	  //! @param bool share states whether a shared copy is constructed
	  //! @return a pointer to an updated copy of the underlying propagator
	virtual
	Propagator*
	copy(Gecode::Space& home, bool share);

	  //! Cost function of the current propagator. The returned cost is responsible
	  //! when the propagator will be scheduled
	  //! @param space home of the alternating cycle propagator
	  //! @param modification event describe how the domain of the variables has changed
	  //! @return cost value of the current propagator
	virtual
	Gecode::PropCost
	cost(const Gecode::Space&, const Gecode::ModEventDelta&) const;

	  //! Executes the current propagator and reports about the propagation
	  //! @param space home of the alternating cycle propagator
	  //! @param modification event describe how the domain of the variables has changed
	  //! @return execution status of the propagator; OK posting was successful or FAILED
	virtual
	Gecode::ExecStatus
	propagate(Gecode::Space& home, const Gecode::ModEventDelta&);

	  //! Destruction
	~AlternatingCycle();

};

#endif /* ALTERNATINGCYCLE_H_ */
