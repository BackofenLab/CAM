/*
 * A class represents a constraint that is responsible for edge degree condition.
 * It describes the fact the each candidate atom for an ITS should gain or loose at
 * most one edge during the reaction, due to alternating ring condition.
 *
 * Created on: 2012
 *  	@author Feras Nahar <naharf@informatik.uni-freiburg.de>
 */

#ifndef EDGEDEGREE_H_
#define EDGEDEGREE_H_

#include <gecode/int.hh>
#include "ReactionGraph.h"

// define a binary propagator since the constraints are binary constraints
typedef Gecode::BinaryPropagator<Gecode::Int::IntView,
								  Gecode::Int::PC_INT_DOM> BinPropagator;

class EdgeDegree: public BinPropagator
{

protected:

	  //! x0 x1 variables of type IntView represent an educt atom and a product atom accordingly
	using BinPropagator::x0;
	using BinPropagator::x1;

	  //! Degree information of educts and products graphs
	const ReactionGraph::EdgeList & eEdgeList;
	const ReactionGraph::EdgeList & pEdgeList;

	  //! the maximally allowed difference in degree
	const size_t maxDegreeDiff;

public:

	  //! Construction
	  //! @param home home space of the current propagator
	  //! @param x0 represents an educt variable
	  //! @param x1 represents a product variable
	  //! @param eduDegrees vector containing the number of out-edges for each atom from the educts
	  //! @param proDegrees vector containing the number of out-edges for each atom from the products
	  //! @param maxDegreeDiff the maximally allowed difference in degree
	EdgeDegree(Gecode::Space & home,
			   Gecode::Int::IntView x0,
			   Gecode::Int::IntView x1,
			   const ReactionGraph::EdgeList& eduDegrees,
			   const ReactionGraph::EdgeList& proDegrees,
			   const size_t maxDegreeDiff );

	  //! Propagator post function of the edge degree constraint
	  //! @param home home space of the current propagator
	  //! @param x0 represents an educt variable
	  //! @param x1 represents a product variable
	  //! @param eduDegrees vector containing the number of out-edges for each atom from the educts
	  //! @param proDegrees vector containing the number of out-edges for each atom from the products
	  //! @param maxDegreeDiff the maximally allowed difference in degree
	static
	Gecode::ExecStatus
	post(Gecode::Space & home, Gecode::Int::IntView x0, Gecode::Int::IntView x1,
			const ReactionGraph::EdgeList& eduDegrees,
			const ReactionGraph::EdgeList& proDegree,
			const size_t maxDegreeDiff );


	  //! Copy constructor: returns a copy of the current space and updates CSP result during search
	  //! i.e. space must be capable to return a copy of itself
	  //! @param space home home of the edge degree propagator
	  //! @param share bool share states whether a shared copy is constructed
	  //! @param edgeDeg is a reference to an object of the this propagator
	EdgeDegree(Gecode::Space& home, bool share, EdgeDegree& edgeDeg);

	  //! The function is responsible for propagator cloning. It calls the previously mentioned copy ctor
	  //! @param home home space of the current propagator
	  //! @param bool share states whether a shared copy is constructed
	  //! @return a pointer to an updated copy of the underlying propagator
	virtual
	Propagator*
	copy(Gecode::Space& home, bool share);

	  //! Cost function of the current propagator. The returned cost controls
	  //! when the propagator will be scheduled
	  //! @param space home of the edge degree propagator
	  //! @param modification event describe how the domain of x0, x1 has changed
	  //! @return cost value of the current propagator
	virtual
	Gecode::PropCost
	cost(const Gecode::Space&, const Gecode::ModEventDelta&) const;

	  //! Executes the current propagator and reports about the propagation
	  //! @param space home of the edge degree propagator
	  //! @param modification event describe how the domain of x0, x1 has changed
	  //! @return execution status of the propagator; OK posting was successful or FAILED
	virtual
	Gecode::ExecStatus
	propagate(Gecode::Space& home, const Gecode::ModEventDelta&);

	  //! Destruction
	~EdgeDegree();

};


#endif /* EDGEDEGREE_H_ */
