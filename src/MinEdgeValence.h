/*
 * A constraint that ensures minimal edge valence. According
 * to minimal required valence, it keeps just the edges which
 * conform to this value and prunes the others.
 *
 * Used to speed up the alternating cycle condition.
 *
 * Created on: 2012
 *  	@author Feras Nahar <naharf@informatik.uni-freiburg.de>
 *  	@author Martin Mann <mmann@informatik.uni-freiburg.de>
 */

#ifndef MINEDGEVALENCE_H_
#define MINEDGEVALENCE_H_

#include <gecode/int.hh>
//#include <vector>
//#include <set>

#include "ValenceMatrix.h"

// define a binary propagator since the constraints are binary constraints
typedef Gecode::BinaryPropagator<Gecode::Int::IntView,
								 Gecode::Int::PC_INT_DOM> BinPropagator;

class MinEdgeValence: public BinPropagator
{

protected:

	  //! x0 x1 variables of type IntView represent an source and target of an edge to constrain
	using BinPropagator::x0;
	using BinPropagator::x1;

	  //! the minimal valences allowed for connecting edges
	const int minEdgeValence;

	  //! edge valence information of graph
	const ValenceMatrix & valences;


public:

	  //! Construction
	  //! @param home home space of the current propagator
	  //! @param x0 represents an educt variable
	  //! @param x1 represents a product variable
	  //! @param minEdgeValence the minimal valence for allowed edges
	  //! @param valences the edge valence information of graph
	MinEdgeValence(Gecode::Space & home,
			   Gecode::Int::IntView x0,
			   Gecode::Int::IntView x1,
 	   	   	   const int minEdgeValence,
			   const ValenceMatrix& valences );

	  //! Propagator post function of the edge degree constraint
	  //! @param home home space of the current propagator
	  //! @param x0 represents an educt variable
	  //! @param x1 represents a product variable
	  //! @param minEdgeValence the minimal valence for allowed edges
	  //! @param valences the edge valence information of graph
	static
	Gecode::ExecStatus
	post(	Gecode::Space & home
			, Gecode::Int::IntView x0
			, Gecode::Int::IntView x1
		    , const int minEdgeValence
			, const ValenceMatrix& valences );


	  //! Copy constructor: returns a copy of the current space and updates CSP result during search
	  //! i.e. space must be capable to return a copy of itself
	  //! @param space home home of the edge degree propagator
	  //! @param share bool share states whether a shared copy is constructed
	  //! @param its is a reference to an object of the this propagator
	MinEdgeValence(Gecode::Space& home, bool share, MinEdgeValence& edgeDeg);

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
	~MinEdgeValence();

};


#endif /* MINEDGEVALENCE_H_ */
