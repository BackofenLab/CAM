/*
 * Atom charge change constraint is responsible for tracing
 * (charged) atoms whose atomic oxidation state changes
 * between educts and products during the reaction.
 * It propagates on the diagonal entries of educts and products
 * representing non-bonding electrons, which are not constant (<> 0).
 *
 * Created on: 2013
 *  	@author Feras Nahar <naharf@informatik.uni-freiburg.de>
 */

#ifndef CHARGECHANGE_H_
#define CHARGECHANGE_H_


#include "ValenceMatrix.h"
#include <gecode/int.hh>

// define a binary propagator since the constraints are binary constraints
typedef Gecode::BinaryPropagator<Gecode::Int::IntView,
								 Gecode::Int::PC_INT_DOM> BinPropagator;

class ChargeChange: public BinPropagator
{

protected:

	  //! x0 x1 variables of type IntView represent an educt atom and a product atom accordingly
	using BinPropagator::x0;
	using BinPropagator::x1;

	  //! Sparse matrices of educts and products graphs
	const ValenceMatrix& matE;
	const ValenceMatrix& matP;

	  //! Value represents charge change between educts and products
	const int chargeVal;

public:

	  //! Construction
	  //! @param home space of the current propagator
	  //! @param x0 represents an educt variable
	  //! @param x1 represents a product variable
	  //! @param sparse matrix of educts graph
	  //! @param sparse matrix of product graph
	ChargeChange( Gecode::Space & home
			   , Gecode::Int::IntView x0
			   , Gecode::Int::IntView x1
			   , const ValenceMatrix& eduMatrixS
			   , const ValenceMatrix& proMatrixS
			   , const int chargeChange);

	  //! Propagator post function of the Charge Conservation constraint
	  //! @param home space of the current propagator
	  //! @param x0 represents an educt variable
	  //! @param x1 represents a product variable
	  //! @param sparse matrix of educts graph
	  //! @param sparse matrix of product graph
	static
	Gecode::ExecStatus
	post(Gecode::Space & home
			, Gecode::Int::IntView x0
			, Gecode::Int::IntView x1
			, const ValenceMatrix& eduMatrixS
			, const ValenceMatrix& proMatrixS
			, const int chargeChange);

	  //! Copy constructor: returns a copy of the current space and updates CSP result during search
	  //! i.e. space must be capable to return a copy of itself
	  //! @param space home of the Charge Conservation propagator
	  //! @param bool share states whether a shared copy is constructed
	  //! @param its is a reference to an object of the this propagator
	ChargeChange(Gecode::Space & home, bool share, ChargeChange & charge);

	  //! The function is responsible for propagator cloning. It calls the previously mentioned copy ctor
	  //! @param home space of the current propagator
	  //! @param bool share states whether a shared copy is constructed
	  //! @return a pointer to an updated copy of the underlying propagator
	virtual
	Propagator*
	copy(Gecode::Space& home, bool share);

	  //! Cost function of the current propagator. The returned cost is controls
	  //! when the propagator will be scheduled
	  //! @param space home of the Charge Conservation propagator
	  //! @param modification event describe how the domain of x0, x1 has changed
	  //! @return cost value of the current propagator
	virtual
	Gecode::PropCost
	cost(const Gecode::Space&, const Gecode::ModEventDelta&) const;

	  //! Executes the current propagator and reports about the propagation
	  //! @param space home of the Charge Conservation propagator
	  //! @param modification event describe how the domain of x0, x1 has changed
	  //! @return execution status of the propagator; OK posting was successful or FAILED
	virtual
	Gecode::ExecStatus
	propagate(Gecode::Space& home, const Gecode::ModEventDelta&);

	  //! Destructor
	~ChargeChange();

};


#endif /* CHARGECHANGE_H_ */
