/*
 * Atom label preservation constraint for the atoms that
 * are part of ITS ring. It ensures the preservation of
 * their types between the educts and the products.
 *
 * The difference between this propagator and AtomLabel
 * propagator is that it posts an extra restriction in
 * respect to atoms label which should formulate the ring.
 *
 * Created on: 2013
 *  	@author Feras Nahar <naharf@informatik.uni-freiburg.de>
 */

#ifndef ITSATOMLABEL_H_
#define ITSATOMLABEL_H_

#include <gecode/int.hh>
#include <set>

// define a binary propagator since the constraints are binary constraints
typedef Gecode::TernaryPropagator<Gecode::Int::IntView,
								 Gecode::Int::PC_INT_DOM> TernaryPropagator;


class ITSAtomLabel: public TernaryPropagator
{

protected:

	  //! x0 x1 variables of type IntView represent an educt atom and a product atom accordingly
	using TernaryPropagator::x0;
	using TernaryPropagator::x1;

	  //! x2 variable holds the label of an ITS atom
	using TernaryPropagator::x2;

	  //! vector containing label ID of each educt resp. product atom
	const std::vector< int > & eduAtomLabelID, & proAtomLabelID;

public:

	  //! Construction
	  //! @param home space of the current propagator
	  //! @param eduVar represents an educt variable
	  //! @param proVar represents a product variable
	  //! @param assigned assigned to 1 as soon as this propagator is subsumed
	  //! @param educts graph
	  //! @param product graph
	ITSAtomLabel(Gecode::Space & home,
			   Gecode::Int::IntView eduVar,
			   Gecode::Int::IntView proVar,
			   Gecode::Int::IntView assigned,
			   const std::vector< int >& eLabelID,
			   const std::vector< int >& pLabelID);

	  //! Propagator post function of the equal label constraint
	  //! @param home space of the current propagator
	  //! @param eduVar represents an educt variable
	  //! @param proVar represents a product variable
	  //! @param assigned assigned to 1 as soon as this propagator is subsumed
	  //! @param educts graph
	  //! @param product graph
	static
	Gecode::ExecStatus
	post(Gecode::Space & home,
		   Gecode::Int::IntView eduVar,
		   Gecode::Int::IntView proVar,
		   Gecode::Int::IntView assigned,
		   const std::vector< int >& eLabelID,
		   const std::vector< int >& pLabelID);

	  //! Copy constructor: returns a copy of the current space and updates CSP result during search
	  //! i.e. space must be capable to return a copy of itself
	  //! @param space home of the equal label propagator
	  //! @param bool share states whether a shared copy is constructed
	  //! @param eqlabel is a reference to an object of the this propagator
	ITSAtomLabel(Gecode::Space& home, bool share, ITSAtomLabel& eqlabel);

	  //! The function is responsible for propagator cloning. It calls the previously mentioned copy ctor
	  //! @param home space of the current propagator
	  //! @param bool share states whether a shared copy is constructed
	  //! @return a pointer to an updated copy of the underlying propagator
	virtual
	Propagator*
	copy(Gecode::Space& home, bool share);

	  //! Cost function of the current propagator. The returned cost is controls
	  //! when the propagator will be scheduled
	  //! @param space home of the equal label propagator
	  //! @param modification event describe how the domain of x0, x1 has changed
	  //! @return cost value of the current propagator
	virtual
	Gecode::PropCost
	cost(const Gecode::Space&, const Gecode::ModEventDelta&) const;

	  //! Executes the current propagator and reports about the propagation
	  //! @param space home of the current propagator
	  //! @param modification event describe how the domain of x0, x1 has changed
	  //! @return execution status of the propagator; OK posting was successful or FAILED
	virtual
	Gecode::ExecStatus
	propagate(Gecode::Space& home, const Gecode::ModEventDelta&);

	  //! Destructor
	~ITSAtomLabel();

};

#endif /* ITSATOMLABEL_H_ */
