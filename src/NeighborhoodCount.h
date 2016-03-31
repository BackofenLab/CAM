/*
 * This Propagator is responsible for finding the neighborhood sets
 * which must appear (already calculated) in the ITS ring.
 *
 * The constraint that this propagator defines is used to accelerate
 * the whole mapping process since it raises immediately a failure
 * if the required neighborhood patterns don't show up in the
 * imaginary subgraph.
 *
 * Created on: 2013
 *  	@author Feras Nahar <naharf@informatik.uni-freiburg.de>
 */

#ifndef NEIGHBORHOODCOUNT_H_
#define NEIGHBORHOODCOUNT_H_

#include <gecode/int.hh>
#include <map>

// Define a  NaryPropagator since the constraints could be more that three.
typedef Gecode::NaryPropagator< Gecode::Int::IntView,
								Gecode::Int::PC_INT_VAL > NProp;

class NeighborhoodCount: public NProp {

protected:

	  //! x array of type IntView represents a list of atoms
	using NProp::x;

	  //! vector holds IDs of the neighborhood sets
	const std::vector< std::size_t > neighID;

	  //! (neighborhood ID, number of appearances) map
	std::map< std::size_t, std::size_t > neighCount;

public:

	  //! Construction
	  //! @param home space of the current propagator
	  //! @param array of int describes the atoms
	  //! @param _neighID vector of neighborhoods identifiers
	  //! @param _neighCount maps between those IDs and their numbers of occurrences
	NeighborhoodCount( Gecode::Space & home,
			    	   Gecode::ViewArray< Gecode::Int::IntView > & atomsArray,
			           const std::vector < size_t > & _neighID,
			           const std::map< std::size_t, std::size_t > & _neighCount);

	  //! Propagator post function of the NeighborhoodCount constraint
	  //! @param home space of the current propagator
	  //! @param array of int describes the atoms
	  //! @param _neighID vector of neighborhoods identifiers
	  //! @param _neighCount maps between those IDs and their numbers of occurrences
	static
	Gecode::ExecStatus
	post( Gecode::Space & home,
		  Gecode::ViewArray< Gecode::Int::IntView > & atomsArray,
		  const std::vector < size_t > & _neighID,
		  const std::map< std::size_t, std::size_t > & _neighCount);


	  //! Copy constructor: returns a copy of the current space and updates CSP result during search
	  //! i.e. space must be capable to return a copy of itself
	  //! @param space home of the NeighborhoodCount propagator
	  //! @param bool share states whether a shared copy is constructed
	  //! @param its is a reference to an object of the this propagator
	NeighborhoodCount( Gecode::Space & home, bool share, NeighborhoodCount & neighborhood );

	  //! The function is responsible for propagator cloning. It calls the previously mentioned copy ctor
	  //! @param home space of the current propagator
	  //! @param bool share states whether a shared copy is constructed
	  //! @return a pointer to an updated copy of the underlying propagator
	virtual
	Propagator*
	copy( Gecode::Space & home, bool share );

	  //! Cost function of the current propagator. The returned cost is controls
	  //! when the propagator will be scheduled
	  //! @param space home of the propagator
	  //! @param modification event describe how the domain of x has changed
	  //! @return cost value of the current propagator
	virtual
	Gecode::PropCost
	cost( const Gecode::Space &, const Gecode::ModEventDelta & ) const;

	  //! Executes the current propagator and reports about the propagation
	  //! @param space home of the propagator
	  //! @param modification event describe how the domain of x has changed
	  //! @return execution status of the propagator; OK posting was successful or FAILED
	virtual
	Gecode::ExecStatus
	propagate( Gecode::Space & home, const Gecode::ModEventDelta & );

	  //! Destructor
	~NeighborhoodCount();

};

#endif /* NEIGHBORHOODCOUNT_H_ */
