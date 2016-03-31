#include "NeighborhoodCount.h"

	NeighborhoodCount::
	NeighborhoodCount( Gecode::Space & home,
			           Gecode::ViewArray< Gecode::Int::IntView > & atomsArray,
			           const std::vector < size_t > & _neighID,
			           const std::map< std::size_t, std::size_t > & _neighCount)
			:
				// calling parent constructor and initializing the data members
				NProp( home, atomsArray ),
				neighID( _neighID ),
				neighCount( _neighCount )
	{
	}
	// ....................................................................................

	Gecode::ExecStatus
	NeighborhoodCount::
	post( Gecode::Space & home,
		  Gecode::ViewArray< Gecode::Int::IntView > & atomsArray,
		  const std::vector < size_t > & _neighID,
		  const std::map< std::size_t, std::size_t > & _neighCount)
	{
		// calling NeighborhoodCount constructor in the post function
		(void) new (home) NeighborhoodCount( home, atomsArray, _neighID, _neighCount );
		return Gecode::ES_OK;
	}
	// ....................................................................................

	// Copy Constructor
	NeighborhoodCount::
	NeighborhoodCount( Gecode::Space & home, bool share, NeighborhoodCount & neighborhood )
			:
				NProp( home, share, neighborhood ),
				neighID( neighborhood.neighID ),
				neighCount( neighborhood.neighCount )
	{
	}
	// ....................................................................................

	// cloning this propagator
	Gecode::Propagator*
	NeighborhoodCount::copy( Gecode::Space & home, bool share )
	{
		return new (home) NeighborhoodCount( home, share, *this );
	}
	// ....................................................................................

	Gecode::PropCost
	NeighborhoodCount::cost( const Gecode::Space &,
						  	 const Gecode::ModEventDelta & ) const
	{
		// assigning a linear low cost to this propagator
		return Gecode::PropCost::linear( Gecode::PropCost::LO, x.size( ) );
	}
	// ....................................................................................

	Gecode::ExecStatus
	NeighborhoodCount::propagate( Gecode::Space & home, const Gecode::ModEventDelta & )
	{
		// go through the atoms and retrieve the corresponding their edge patterns
		// decrement their occurrences-numbers i.e. mark them as seen
		if ( x.assigned() ) {
			for ( int i = 0; i < x.size(); i++ ) {
//				std::cout << "Val: " << x[i].val() << ", PattID: " << neighID[x[i].val()]
//                        << ", Count: " << neighCount[ neighID[x[i].val()] ] << std::endl;
				if ( neighCount[ neighID[x[i].val()] ] != 0)
					neighCount[ neighID[x[i].val()] ]--;
			}

			 /*!
			  * if still exists any edge pattern which is not already scanned i.e.
			  * occurrence-number greater than 0, return a fail. That means, not
			  * all required neighborhood candidates are present.
			  */
			std::map< std::size_t, std::size_t >::const_iterator iter;
			for ( iter = neighCount.begin(); iter != neighCount.end(); ++iter ) {
				if ( (*iter).second > 0 ) {
//					 std::cout <<" neighcount : "<< (*iter).second << " !!Fail!! " << std::endl;
					return Gecode::ES_FAILED;
				}
			}
		}// if
		// std::cout << "Fix" << std::endl;
		return Gecode::ES_FIX;

	}
	// ....................................................................................

	NeighborhoodCount::
	~NeighborhoodCount()
	{
	}
	// ....................................................................................

