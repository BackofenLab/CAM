#include "ConnectedCoverage.h"
#include "GC_StlSetRangeIterator.h"

	ConnectedCoverage::
	ConnectedCoverage( Gecode::Space & home
						 , Gecode::ViewArray< Gecode::Int::IntView > & _x
						 , const sgm::Graph_Interface::CompLabel & _complabel
						 , const size_t connNum )
			:
				  // calling parent constructor and initializing the data members
				NPropagator( home, _x )
				, complabel( _complabel )
				, connectedNum( connNum )
	{
	}

/////////////////////////////////////////////////////////////////////////////////////////////////////

	Gecode::ExecStatus
	ConnectedCoverage::
	post( Gecode::Space & home
			, Gecode::ViewArray< Gecode::Int::IntView > & _x
			, const sgm::Graph_Interface::CompLabel & _complabel
			, const size_t connNum )
	{
		  // calling ConnectedCoverage constructor in the post function
		(void) new (home) ConnectedCoverage( home, _x, _complabel, connNum );
	    return Gecode::ES_OK;
	}

/////////////////////////////////////////////////////////////////////////////////////////////////////

	// copy constructor
	ConnectedCoverage::
	ConnectedCoverage( Gecode::Space & home
						, bool share
						, ConnectedCoverage & connCover )
	    :
		    NPropagator( home, share, connCover )
		    , complabel( connCover.complabel )
			, connectedNum( connCover.connectedNum )
	{

	}

/////////////////////////////////////////////////////////////////////////////////////////////////////

	  // cloning this propagator
	Gecode::Propagator*
	ConnectedCoverage::
	copy( Gecode::Space & home, bool share )
	{
	    return new (home) ConnectedCoverage( home, share, *this );
	}

/////////////////////////////////////////////////////////////////////////////////////////////////////

	Gecode::PropCost
	ConnectedCoverage::
	cost( const Gecode::Space &
		  , const Gecode::ModEventDelta & ) const
	{
		  // assigning a linear low cost to this propagator
		return Gecode::PropCost::linear( Gecode::PropCost::LO, x.size() );
	}

/////////////////////////////////////////////////////////////////////////////////////////////////////

	Gecode::ExecStatus
	ConnectedCoverage::
	propagate( Gecode::Space & home, const Gecode::ModEventDelta & )
	{
		  // store labels of connected components in a set
		std::set< int > compSet ;

		  // run over all domains and check if all connected components present
		for ( int i = 0; i < x.size(); i++ ) {
			  // get iterator on values
			Gecode::Int::ViewValues<Gecode::Int::IntView> xiVal(x[i]);
			while (xiVal()) {
				compSet.insert( complabel.at( xiVal.val() ) );
				  // fast propagation end
				if (compSet.size() == connectedNum) {
					return Gecode::ES_FIX;
				}
				  // next value;
				++xiVal;
			}
		}

		  /*
		   * if not the same that means not all connected
		   * components were covered --> FAIL
		   */
		if ( compSet.size() < connectedNum ) {
			return Gecode::ES_FAILED;
		}

		return Gecode::ES_FIX;

	}

/////////////////////////////////////////////////////////////////////////////////////////////////////

	ConnectedCoverage::~ConnectedCoverage(){

	}

/////////////////////////////////////////////////////////////////////////////////////////////////////
