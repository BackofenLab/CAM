#include "MinEdgeValence.h"
#include "GC_StlSetRangeIterator.h"

//#include <algorithm>

	MinEdgeValence::MinEdgeValence(Gecode::Space & home,
			   	   	   	   	 Gecode::Int::IntView from,
			   	   	   	   	 Gecode::Int::IntView to,
			   	   	   	   	 const int minEdgeValence,
			   	   	   	   	 const ValenceMatrix& valences )
			:
			   	// calling parent constructor and initializing the data members
			   	BinPropagator(home, from, to),
				minEdgeValence(minEdgeValence),
				valences(valences)
	{}
	// ....................................................................................

	Gecode::ExecStatus
	MinEdgeValence::post(Gecode::Space & home,
			Gecode::Int::IntView x0,
			Gecode::Int::IntView x1,
			const int minEdgeValence,
			const ValenceMatrix& valences )
	{
		// calling minimal edge valence constructor in the post function
		(void) new (home) MinEdgeValence(home, x0, x1, minEdgeValence, valences );
	    return Gecode::ES_OK;
	}
	// ....................................................................................

	// copy constructor
	MinEdgeValence::MinEdgeValence(Gecode::Space& home, bool share, MinEdgeValence& toCopy)
	    :
		    BinPropagator(home, share, toCopy),
			minEdgeValence(toCopy.minEdgeValence),
			valences(toCopy.valences)

	{}
	// ....................................................................................

	// cloning this propagator
	Gecode::Propagator*
	MinEdgeValence::copy(Gecode::Space& home, bool share) {
	    return new (home) MinEdgeValence(home, share, *this);
	}
	// ....................................................................................

	Gecode::PropCost
	MinEdgeValence::cost(const Gecode::Space&,
			const Gecode::ModEventDelta&) const
	{
		// assigning a quadratic high cost to this propagator
	    return Gecode::PropCost::quadratic(Gecode::PropCost::LO, 2);
	}
	// ....................................................................................

	Gecode::ExecStatus
	MinEdgeValence::propagate(Gecode::Space& home, const Gecode::ModEventDelta&) {
		// Iterators to go through the values of educt and products variables
		Gecode::Int::ViewValues<Gecode::Int::IntView> i, j;
		i.init(x0);
		std::set<int> fromSet, toSet;
		bool allEdgesValid = true;
		while (i()) {
		  j.init(x1);
		  while (j()) {
			/*
			 * edges that show minimal valence
			 */
			if ( valences.at( i.val(), j.val() ) >= minEdgeValence )
			{
				// std::cout << "holds for " << " i = " << i.val() << ", j = " << j.val() << std::endl;
				/*
				 * if holds then store the values of the connected atoms
				 * in two sets to prune their domains to these values
				 */
				fromSet.insert(i.val());
				toSet.insert(j.val());
			} else {
				allEdgesValid = false;
//				std::cout << "NOT holds for " << " i = " << i.val() << ", j = " << j.val() << std::endl;
			}
			++j;
		  }
		  ++i;
		}// end outer while

		if (!allEdgesValid) {
			GC_StlSetRangeIterator dataFrom(&fromSet);
			GC_StlSetRangeIterator dataTo(&toSet);

			// restrict both domains to supported values
			GECODE_ME_CHECK(x0.narrow_r(home, dataFrom, false));
			GECODE_ME_CHECK(x1.narrow_r(home, dataTo, false));
		}

		/*
		// test values after pruning
		Gecode::Int::ViewValues<Gecode::Int::IntView> n, m;
		n.init(x0);
		while (n()) {
			m.init(x1);
			while (m()) {
				std::cout << "vals remained: "<< "(i = " << n.val() << ", "
						  << " j = " << m.val() << ")" << std::endl;
				++m;
			}
			++n;
		}
		*/

		if (allEdgesValid || x0.assigned() || x1.assigned()) {
			// std::cout << "subsumed" << std::endl;
			return home.ES_SUBSUMED(*this); // dispose the propagator if one of its views is assigned
		} else {
			// std::cout << "fix" << std::endl;
			return Gecode::ES_FIX;
		}
	}
	// ....................................................................................

	MinEdgeValence::~MinEdgeValence() {

	}
	// ....................................................................................
