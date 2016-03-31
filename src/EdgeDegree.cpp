#include "EdgeDegree.h"
#include "GC_StlSetRangeIterator.h"

	EdgeDegree::EdgeDegree(Gecode::Space& home,
			   	   	   	   	 Gecode::Int::IntView x0,
			   	   	   	   	 Gecode::Int::IntView x1,
			   	   	   	   	 const ReactionGraph::EdgeList& eduDegrees,
			   	   	   	   	 const ReactionGraph::EdgeList& proDegrees,
			   	   	   	   	 const size_t maxDegreeDiff )
			:
			   	  // calling parent constructor and initializing the data members
			   	BinPropagator(home, x0, x1),
				eEdgeList(eduDegrees),
				pEdgeList(proDegrees),
				maxDegreeDiff(maxDegreeDiff)
	{

	}
	// ....................................................................................

	Gecode::ExecStatus
	EdgeDegree::post(Gecode::Space& home,
			Gecode::Int::IntView x0,
			Gecode::Int::IntView x1,
			const ReactionGraph::EdgeList& eduDegrees,
			const ReactionGraph::EdgeList& proDegrees,
			const size_t maxDegreeDiff)
	{
		  // calling edge degree constructor in the post function
		(void) new (home) EdgeDegree(home, x0, x1, eduDegrees, proDegrees, maxDegreeDiff);
	    return Gecode::ES_OK;
	}
	// ....................................................................................

	  // copy constructor
	EdgeDegree::EdgeDegree(Gecode::Space& home, bool share, EdgeDegree& toCopy)
	    :
		    BinPropagator(home, share, toCopy),
			eEdgeList(toCopy.eEdgeList),
			pEdgeList(toCopy.pEdgeList),
			maxDegreeDiff(toCopy.maxDegreeDiff)

	{}
	// ....................................................................................

	  // cloning this propagator
	Gecode::Propagator*
	EdgeDegree::copy(Gecode::Space& home, bool share) {
	    return new (home) EdgeDegree(home, share, *this);
	}
	// ....................................................................................

	Gecode::PropCost
	EdgeDegree::cost(const Gecode::Space&,
			const Gecode::ModEventDelta&) const
	{
		  // assigning a quadratic high cost to this propagator
	    return Gecode::PropCost::quadratic(Gecode::PropCost::HI, 2);
	}
	// ....................................................................................

	Gecode::ExecStatus
	EdgeDegree::propagate(Gecode::Space& home, const Gecode::ModEventDelta&) {
		  // Iterators to go through the values of educt and products variables
		Gecode::Int::ViewValues<Gecode::Int::IntView> i, j;
		i.init(x0);
		std::set<int> eduSet, proSet, observedDiffs;
		while (i()) {
		  j.init(x1);
		  while (j()) {

			  // each atom can loose or gain at most one edge during the reaction
			size_t diff = (size_t)abs(eEdgeList.at(i.val()).size() - pEdgeList.at(j.val()).size());
			  // store observed difference to check whether not another pruning is possible
			observedDiffs.insert(diff);
			if ( diff <= maxDegreeDiff) {
				// std::cout << "holds for " << " i = " << i.val() << ", j = " << j.val() << std::endl;
				/*
				 * if holds then store the values of the educts and products
				 * in two set to prune their domains to these values
				 */
				eduSet.insert(i.val());
				proSet.insert(j.val());
			} /*else {
				std::cout << "NOT holds for " << " i = " << i.val() << ", j = " << j.val() << std::endl;
			  }*/
			++j;
		  }
		  ++i;
		}// end outer while

		GC_StlSetRangeIterator dataE(&eduSet);
		GC_StlSetRangeIterator dataP(&proSet);

		  // restrict both domains to supported values
		GECODE_ME_CHECK(x0.narrow_r(home, dataE, false));
		GECODE_ME_CHECK(x1.narrow_r(home, dataP, false));


		  // no further propagation possible if
		  // - only diff values <= maxDegreeDiff were observed (last entry <= maxDegreeDiff)
		  // OR - x0/x1 is assigned
		if ((!observedDiffs.empty() && (size_t)*(observedDiffs.rbegin()) <= maxDegreeDiff )
				|| x0.assigned() || x1.assigned())
		{
			// std::cout << "subsumed edge degree " << std::endl;
			return home.ES_SUBSUMED(*this); // dispose the propagator if one of its views is assigned
		} else {
			// std::cout << "fix edge degree" << std::endl;
			return Gecode::ES_FIX;
		}
	}

	// ....................................................................................

	EdgeDegree::~EdgeDegree() {

	}
	// ....................................................................................
