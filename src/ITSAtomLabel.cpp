#include "ITSAtomLabel.h"
#include "GC_StlSetRangeIterator.h"

	ITSAtomLabel::ITSAtomLabel(Gecode::Space & home,
						   Gecode::Int::IntView x0,
						   Gecode::Int::IntView x1,
						   Gecode::Int::IntView label,
						   const std::vector< int >& eLabelID,
						   const std::vector< int >& pLabelID)
			:
				// calling parent constructor and initializing the data members
				TernaryPropagator(home, x0, x1, label),
				eduAtomLabelID(eLabelID),
				proAtomLabelID(pLabelID)
	{
	}
	// ....................................................................................

	Gecode::ExecStatus
	ITSAtomLabel::post(Gecode::Space & home,
			Gecode::Int::IntView x0,
			Gecode::Int::IntView x1,
			Gecode::Int::IntView label,
			const std::vector< int > &eLabelID,
			const std::vector< int > &pLabelID)
	{
		// calling ITS atom label constructor in the post function
		(void) new (home) ITSAtomLabel(home, x0, x1, label, eLabelID, pLabelID);
	    return Gecode::ES_OK;
	}
	// ....................................................................................

	// copy constructor
	ITSAtomLabel::ITSAtomLabel(Gecode::Space& home, bool share, ITSAtomLabel& toCopy)
	    :
	    	TernaryPropagator(home, share, toCopy),
	    	eduAtomLabelID(toCopy.eduAtomLabelID),
	    	proAtomLabelID(toCopy.proAtomLabelID)
	{
	}
	// ....................................................................................

	// cloning this propagator
	Gecode::Propagator*
	ITSAtomLabel::copy(Gecode::Space& home, bool share) {
	    return new (home) ITSAtomLabel(home, share, *this);
	}
	// ....................................................................................

	Gecode::PropCost
	ITSAtomLabel::cost(const Gecode::Space&,
			const Gecode::ModEventDelta&) const
	{
		// assigning a linear low cost to this propagator
	    return Gecode::PropCost::linear(Gecode::PropCost::LO, 3);
	}
	// ....................................................................................

	Gecode::ExecStatus
	ITSAtomLabel::propagate(Gecode::Space& home, const Gecode::ModEventDelta&) {
		std::set< int > labelsE, labelsP, labelsAll;
		  // Iterators to go through the values of educt and products variables
		Gecode::Int::ViewValues<Gecode::Int::IntView> i(x0), j(x1), all(x2);

		while (i()) { // inserting the labels of educt atoms in a set
			// std::cout << " i = " << i.val() << std::endl;
			labelsE.insert(eduAtomLabelID.at(i.val())); ++i;
		}

		while (j()) { // inserting the labels of products atoms in a set
			// std::cout << " j = " << j.val() << std::endl;
			labelsP.insert(proAtomLabelID.at(j.val())); ++j;
		}

		while (all()) {
			labelsAll.insert(all.val()); ++all;
		}

		  // performing set intersection to find out the preserved atom labels
		std::set< int > labelsIntersectEP;
		std::set_intersection(labelsE.begin(), labelsE.end(),
							  labelsP.begin(), labelsP.end(),
							  inserter(labelsIntersectEP, labelsIntersectEP.end()));
		std::set< int > labelsIntersectAll;
		std::set_intersection(labelsIntersectEP.begin(), labelsIntersectEP.end(),
								labelsAll.begin(), labelsAll.end(),
								inserter(labelsIntersectAll, labelsIntersectAll.end()));

		  // propagate labels for these atoms
		GC_StlSetRangeIterator allLabels( &labelsIntersectAll );
		GECODE_ME_CHECK(x2.inter_r(home, allLabels, false));

		if ( labelsE.size() > labelsIntersectAll.size() ) {
			i.init(x0);
			std::set<int> delEdu;
			while (i()) {
				  // check for label existence in the products domain
				if (labelsIntersectAll.find(eduAtomLabelID.at(i.val()))
						== labelsIntersectAll.end())
				{
					delEdu.insert( i.val() );
				}
				++i;
			}
			GC_StlSetRangeIterator dataE( &delEdu);

			  // prune the variables to educts values which appear in the intersection
			GECODE_ME_CHECK(x0.minus_r(home, dataE, false));
		}

		if ( labelsP.size() > labelsIntersectAll.size() ) {
			j.init(x1);
			std::set<int> delPro;
			while (j()) {
				// check for label existence in the products domain
				if (labelsIntersectAll.find(proAtomLabelID.at(j.val()))
						== labelsIntersectAll.end())
				{
					delPro.insert(j.val());
				}
				++j;
			}

			GC_StlSetRangeIterator dataP( &delPro );
			  // prune the variables to products values which appear in the intersection
			GECODE_ME_CHECK(x1.minus_r(home, dataP, false));
		}

		  // no further propagation possible if only one label left
		if ( x2.assigned() ) {
			// std::cout << "Subsumed ITSAtomLabel" << x0.val() << " " << x1.val() << " " << x2.val() << std::endl;
			return home.ES_SUBSUMED(*this);
		}
		else {
			// std::cout << "Fix ITSAtomLabel" << std::endl;
			return Gecode::ES_FIX;
		}

	}
	// ....................................................................................

	ITSAtomLabel::~ITSAtomLabel(){

	}
	// ....................................................................................
