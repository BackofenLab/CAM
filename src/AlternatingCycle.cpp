#include "AlternatingCycle.h"
#include "GC_StlSetRangeIterator.h"

using namespace Gecode;

	AlternatingCycle::
	AlternatingCycle( Gecode::Space & home
			, Int::IntView _eduFirst
			, Int::IntView _eduSecond
			, Int::IntView _proFirst
			, Int::IntView _proSecond
			, const ValenceMatrix& eduMatrixS
			, const ValenceMatrix& proMatrixS
			, const int valenceDifference
): Propagator(home)
	, eduFirst(_eduFirst)
	, eduSecond(_eduSecond)
	, proFirst(_proFirst)
	, proSecond(_proSecond)
	, matE(eduMatrixS)
	, matP(proMatrixS)
	, valenceDifference(valenceDifference)

	{
		eduFirst.subscribe( home, *this, Int::PC_INT_DOM );
		eduSecond.subscribe( home, *this, Int::PC_INT_DOM );
		proFirst.subscribe( home, *this, Int::PC_INT_DOM );
		proSecond.subscribe( home, *this, Int::PC_INT_DOM );
	}
	// ....................................................................................

	size_t
	AlternatingCycle::
	dispose( Space& home )
	{
		eduFirst.cancel( home, *this, Int::PC_INT_DOM );
		eduSecond.cancel( home, *this, Int::PC_INT_DOM );
	    proFirst.cancel( home, *this, Int::PC_INT_DOM );
	    proSecond.cancel( home, *this, Int::PC_INT_DOM );
	    (void) Gecode::Propagator::dispose(home);
	    return sizeof(*this);
	}
	// ....................................................................................

	ExecStatus
	AlternatingCycle::
	post( Space & home
			, Int::IntView _eduFirst
			, Int::IntView _eduSecond
			, Int::IntView _proFirst
			, Int::IntView _proSecond
			, const ValenceMatrix& eduMatrixS
			, const ValenceMatrix& proMatrixS
			, const int valenceDifference )
	{
		  // calling alternating cycle constructor in the post function
		(void) new (home) AlternatingCycle(home, _eduFirst, _eduSecond,
				_proFirst, _proSecond, eduMatrixS, proMatrixS, valenceDifference);
	    return Gecode::ES_OK;
	}
	// ....................................................................................

	  // copy constructor
	AlternatingCycle::
	AlternatingCycle(Space& home
			, bool share
			, AlternatingCycle& toCopy )

	: Propagator(home, share, toCopy)
	, matE(toCopy.matE)
	, matP(toCopy.matP)
	, valenceDifference(toCopy.valenceDifference)
	{
		eduFirst.update( home, share, toCopy.eduFirst );
		eduSecond.update( home, share, toCopy.eduSecond );
		proFirst.update( home, share, toCopy.proFirst );
		proSecond.update( home, share, toCopy.proSecond );
	}
	// ....................................................................................

	Propagator*
	AlternatingCycle::copy(Space& home, bool share) {
	    return new (home) AlternatingCycle(home, share, *this);
	}
	// ....................................................................................

	PropCost
	AlternatingCycle::cost(const Space&, const ModEventDelta&) const
	{
		 /*
		  * assigning a crazy low cost to this propagator since its complexity greater
		  * than cubic (large polynomial).
		  */
	    return PropCost::crazy(PropCost::LO, 4);
	}
	// ....................................................................................

	ExecStatus
	AlternatingCycle::propagate(Space& home, const ModEventDelta&) {
		// Iterators to go through the values of educts and products variables
		Gecode::Int::ViewValues<Gecode::Int::IntView> pIterFirst, pIterSecond, eIterFirst, eIterSecond;
		std::set<int> eSetFirst, eSetSecond, pSetFirst, pSetSecond;
		pIterFirst.init(proFirst);
		while (pIterFirst()) {
		  pIterSecond.init(proSecond);
		  while (pIterSecond()) {
			  eIterFirst.init(eduFirst);
			  while (eIterFirst()) {
				eIterSecond.init(eduSecond);
			  	while (eIterSecond()) {
//			  		std::cout << "Pi = " << pIterFirst.val() << ", "
//			  				  << "Pj = " << pIterSecond.val() << " - ";
//			  		std::cout << "Ei = " << eIterFirst.val() << ", "
//			  				  << "Ej = " << eIterSecond.val() << std::endl;
					/*
					 * bond formation (val = 1) in case of an even index
					 * if the subtraction subscribes to bond formulation add the
					 * corresponding domain values of each variable to the support sets
					 */
					if (matP.at(pIterFirst.val(), pIterSecond.val()) -
							matE.at(eIterFirst.val(), eIterSecond.val()) == valenceDifference)
					{
						eSetFirst.insert(eIterFirst.val());
						eSetSecond.insert(eIterSecond.val());
						pSetFirst.insert(pIterFirst.val());
						pSetSecond.insert(pIterSecond.val());
					}
			  		++eIterSecond;
			  	}
			  	++eIterFirst;
			  }
			  ++pIterSecond;
		  }
		  ++pIterFirst;
		}// end most outer while

		GC_StlSetRangeIterator eDataFrist(&eSetFirst);
		GC_StlSetRangeIterator eDataSecond(&eSetSecond);
		GC_StlSetRangeIterator pDataFirst(&pSetFirst);
		GC_StlSetRangeIterator pDataSecond(&pSetSecond);

		// pruning domains to supported values
		GECODE_ME_CHECK( eduFirst.narrow_r(home, eDataFrist, false) );
		GECODE_ME_CHECK( eduSecond.narrow_r(home, eDataSecond, false) );
		GECODE_ME_CHECK( proFirst.narrow_r(home, pDataFirst, false) );
		GECODE_ME_CHECK( proSecond.narrow_r(home, pDataSecond, false) );

		if (eduFirst.assigned() && eduSecond.assigned() && proFirst.assigned() && proSecond.assigned()) {
			// std::cout << "subsumed alternating cycle" <<  std::endl;
			return home.ES_SUBSUMED(*this); // dispose the propagator if its all views are assigned
		} else {
		    // std::cout << "fix alternating cycle" << std::endl;
			return Gecode::ES_FIX;
		}


	}
	// ....................................................................................

	AlternatingCycle::~AlternatingCycle() {

	}
	// ....................................................................................
