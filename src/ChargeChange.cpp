#include "ChargeChange.h"
#include "GC_StlSetRangeIterator.h"

	ChargeChange::ChargeChange(Gecode::Space & home
			   	   	   	   	 , Gecode::Int::IntView x0
			   	   	   	   	 , Gecode::Int::IntView x1
			   	   	   	   	 , const ValenceMatrix& eduMatrixS
			   	   	   	   	 , const ValenceMatrix& proMatrixS
			   	   	   	   	 , const int chargeChange)
			:
			   	  // calling parent constructor and initializing the data members
			   	BinPropagator( home, x0, x1 )
				, matE( eduMatrixS )
				, matP( proMatrixS )
				, chargeVal( chargeChange )
	{

	}
	// ....................................................................................

	Gecode::ExecStatus
	ChargeChange::post(Gecode::Space & home
			, Gecode::Int::IntView x0
			, Gecode::Int::IntView x1
			, const ValenceMatrix& eduMatrixS
			, const ValenceMatrix& proMatrixS
			, const int chargeChange)
	{
		  // calling charge conservation constructor in the post function
		(void) new (home) ChargeChange(home, x0, x1, eduMatrixS, proMatrixS, chargeChange);
	    return Gecode::ES_OK;
	}
	// ....................................................................................

	  // copy constructor
	ChargeChange::ChargeChange(Gecode::Space& home, bool share, ChargeChange& charge)
	    :
		    BinPropagator(home, share, charge)
	    	, matE(charge.matE)
	    	, matP(charge.matP)
	    	, chargeVal( charge.chargeVal)

	{

	}
	// ....................................................................................

	  // cloning this propagator
	Gecode::Propagator*
	ChargeChange::copy(Gecode::Space& home, bool share) {
	    return new (home) ChargeChange(home, share, *this);
	}
	// ....................................................................................

	Gecode::PropCost
	ChargeChange::cost(const Gecode::Space&,
			const Gecode::ModEventDelta&) const
	{
		  // assigning a quadratic low cost to this propagator
		return Gecode::PropCost::quadratic(Gecode::PropCost::LO, 2);
	}
	// ....................................................................................

	Gecode::ExecStatus
	ChargeChange::propagate(Gecode::Space& home, const Gecode::ModEventDelta&) {
		  // Iterators to go through the values of educt and products variables
		Gecode::Int::ViewValues<Gecode::Int::IntView> i, j;
		i.init(x0);
		std::set<int> eSupport, pSupport;
		while (i()) {
		  j.init(x1);
		  while (j()) {
			// std::cout << "i= " << i.val() << " , j= " << j.val() << std::endl;
			/*
			 * Charge Conservation holds when the difference between diagonal values of the educts
			 * and products equals +- 1
			 */
			// if (matE.at(i.val(), i.val()) - matP.at(j.val(), j.val()) != 0)
			if ( abs(matE.at(i.val(), i.val()) - matP.at(j.val(), j.val())) == abs(chargeVal) ) {
				eSupport.insert(i.val());
				pSupport.insert(j.val());
			}
			++j;
		  }
		  ++i;
		}// end outer while

		GC_StlSetRangeIterator dataE(&eSupport);
		GC_StlSetRangeIterator dataP(&pSupport);

		  // prune the domain of educt variables to supported values
		GECODE_ME_CHECK(x0.narrow_r(home, dataE, false));
		  // prune the domain of products variables to supported values
		GECODE_ME_CHECK(x1.narrow_r(home, dataP, false));

		if (x0.assigned() || x1.assigned()) {
			 // std::cout << "subsumed charge conservation"  << x0.val() << " - " << x1.val() << std::endl;
			  // dispose the propagator if one of its views is assigned
			return home.ES_SUBSUMED(*this);
		} else {
				// std::cout << "fix charge conservation" << std::endl;
				  // neither subsumed nor failed the propagator won't be scheduled any more
				return Gecode::ES_FIX;
		  }
	}
	// ....................................................................................

	ChargeChange::~ChargeChange() {

	}
	// ....................................................................................




