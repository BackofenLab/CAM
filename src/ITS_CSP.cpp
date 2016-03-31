
#include "ITS_CSP.h"

using namespace Gecode;


//////////////////////////////////////////////////////////////////////////////


ITS_CSP::
ITS_CSP( const ITS & its
		, const ReactionGraph & educts
		, const ReactionGraph & products
) :	its(its)
	, educts(educts)
	, products(products)
	, eduITS( *this, its.getSize(), 0, educts.getGraphToSearchSize()-1 )
	, proITS( *this, its.getSize(), 0, products.getGraphToSearchSize()-1 )
	, itsAtomLabel( *this, its.getSize(), 0, ReactionGraph::getAtomLabelNumber()-1 )

{

	//////////////   CONSTRAINT SETUP  //////////////////////

	  // ensure bijective mapping i.e. distinct assignment of atoms
	distinct( *this, eduITS );
	distinct( *this, proITS );

	const int eductsFirstProtonIndex = educts.getGraphToSearchSize() - educts.getGraphToSearchProtonNumber();
	assert(eductsFirstProtonIndex > 0);
	for ( size_t i = 0; i < its.getSize(); i++ )
	{
	      // atom label preservation between educts and products
		ensureLabelITS( *this, eduITS[i], proITS[i], itsAtomLabel[i] );

		  // ensure gaining or loosing a specific number of edges
		  // (partly depends on atom charge change)
		edgeDegree( *this, eduITS[i], proITS[i], std::max(1, abs(its.getChargeChange().at(i))) );

		  // ensure a change in atom charge
		chargeChange( *this, eduITS[i], proITS[i], its.getChargeChange().at(i) );

		  // ITS ring symmetry breaking
		if (i>0) {
			rel(*this, eduITS[0], IRT_LE, eduITS[i]);
		}

	}


	for ( size_t i = 0; i < (its.getSize()-1); i++ )
	{
		  // ensure alternating cycle structure of the ITS in the mapping
		alternateCycle( *this, eduITS[i], eduITS[i+1], proITS[i], proITS[i+1], its.getBondChange().at(i) );
	}

	  // ensure alternating cycle for ring closure
	alternateCycle( *this, eduITS[its.getSize()-1], eduITS[0], proITS[its.getSize()-1], proITS[0],
			its.getBondChange().at(its.getSize()-1) );

	/////////////////////////  LOCAL NEIGHBORHOOD PREPROCESSING  /////////////////////////

	  // Constraining ITS-participated atom to conform to the pre-computed number of occurrences
	std::map< size_t, size_t> itsLabelCount = ReactionGraph::getITS_Members( educts.getNeighborhood(), products.getNeighborhood() );
	for ( std::map< size_t, size_t>::const_iterator it = itsLabelCount.begin();
			it != itsLabelCount.end(); ++it )
	{
		Gecode::count( *this, itsAtomLabel, (*it).first, IRT_GQ, (*it).second );
	}

	  // Constraint to guarantee the required neighborhood candidates are available in the ITS
	countNeighborhood( *this, eduITS, educts.getNeighIDs(), ReactionGraph::
			getNeighCount( educts.getNeighborhood(), products.getNeighborhood() ) );
	countNeighborhood( *this, proITS, products.getNeighIDs(), ReactionGraph::
			getNeighCount( products.getNeighborhood(), educts.getNeighborhood() ) );

	/////////////////////////  SYMMETRY BREAKING  /////////////////////////

	  // branching on educts
	branch(*this, eduITS, INT_VAR_SIZE_MIN(), INT_VAL_MIN());
	  // branching on products
	branch(*this, proITS, INT_VAR_SIZE_MIN(), INT_VAL_MIN());




}

/////////////////////////////////////////////////////////////////////////

ITS_CSP::~ITS_CSP() {

}

/////////////////////////////////////////////////////////////////////////

  // copy constructor used by copy function
ITS_CSP::
ITS_CSP( bool share, ITS_CSP & itsCSP )
	:
		  // call parent class constructor
		Space( share, itsCSP )
		, its( itsCSP.its )
		, educts( itsCSP.educts )
		, products (itsCSP.products )

{
	  // modifying data members during the copying
	eduITS.update( *this, share, itsCSP.eduITS );
	proITS.update( *this, share, itsCSP.proITS );
	itsAtomLabel.update( *this, share, itsCSP.itsAtomLabel );
}

/////////////////////////////////////////////////////////////////////////

Gecode::Space*
ITS_CSP::
copy( bool share )
{
	  // call copy constructor
	return new ITS_CSP( share, *this );
}

/////////////////////////////////////////////////////////////////////////

sgm::Match
ITS_CSP::
getEductITS() const {
	  // get assignment
	sgm::Match assignment(eduITS.size(), 0);
	for (size_t i=0; i<assignment.size(); ++i) {
		assert(eduITS[i].assigned());
		assignment[i] = (size_t) eduITS[i].val();
	}
	return assignment;
}

/////////////////////////////////////////////////////////////////////////

sgm::Match
ITS_CSP::
getProductITS() const {
	  // get assignment
	sgm::Match assignment(proITS.size(), 0);
	for (size_t i=0; i<assignment.size(); ++i) {
		assert(proITS[i].assigned());
		assignment[i] = (size_t) proITS[i].val();
	}
	return assignment;
}


/////////////////////////////////////////////////////////////////////////

void
ITS_CSP::
edgeDegree( Home home
		, IntVar x0
		, IntVar x1
		, const size_t maxDegreeDiff )
{
	  // checking for failure if not then posting edge degree constraint
	if ( home.failed() )
		return;

	GECODE_ES_FAIL(EdgeDegree::post(home, x0, x1, educts.getEdgeList(),
			products.getEdgeList(), maxDegreeDiff));
}

/////////////////////////////////////////////////////////////////////////

void
ITS_CSP::
chargeChange( Home home
		, IntVar x0
		, IntVar x1
		, const int chargeChangeVal)
{
	if ( home.failed() )
		return;
	GECODE_ES_FAIL( ChargeChange::post( home, x0, x1, educts.getValenceMatrix(),
			products.getValenceMatrix(), chargeChangeVal ) );
}

/////////////////////////////////////////////////////////////////////////

void
ITS_CSP::
alternateCycle( Home home
		, IntVar x0
		, IntVar x1
		, IntVar x2
		, IntVar x3
		, const int valenceDiff )
{
	if ( home.failed() )
		return;
	GECODE_ES_FAIL(AlternatingCycle::post(home, x0, x1, x2, x3,
			educts.getValenceMatrix(), products.getValenceMatrix(), valenceDiff ));
}

/////////////////////////////////////////////////////////////////////////

void
ITS_CSP::
minEdgeValence( Home home, IntVar x0, IntVar x1
				 , const int minEdgeValence
				 , const ValenceMatrix & valences )
{
	  // checking for failure if not then posting edge degree constraint
	if ( home.failed() )
		return;
	GECODE_ES_FAIL(MinEdgeValence::post(home, x0, x1, minEdgeValence, valences));
}

/////////////////////////////////////////////////////////////////////////

void
ITS_CSP::
ensureLabelITS( Home home, IntVar x0, IntVar x1, IntVar x2 )
{
	if ( home.failed() )
		return;
	GECODE_ES_FAIL(ITSAtomLabel::post(home, x0, x1, x2,
		  				  educts.getGraphToSearchLabel(), products.getGraphToSearchLabel()));
}

/////////////////////////////////////////////////////////////////////////

void
ITS_CSP::
countNeighborhood( Home home
		, const IntVarArray & atomsArray
		, const std::vector < size_t > & neighID
		, const std::map< std::size_t, std::size_t > & neighCount )
{
	if ( home.failed() ) {
		return;
	}

	ViewArray< Int::IntView > atomViewArray( home, IntVarArgs(atomsArray) );
	GECODE_ES_FAIL( NeighborhoodCount::post( home, atomViewArray, neighID, neighCount ) );
}

/////////////////////////////////////////////////////////////////////////

void
ITS_CSP::
coverConnectedComps( Gecode::Home home
		, const IntVarArray & x
		, const sgm::Graph_Interface::CompLabel & compLabel
		, const size_t connNum )
{
	if ( home.failed() ) {
		return;
	}
	ViewArray< Int::IntView > atomViewArray( home, IntVarArgs( x ) );
	GECODE_ES_FAIL( ConnectedCoverage::post( home, atomViewArray, compLabel, connNum ) );
}

/////////////////////////////////////////////////////////////////////////

void
ITS_CSP::
print( void ) const
{
	std::cout << "Educts   ITS: " << eduITS << std::endl;
	std::cout << "Products ITS: " << proITS << std::endl;
	std::cout << std::endl;
}

/////////////////////////////////////////////////////////////////////////

const ReactionGraph &
ITS_CSP::
getEduGraph() const
{
	return educts;
}

/////////////////////////////////////////////////////////////////////////

const ReactionGraph &
ITS_CSP::
getProGraph() const
{
	return products;
}

/////////////////////////////////////////////////////////////////////////

const ITS&
ITS_CSP::
getITS() const
{
	return its;
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
