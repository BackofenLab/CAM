
#include "ITS_CSP.h"

using namespace Gecode;


//////////////////////////////////////////////////////////////////////////////


ITS_CSP::
ITS_CSP( const ITS & its
		, const ReactionGraph & educts
		, const ReactionGraph & products
		, ITS_CSP::SymmetryHandler & symHandler
) :	its(its)
	, educts(educts)
	, products(products)
	, symHandler(symHandler)
	, eduITS( *this, its.getSize(), 0, educts.getGraphToSearchSize()-1 )
	, eduITSproton( *this, its.getSize(), 0, 1 )
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

		  // store whether or not a mapping includes protons
		rel( *this, eduITS[i], IRT_GQ, eductsFirstProtonIndex, eduITSproton[i] );
	}

	// restrict number of protons within ITS to be less than half
	const size_t maxProtonNumber = (its.getSize() -1) / 2;
	linear( *this, eduITSproton, IRT_LQ, maxProtonNumber );


	for ( size_t i = 0; i < (its.getSize()-1); i++ )
	{
		  // ensure alternating cycle structure of the ITS in the mapping
		alternateCycle( *this, eduITS[i], eduITS[i+1], proITS[i], proITS[i+1], its.getBondChange().at(i) );

		  // ensure minimal edge valence to speedup alternating cycle
		if ( its.getBondChange().at(i) > 0 )
			minEdgeValence( *this, proITS[i], proITS[i+1], abs(its.getBondChange().at(i)), products.getGraphToSearchValences() );
		else if ( its.getBondChange().at(i) < 0 )
			minEdgeValence( *this, eduITS[i], eduITS[i+1], abs(its.getBondChange().at(i)), educts.getGraphToSearchValences() );
	}

	  // ensure alternating cycle for ring closure
	alternateCycle( *this, eduITS[its.getSize()-1], eduITS[0], proITS[its.getSize()-1], proITS[0],
			its.getBondChange().at(its.getSize()-1) );

	  // ensure minimal edge valence to speedup alternating cycle for the last bond pair
	if ( its.getBondChange().at(its.getSize()-1) > 0 )
		minEdgeValence( *this, proITS[its.getSize()-1], proITS[0], abs(its.getBondChange().at(its.getSize()-1)), products.getGraphToSearchValences() );
	else if ( its.getBondChange().at(its.getSize()-1) < 0 )
		minEdgeValence( *this, eduITS[its.getSize()-1], eduITS[0], abs(its.getBondChange().at(its.getSize()-1)), educts.getGraphToSearchValences() );

	  // ensure that all connected components are included in the mapping
	if ( educts.getGraphToSearchCompNum() > 1 )
	{
		coverConnectedComps( *this, eduITS, educts.getGraphToSearchCompLabel(), educts.getGraphToSearchCompNum() );
	}
	if ( products.getGraphToSearchCompNum() > 1 )
	{
		coverConnectedComps( *this, proITS, products.getGraphToSearchCompLabel(), products.getGraphToSearchCompNum() );
	}

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

	  // order constraints to break ITS symmetries
	sgm::PA_OrderCheck::CheckList::const_iterator iter;
	for ( iter = its.getNeededOrderChecks().begin(); iter != its.getNeededOrderChecks().end(); iter++ ) {
		rel(*this, eduITS[(*iter).first], IRT_LE, eduITS[(*iter).second]);
	}
	  // branching on educts
	branch(*this, eduITS, INT_VAR_SIZE_MIN(), INT_VAL_MIN());
	  // post symmetry checker for educt assignments
	branch( *this, &checkEduSymmetries );
	  // branching on products
	branch(*this, proITS, INT_VAR_SIZE_MIN(), INT_VAL_MIN());
	  // post symmetry checker for product assignments (educt assignment assumed)
	branch( *this, &checkProSymmetries );

	/*
	//////////////////////////////////////////////////////////
	// POSSIBLE HEURISTICS
	 * # if multi molecule : (i,j) of bond break/formation to be tried first to
	 *   span different molecules
	 * # identify Exact Subgraphs Matches (ESM) that cover >50% of the molecules
	 *   for large reactions
	 *   - shrink ESM graphs by one bond and avoid resulting subgraph in CAM
	 * # ensure that "reactive" atoms as "O,N,S" are involved (not only C,H) : ! NOT ALWAYS TRUE : R00930, R01344, R00353, R00743
	 */



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
		, symHandler( itsCSP.symHandler )

{
	  // modifying data members during the copying
	eduITS.update( *this, share, itsCSP.eduITS );
	eduITSproton.update( *this, share, itsCSP.eduITSproton );
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

sgm::Match
ITS_CSP::
getEduITSsolution() const {
	  // get assignment
	sgm::Match curSol = getEductITS();
	SymmetricMatchSet allSym;
	  // generate all symmetries including identity
	storeSymmetries( curSol, educts, allSym, false );
	assert(!allSym.empty());
	  // return smallest symmetry of all
	return *(allSym.begin());
}

/////////////////////////////////////////////////////////////////////////

sgm::Match
ITS_CSP::
getProITSsolution() const {
	  // get assignment
	sgm::Match curSol = getProductITS();
	SymmetricMatchSet allSym;
	  // generate all symmetries including identity
	storeSymmetries( curSol, products, allSym, false );
	assert(!allSym.empty());
	  // return smallest symmetry of all
	return *(allSym.begin());
}

/////////////////////////////////////////////////////////////////////////

void
ITS_CSP::
storeSymmetries( const sgm::Match& curSol
			, const ReactionGraph & graph
			, ITS_CSP::SymmetricMatchSet & symSolStorage
			, const bool removeIdentity ) const
{
	  // set of all unique symmetric solutions
	SymmetricMatchSet symSol;

	  // add current solution if not to be ignored
	if (!removeIdentity) {
		symSol.insert(curSol);
	}

	  // get symmetries according to ITS variable shuffling
	getITSSymmetries( curSol, its, symSol );

	  // generate all value symmetric solutions and store
	sgm::Match nextSym(curSol.size());
	for (ReactionGraph::SymmetricMatchSet::const_iterator
			curSym = graph.getGraphToSearchSymmetries().begin();
			curSym != graph.getGraphToSearchSymmetries().end(); ++curSym)
	{
		  // get symmetric "sub-match", ie. symmetric matching positions of current assignment
		for (size_t i=0; i<curSol.size(); ++i) {
			assert( curSol.at(i) < curSym->size() );
			nextSym[i] = curSym->at( curSol.at(i) );
		}
		  // insert to set (ensures uniqueness of symmetric solutions)
		symSol.insert(nextSym);
		  // get symmetries according to ITS variable shuffling
		getITSSymmetries( nextSym, its, symSol );
	}
	  // store symmetric solutions
	for (SymmetricMatchSet::const_iterator curSym = symSol.begin(); curSym != symSol.end(); ++curSym) {
		  // check if this is the identity symmetry
		bool identic = removeIdentity;
		for (size_t i=0; identic && i<curSym->size(); ++i) {
			identic = curSym->at(i) == curSol.at(i);
		}
		  // store only non-identity symmetric solutions if requested
		if (!identic) {
			symSolStorage.insert( *curSym );
		}
	}
}


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

void
ITS_CSP::
checkEduSymmetries( Gecode::Space & home )
{

	ITS_CSP & csp = dynamic_cast< ITS_CSP &>(home);
	  // ensure propagation and check for failure
	if ( csp.status() == SS_FAILED ) {
		return;
	}

	  // if this a symmetric sub solution -> fail the search here
	if (csp.symHandler.eduSymSolutions.find( csp.getEductITS() ) != csp.symHandler.eduSymSolutions.end()) {
		csp.fail();
	}
	  // clear previous product symmetries, since we swapped to another educt assignment
	if (csp.symHandler.eduLastSolution != csp.getEductITS()) {
		csp.symHandler.proSymSolutions.clear();
	}
}

/////////////////////////////////////////////////////////////////////////

void
ITS_CSP::
checkProSymmetries( Gecode::Space & home )
{
	ITS_CSP & csp = dynamic_cast< ITS_CSP &>(home);
	  // ensure propagation and check for failure
	if ( csp.status() == SS_FAILED ) {
		return;
	}

	  // if this a symmetric sub solution -> fail the search here
	if (csp.symHandler.proSymSolutions.find( csp.getProductITS() ) != csp.symHandler.proSymSolutions.end()) {
		csp.fail();
	} else {
		  // try to generate symmetric educt solution only for first product assignment
		if (csp.symHandler.eduLastSolution != csp.getEductITS()) {
			  // store current educt assignment
			csp.symHandler.eduLastSolution = csp.getEductITS();
			  // store all symmetric solutions of educts
			csp.storeSymmetries( csp.symHandler.eduLastSolution, csp.educts, csp.symHandler.eduSymSolutions, true );
		}
		  // store all symmetric solutions of products
		csp.storeSymmetries( csp.getProductITS(), csp.products, csp.symHandler.proSymSolutions, true );
	}
}


/////////////////////////////////////////////////////////////////////////

void
ITS_CSP::
getITSSymmetries( const sgm::Match & curSol
					, const ITS & its
					, SymmetricMatchSet & symSol )
{
	  // generate ITS-based symmetries of this match symmetry
	const ITS::MatchSet & itsSym = its.getSymmetricMatches();
	sgm::Match nextSym(curSol.size());
	for(ITS::MatchSet::const_iterator sym = itsSym.begin(); sym != itsSym.end(); ++sym ) {
		assert(sym->size() == curSol.size());
		  // create symmetric solution via variable swapping
		for (size_t i=0; i<sym->size(); ++i) {
			nextSym[i] = curSol.at( sym->at(i) );
		}
		  // store symmetry of the symmetric solution
		symSol.insert(nextSym);
	}

}

/////////////////////////////////////////////////////////////////////////

const
ITS_CSP::SymmetricMatchSet &
ITS_CSP::getProSymmetries() const
{
	return symHandler.proSymSolutions;
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
