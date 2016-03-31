
#include "ReactionGraph.h"

#include <climits>

#include <boost/lexical_cast.hpp>

#include <sgm/MR_Storing.hh>
#include <sgm/GM_vf2.hh>

#include <ggl/chem/MoleculeUtil.hh>

using namespace ggl::chem;

///////////////////////////////////////////////////////////////////////////

std::map< std::string, size_t > ReactionGraph::atomLabel2ID;
std::map< size_t, std::string > ReactionGraph::id2AtomLabel;
std::map< ReactionGraph::ValenceNeighPair, std::size_t > ReactionGraph::mapNeighborhood2ID;

///////////////////////////////////////////////////////////////////////////

ReactionGraph::ReactionGraph( const sgm::Graph_Interface & molGraphs
			, const bool createComponenLabeling
			, const bool _verbose)
	: molecules()
	, molGraph(NULL)
	, molIndexOrder( molGraphs.getNodeNumber() )
	, molIndexProtons( molIndexOrder.size() )
	, molIndexObsoleteProtons( molIndexOrder.size() )
    , searchGraph(NULL)
	, graphToSearchSymmetries()
	, graphToSearchValences()
	, graphToSearchLabel()
	, graphToSearchCompLabel()
	, graphToSearchCompNum(0)
	, verbose(_verbose)
{
	  // make a copy of the molecule graphs
	MoleculeUtil::copy( molGraphs, molecules );
	molGraph = new Molecule_Graph( molecules );

	  // list of obsolete protons indices with an according mapping to the
	  // "master proton" to be used instead
	NodeMap obsoleteProtons = getObsoleteProtons( *molGraph );

	  // generate node-ordered graph where all H-atoms are with highest index
	  // and all obsolete protons are at the end
	molIndexObsoleteProtons = molIndexOrder.size() - obsoleteProtons.size();
	size_t posNonH = 0;
	size_t posHobs = molIndexObsoleteProtons;
	size_t posH = posHobs - 1;

	for (size_t i = 0; i < molIndexOrder.size(); ++i) {
		  // check if proton
		std::string atomLabel = ggl::chem::MoleculeUtil::getAtom(molGraph->getNodeLabel(i));
		  // store mapping
		if (atomLabel2ID.find(atomLabel) == atomLabel2ID.end()) {
			size_t graphToSearchLabelNumber = atomLabel2ID.size();
			atomLabel2ID[atomLabel] = graphToSearchLabelNumber;
			id2AtomLabel[graphToSearchLabelNumber] = atomLabel;

		}

		if (atomLabel.compare("H") == 0) {
			  // check if obsolete proton
			if (obsoleteProtons.find(i) != obsoleteProtons.end()) {
				molIndexOrder[posHobs] = i;
				++posHobs;
			} else {
				  // master proton
				molIndexOrder[posH] = i;
				--posH;
			}
		} else {
			  // no proton
			molIndexOrder[posNonH] = i;
			++posNonH;
		}
	}
	// store first index of a proton index
	molIndexProtons = posNonH;

	  // store label encoding
	graphToSearchLabel.resize( molIndexObsoleteProtons, 0 );
	for ( size_t i = 0; i < molIndexObsoleteProtons; ++i )
	{
		graphToSearchLabel[i] = (int) atomLabel2ID.find(ggl::chem::MoleculeUtil::
				getAtom(molGraph->getNodeLabel(molIndexOrder.at(i))))->second;
	}

	  // get symmetries of graph to match and store
	sgm::MR_StoringInsertT< SymmetricMatchSet > matchStoring(graphToSearchSymmetries);
	sgm::GM_vf2 matcher;
	sgm::SubGraph::NodeList matchNodes(molIndexOrder.begin(), molIndexOrder.begin() + molIndexObsoleteProtons);
	// sgm::SubGraph searchGraph( getGraphToSearch() );
	searchGraph = new sgm::SubGraph( *molGraph, matchNodes );
	matcher.findMatches( sgm::Pattern( *searchGraph ), *searchGraph, matchStoring, UINT_MAX );

	  // get charge and bond valence information encoding
    graphToSearchValences = ValenceMatrix( *searchGraph );

	if ( verbose ) {
		std::cout << *searchGraph << std::endl;
		graphToSearchValences.printMatrix();
	}

	  // check if connected component information is needed
	if (createComponenLabeling) {
		  // get connected component labeling and number
		graphToSearchCompNum = sgm::Graph_Interface::connectedComponents( *searchGraph, graphToSearchCompLabel );
	} else {
		  // set dummy connected component labeling and number
		graphToSearchCompNum = 1;
		graphToSearchCompLabel.resize(searchGraph->getNodeNumber(), 0);
	}

	  // construct list of out-edges for each atom in the graph
	edgeList.resize( graphToSearchValences.numRows() );
	  // store the numbers of out-edges in a vector
	for (size_t r = 0; r < graphToSearchValences.numRows(); r++)
	{
		for (size_t c = 0; c < graphToSearchValences.numColumns(); c++)
		{	 /*!
			  * if a cell is not null that means, there is a valence
			  * and accordingly an edge exists. exclude charge at
			  * diagonal entries (r = c) if exists.
			  */
			if ( graphToSearchValences.at(r, c) != 0 && (r != c) )
				edgeList[r].insert( graphToSearchValences.at( r, c ) );

		}// inner for
	}// outer for

	fillNeighborhoodSet();
}

///////////////////////////////////////////////////////////////////////////

ReactionGraph::~ReactionGraph() {
	if (molGraph != NULL) {
		delete molGraph;
		molGraph = NULL;
	}

	if ( searchGraph != NULL ) {
		delete searchGraph;
		searchGraph = NULL;
	}
}


///////////////////////////////////////////////////////////////////////////

const ggl::chem::Molecule_Graph &
ReactionGraph::
getMoleculeGraph() const
{
	return *molGraph;
}

///////////////////////////////////////////////////////////////////////////

const
sgm::SubGraph &
ReactionGraph::
getGraphToSearch() const
{
	  // get list of nodes to match
//	sgm::SubGraph::NodeList matchNodes(molIndexOrder.begin(), molIndexOrder.begin() + molIndexObsoleteProtons);
	  // return according subgraph
//	return sgm::SubGraph( *molGraph, matchNodes );
	return *searchGraph;
}

///////////////////////////////////////////////////////////////////////////

const size_t
ReactionGraph::
getGraphToSearchSize() const {
	return molIndexObsoleteProtons;
}

///////////////////////////////////////////////////////////////////////////

const size_t
ReactionGraph::
getGraphToSearchProtonNumber() const {
	return molIndexObsoleteProtons-molIndexProtons;
}

///////////////////////////////////////////////////////////////////////////

const ReactionGraph::SymmetricMatchSet&
ReactionGraph::
getGraphToSearchSymmetries() const {
	return graphToSearchSymmetries;
}

///////////////////////////////////////////////////////////////////////////

const ValenceMatrix&
ReactionGraph::
getGraphToSearchValences() const {
	return graphToSearchValences;
}

///////////////////////////////////////////////////////////////////////////

const std::vector<int> &
ReactionGraph::
getGraphToSearchLabel() const {

	return graphToSearchLabel;
}

///////////////////////////////////////////////////////////////////////////

const size_t
ReactionGraph::
getAtomLabelNumber() {
	return atomLabel2ID.size();
}

///////////////////////////////////////////////////////////////////////////

ggl::chem::Molecule
ReactionGraph::
getGraphToMatch( const sgm::Match & itsMatch, sgm::Match & newITSindices ) const
{
	   // get molecule with specially labeled ITS-atoms
	 Molecule mol;
	 // MoleculeUtil::copy( getGraphToSearch(), mol);
	 MoleculeUtil::copy( *searchGraph, mol);

	 // TODO : add class ID with fixed length and prune again to support given class information

	 boost::property_map< Molecule, PropNodeLabel >::type
		 label = boost::get( PropNodeLabel(), mol );
	   // relabel ITS participating atoms
	 for ( size_t i = 0; i < itsMatch.size(); ++i ) {
		 std::string nodeLabel = label[boost::vertex( itsMatch.at(i), mol )];
		 assert(nodeLabel.find(':')==std::string::npos);
		 label[boost::vertex( itsMatch.at(i), mol )] = nodeLabel + ":" + boost::lexical_cast< std::string >( i + 1 );
	 }

	   // compress non-ITS protons into adjacent atom labels
	 MoleculeUtil::compressHnodes( mol );

	   // get new mapping and undo special marking
	 newITSindices.resize(itsMatch.size());
	 for (size_t i=0; i < boost::num_vertices(mol); ++i) {
		 std::string nodeLabel = label[boost::vertex( i, mol )];
		 int curClass = MoleculeUtil::getClass( nodeLabel );
		 if (curClass > 0) {
			   // cut off class information
			 label[boost::vertex( i, mol )] = nodeLabel.substr(0,nodeLabel.find(':'));
			   // store ITS index information
			 newITSindices[curClass-1] = i;
		 }
	 }

	 return mol;
}

///////////////////////////////////////////////////////////////////////////

ReactionGraph::NodeMap
ReactionGraph::
getObsoleteProtons( const sgm::Graph_Interface& mol )
{
	NodeMap obsolete;
	 // stores the information for non-protons if a proton was already added or not
	std::vector< int > alreadyAdded( mol.getNodeNumber(), -1);

	for (size_t i=0; i < mol.getNodeNumber(); ++i) {
		if (mol.getNodeLabel(i).compare("H") == 0) {
			sgm::Graph_Interface::OutEdge_iterator bond = mol.getOutEdgesBegin(i);
			assert( bond != mol.getOutEdgesEnd(i) );
			  // check if for the adjacent node a proton was already touched
			if ( alreadyAdded.at(bond->getToIndex()) >= 0 ) {
				  // note that this proton can be ignored and store replacement
				obsolete[i] = (size_t)alreadyAdded.at(bond->getToIndex());
			} else {
				  // mark that one proton was already ignored for this node
				alreadyAdded[bond->getToIndex()] = (int)i;
			}
		}
	}
	return obsolete;
}

///////////////////////////////////////////////////////////////////////////

void
ReactionGraph::
fillNeighborhoodSet()
{
	neighboursID.resize( graphToSearchValences.numRows(), 0 );
	 /*!
	  * locate neighbor atoms and their bond valences from each atom in
	  * the educts and products:
	  * (Valence, Neighbor) pair for each educts/products atom
	  */
	ValenceNeighPair valenceNeigh;
	size_t atomLabID, neigh ;
	std::vector< ValenceNeighPair > neighVect;
	for ( size_t r = 0; r < graphToSearchValences.numRows(); r++ )
	{
		atomLabID = atomLabel2ID.find( ggl::chem::MoleculeUtil::getAtom( searchGraph->getNodeLabel( r )) )->second;
		for ( size_t c = 0; c < graphToSearchValences.numColumns(); c++ )
		{
			 if ( r != c && graphToSearchValences.at(r, c) != 0 ) { // if holds, there is an edge-valence
			 	  // store neighbor's atom label ID
				neigh = atomLabel2ID.find( ggl::chem::MoleculeUtil::getAtom( searchGraph->getNodeLabel( c )) )->second;
			 	  // saving the neighbor atom label together with the valence
			 	valenceNeigh.insert( std::pair< int, size_t > ( graphToSearchValences.at(r, c), neigh ) );
			 }
		}

		/*!
		 * constructing the local neighborhood set of the currently explored atom
		 * i.e. inserting its neighbors and edge valences into it
		 */
		neighborhood.insert( EdgePattern( atomLabID, valenceNeigh ) );
		  // assigning ID to each edge pattern
		if (mapNeighborhood2ID.find(valenceNeigh) == mapNeighborhood2ID.end()) {
			size_t valenceNeighNumber = mapNeighborhood2ID.size();
			mapNeighborhood2ID[valenceNeigh] = valenceNeighNumber;
		}
		  // storing the edge pattern for each atom to replace them later by their IDs
		neighVect.push_back( valenceNeigh );
		valenceNeigh.clear();
	}

	  // designate each atom now with ID of its edge pattern
	for ( size_t i = 0; i < graphToSearchValences.numRows(); i++ )
		neighboursID[i] = mapNeighborhood2ID[ neighVect[i] ];

	// printNeighborhood( neighborhood );
}

/////////////////////////////////////////////////////////////////////////

const
ReactionGraph::LocalNeighborhood &
ReactionGraph::
getNeighborhood() const
{
	return neighborhood;
}

/////////////////////////////////////////////////////////////////////////

ReactionGraph::LocalNeighborhood
ReactionGraph::
substract( const ReactionGraph::LocalNeighborhood & neighborhood1
		   , const ReactionGraph::LocalNeighborhood & neighborhood2 )
{
	ReactionGraph::LocalNeighborhood neighDiff;
	std::set_difference( neighborhood1.begin(), neighborhood1.end()
						, neighborhood2.begin(), neighborhood2.end()
						, std::inserter(neighDiff, neighDiff.end()) );
	return neighDiff;
}

///////////////////////////////////////////////////////////////////////////

std::map< std::size_t, std::size_t >
ReactionGraph::
getITS_Members( const ReactionGraph::LocalNeighborhood & neighborhood1
	    		 , const ReactionGraph::LocalNeighborhood & neighborhood2 )
{
	  // store label IDs of ITS-participated atoms with corresponding number of appearances
	std::map< std::size_t, std::size_t > itsLabelCount ;
	ReactionGraph::LocalNeighborhood neighDiff = ReactionGraph::substract( neighborhood1, neighborhood2 );
	ReactionGraph::LocalNeighborhood::const_iterator it;
	for ( it = neighDiff.begin(); it != neighDiff.end(); ++it ) {
		itsLabelCount[(*it).first]++;
	}

//std::map< std::size_t, std::size_t >::iterator i;
//for ( i = itsLabelCount.begin(); i != itsLabelCount.end(); ++i ) {
//	std::cout << (*i).first << "-" << (*i).second << std::endl;
//}

	return itsLabelCount;
}

///////////////////////////////////////////////////////////////////////////

std::map< std::size_t, std::size_t >
ReactionGraph::
getNeighCount( const ReactionGraph::LocalNeighborhood & neighborhood1
	    	    , const ReactionGraph::LocalNeighborhood & neighborhood2 )
{
	  // Maps each neighborhood ID to its number of appearances
    std::map< std::size_t, std::size_t > neighCount ;
	ReactionGraph::LocalNeighborhood neighDiff = ReactionGraph::substract( neighborhood1, neighborhood2 );

	LocalNeighborhood::const_iterator it1;
	for ( it1 = neighDiff.begin(); it1 != neighDiff.end(); ++it1 ) {
	      // count up the number of resulted edge-patterns based on their IDs
		neighCount[ mapNeighborhood2ID[(*it1).second] ]++;
	}

//std::map< std::size_t, std::size_t >::iterator i;
//for ( i = neighCount.begin(); i != neighCount.end(); ++i ) {
//	std::cout << (*i).first << "-" << (*i).second << std::endl;
//}

	return neighCount;
}

///////////////////////////////////////////////////////////////////////////
void
ReactionGraph::
printNeighborhood ( const LocalNeighborhood & neighborhood )
{
	LocalNeighborhood::const_iterator iter;
	ValenceNeighPair::const_iterator vnIter;
	for ( iter = neighborhood.begin(); iter != neighborhood.end(); ++iter ) {

		std::cout << id2AtomLabel.find((*iter).first)->second << ", { ";
		ValenceNeighPair vnSet = (*iter).second;
		for (vnIter = vnSet.begin(); vnIter != vnSet.end(); ++vnIter) {
			std::cout << "(" << (*vnIter).first << ", " << id2AtomLabel.find((*vnIter).second)->second << ") ";
		}
		std::cout << "}" << std::endl;
	}
	std::cout << std::endl;
}

///////////////////////////////////////////////////////////////////////////

const ReactionGraph::EdgeList &
ReactionGraph::
getEdgeList() const
{
   	return edgeList;
}

///////////////////////////////////////////////////////////////////////////

const ValenceMatrix &
ReactionGraph::
getValenceMatrix() const
{
	return graphToSearchValences;
}

///////////////////////////////////////////////////////////////////////////

bool
ReactionGraph::
isHomovalent() const
{
	for ( size_t r = 0; r < graphToSearchValences.numRows(); r++ )
		if ( graphToSearchValences.at(r, r) != 0 )
			return false;
	return true;
}

///////////////////////////////////////////////////////////////////////////

const std::vector< std::size_t > &
ReactionGraph::
getNeighIDs() const
{
	return neighboursID;
}

///////////////////////////////////////////////////////////////////////////

const ReactionGraph::CompLabel &
ReactionGraph::
getGraphToSearchCompLabel() const {
	return graphToSearchCompLabel;
}

///////////////////////////////////////////////////////////////////////////

const size_t
ReactionGraph::
getGraphToSearchCompNum() const {
	return graphToSearchCompNum;
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////



