#ifndef REACTIONGRAPH_H_
#define REACTIONGRAPH_H_

#include <map>
#include <set>

#include <sgm/SubGraph.hh>

#include <ggl/chem/Molecule.hh>

#include "ValenceMatrix.h"

/*!
 * This class represents either an educt or a product graph participating in
 * a reaction.
 *
 * Furthermore, it provides additional information.
 *
 */
class ReactionGraph {

public:

	typedef std::map<size_t, size_t> NodeMap;

	  //! vector of symmetries
	typedef std::set< sgm::Match > SymmetricMatchSet;

	  //! Vector containing out-edges of a node in the underlying graph
	typedef std::vector< std::multiset< int > > EdgeList;

	  //! Pair defines valence and atom label ID representing neighbor atom
	typedef std::multiset< std::pair < int, std::size_t > > ValenceNeighPair;

	  //! Pair defines an atom with its associated neighbor atom and corresponding valence
	typedef std::pair < std::size_t, ValenceNeighPair > EdgePattern;

	  //! Set which stores the neighbor atoms together with their linkages
	typedef std::multiset< EdgePattern > LocalNeighborhood;

	typedef sgm::Graph_Interface::CompLabel CompLabel;


protected:

	  //! global mapping of atom strings to IDs --> valid for all instances
	static std::map< std::string, size_t > atomLabel2ID;

	 //! reverse mapping of reaction graph AtomLabels --> IDs; used for print purposes
	static std::map < std::size_t, std::string> id2AtomLabel;

	  //! a single graph representation of all molecules
	ggl::chem::Molecule molecules;

	  //! wrapped representation of molecules
	ggl::chem::Molecule_Graph * molGraph;

	  /*! index shuffling of the molecules' nodes where protons are
	   * shifted to the end.
	   */
	sgm::SubGraph::NodeList molIndexOrder;

	  /*! first index of a proton within molIndexOrder. All following
	   * indices are protons too
	   */
	size_t molIndexProtons;

	  /*! first index of an obsolete proton within molIndexOrder. All following
	   * indices are obsolete protons too
	   */
	size_t molIndexObsoleteProtons;

	  /*!
	   * Graph representation excluding obsolete protons
	   */
	sgm::SubGraph * searchGraph;

	  //! Symmetries of the graph accessible via getGraphToSearch()
	SymmetricMatchSet graphToSearchSymmetries;

	  //! charge and bond valence information encoded as adjacency information
	ValenceMatrix graphToSearchValences;

	  //! label encoding for each atom
	std::vector< int > graphToSearchLabel;

	  //! connected component labeling
	CompLabel graphToSearchCompLabel;

	  //! number of connected components
	size_t graphToSearchCompNum;

	  //! Sorted lists of edge valences for each node
	EdgeList edgeList;

	  //! Neighborhood set of the underlying reaction graph
	LocalNeighborhood neighborhood;

	  //! Associate all neighborhood sets with IDs
	static std::map< ValenceNeighPair, std::size_t > mapNeighborhood2ID;

	  //! Associate each atom with its neighborhood set ID from the map above
	std::vector< std::size_t > neighboursID;

	  //! whether or not to do verbose output
	const bool verbose;

public:

	 /*!
	  * Construction
	  *
	  * @param moleculeGraph the molecule graph covering all molecules
	  * @param createComponenLabeling whether or not a connected component
	  *        labeling is to be constructed.
	  * @param verbose whether or not to do verbose reporting
	  */
	ReactionGraph( const sgm::Graph_Interface & moleculeGraphs
					, const bool createComponenLabeling
					, const bool verbose );


	 //! destruction
	virtual ~ReactionGraph();


	  /*!
	   * Access to the original graph containing all atoms etc.
	   * @return the molecule graph representing one side of a reaction
	   */
	const ggl::chem::Molecule_Graph &
	getMoleculeGraph() const;


	  /*!
	   * The nodes and edges of the molecules to be considered for an ITS
	   * matching. Obsolete protons are not part of this graph.
	   * @return the graph to be matched
	   */
	const
	sgm::SubGraph &
	getGraphToSearch() const;


	  /*!
	   * Number of nodes of the graph accessible via getGraphToSearch().
	   * @return the number of nodes of the graph
	   */
	const size_t
	getGraphToSearchSize() const;

	  /*!
	   * Number of protons within the graph accessible via geGraphToSearch().
	   * @return the number of protons within the graph
	   */
	const size_t
	getGraphToSearchProtonNumber() const;

	  /*!
	   * Symmetries of the graph accessible via getGraphToSearch() excluding
	   * identity.
	   * @return the set of symmetries
	   */
	const SymmetricMatchSet&
	getGraphToSearchSymmetries() const;


	  /*!
	   * Charge and bond valence information of the graph accessible via
	   * getGraphToSearch().
	   * @return the charge and bond valence information
	   */
	const ValenceMatrix&
	getGraphToSearchValences() const;


	  /*!
	   * Atom label ID encoding in range [0,getGraphToSearchLabelNumber())
	   * @return vector of label IDs for each node
	   */
	const std::vector<int> &
	getGraphToSearchLabel() const;


	  /*!
	   * Connected component ID for each atom node.
	   * @return Connected component assignment for each node index.
	   */
	const CompLabel &
	getGraphToSearchCompLabel() const;

	  /*!
	   * Number of connected components.
	   * @return Number of connected components
	   */
	const size_t
	getGraphToSearchCompNum() const;


	  /*!
	   * Number of different atom labels within atomLabel2ID
	   * @return the number of known atom labels
	   */
	static
	const size_t
	getAtomLabelNumber();


	  /*!
	   * Returns reaction graph encoding for a given ITS match where all non-ITS
	   * protons are compressed into atom labels.
	   * @param itsMatch the nodes participating in the ITS
	   * @param newITSindices the indices of the ITS nodes in the returned graph
	   * @return the compressed graph
	   */
	ggl::chem::Molecule
	getGraphToMatch( const sgm::Match & itsMatch, sgm::Match & newITSindices ) const;


	   /*!
	    * Get the degree/edge information of each node in the reaction graph
	    * @return vector of ordered edge valence information for each atom
	    */
	const EdgeList &
	getEdgeList() const;


	  /*!
	   * Access to the sparse matrix representation of the graph
	   * @return graphToSearchValences sparse matrix
	   */
	const ValenceMatrix &
	getValenceMatrix() const;


	   /*!
	    * Check the valence matrix for homovalence in the underlying reaction graph
	    * @return true if homovalent, false if not
	    */
	bool
	isHomovalent() const;


       /*!
        * Access to the local neighborhood set of the underlying reaction graph
        * @return local neighborhood set
        */
	const
	LocalNeighborhood &
	getNeighborhood() const;


	  /*!
	   * Print out the overall neighborhood-set (all edge patterns)
	   * @param LocalNeighborhood to print (educts/products neighborhoods)
	   */
	static
	void
	printNeighborhood ( const LocalNeighborhood & neighborhood ) ;


	 /*!
	  * Compute ITS-participated atoms by subtracting neighborhood sets: educts\products or products\educts
	  * @param neighborhood1 neighborhood set of educts/products
	  * @param neighborhood2 neighborhood set of educts/products
	  * @return map each neighborhood ID to its number of occurrences
	  */
	static
	LocalNeighborhood
	substract( const ReactionGraph::LocalNeighborhood & neighborhood1
			  , const ReactionGraph::LocalNeighborhood & neighborhood2 );


	  /*!
	   * Access to the vector holding the identifiers of educts neighborhood sets
	   * @return vector of neighborhood IDs of each atom in the reaction graph
	   */
	const
	std::vector< std::size_t > &
	getNeighIDs() const;

	  /*!
	   * Access to the map that associates each label of ITS-membered atoms
	   * with its number of appearances in the ITS ring.
	   * @param neighborhood1 neighborhood set of educts/products
	   * @param neighborhood2 neighborhood set of educts/products
	   * @return map (label id, its count) in the ITS ring.
	   */
	static
	std::map< std::size_t, std::size_t >
	getITS_Members( const ReactionGraph::LocalNeighborhood & neighborhood1
		    		 , const ReactionGraph::LocalNeighborhood & neighborhood2 );


	 /*!
	  * Access to the map that associates each neighborhood identifier
	  * with its number of occurrences.
	  * @param neighborhood1 neighborhood set of educts/products
	  * @param neighborhood2 neighborhood set of educts/products
	  * @return map (neighborhood id, its count) for the underlying graph.
	  */
	static
	std::map< std::size_t, std::size_t >
	getNeighCount( const ReactionGraph::LocalNeighborhood & neighborhood1
			       , const ReactionGraph::LocalNeighborhood & neighborhood2 );

protected:

	  /*!
	   * computes the set of proton node indices that are obsolete for the
	   * ITS identification since there is at least one other "master" proton
	   * (not in the list) that is adjacent to the same non-proton.
	   * @param molecule the molecule graph to analyze
	   * @return the list of obsolete proton indices and their according master
	   * proton
	   */
	static
	NodeMap
	getObsoleteProtons( const sgm::Graph_Interface& mol );

	 /*!
	  * For each atom in the educts and the products, determines its neighbor
	  * atoms (label ID of them) in addition to the edge valences of the explored
	  * atom to them. It stores the resulted local neighborhood in form of sets.
	  */
	void
	fillNeighborhoodSet( );

};

#endif /* REACTIONGRAPH_H_ */
