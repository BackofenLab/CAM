/*
 * This class represents a graph with virtual relabeled edges. It inherits
 * graph information from the graph interface and additionally controls
 * the nodes and the edges of the ITS. They form a subset of the underlying
 * full educts, products graphs. Virtual edges with the label "ITS-edge"
 * are constructed between ITS vertices-subset allowing a later matching
 * between educts and products graph.
 *
 * Created on: 2013
 *		@author Martin Mann http://www.bioinf.uni-freiburg.de/~mmann/
 *  	@author Feras Nahar <naharf@informatik.uni-freiburg.de>
 */

#ifndef ITS_GRAPH_INTERFACE_H_
#define ITS_GRAPH_INTERFACE_H_

#include "sgm/Graph_Interface.hh"
#include "sgm/Match.hh"

#include <map>
#include <set>

class ITS_Graph_Interface : public sgm::Graph_Interface {

	protected:

		  //! the label postfix added to each ITS participating atom label
		static const std::string ITSnodeLabelAdd;

		  //! The original graph of this graph object.
		const Graph_Interface & fullGraph;

		  //! all nodes participating in the ITS ring
		const std::set<int> nodesITS;

		  //! holds for each ITS nodes its neighbors in the ITS ring
		std::map<size_t, std::set<size_t> > neighInITS;

		  /*!
		   * Special edge descriptor to enable the iteration over the edges
		   * covered by this graph.
		   *
		   */
		class EdgeDescriptor : public sgm::Graph_Interface::EdgeDescriptor {

		protected:

			  //! the source of the edge described
			using Graph_Interface::EdgeDescriptor::from;

			  //! the target of the edge described
			using Graph_Interface::EdgeDescriptor::to;

			  //! the iterator to the current edge
			OutEdge_iterator curEdge;

			  //! the iterator to the edge iteration end
			OutEdge_iterator edgeEnd;

			  //! access to the parent class for full graph information
			const ITS_Graph_Interface & parent;

			  //! whether or not we are to return a virtual bond
			bool returnVirtualBondNeigh;

			  //! points to a neighbor within parent.neighITSlist
			std::set<size_t>::const_iterator curVirtualBondNeigh;

			  //! virtual label for ITS edges
			const static std::string ITSedgeLabel;

			  /*!
			   * Proceeds to next available virtual neighbor and returns whether
			   * or not this was successful.
			   * @return true = next neighbor was set and is available; false otherwise
			   */
			bool
			setNextVirtualNeighbor();

		public:
			  //! Construction
			  //! @param fromIndex the edge source index or -1 if no edge is to be described
			  //! @param curEdge the iterator to the current edge
			  //! @param edgeEnd the iterator to the iteration end
			  //! @param parent access to the parent full graph information
			EdgeDescriptor(	const int fromIndex
							, const Graph_Interface::OutEdge_iterator & curEdge
							, const Graph_Interface::OutEdge_iterator & edgeEnd
							, const ITS_Graph_Interface & parent );

			  //! Destruction
			virtual
			~EdgeDescriptor();

			  //! Access to the label of the edge.
			  //! @return the edge label
			virtual
			const std::string&
			getEdgeLabel(void) const;

			  //! Equality comparison
			  //! @param ed the edge to compare to
			  //! @return true if both descriptors describe the same edge
			virtual
			bool
			operator==(const EdgeDescriptor& ed ) const;

			  //! Inequality comparison
			  //! @param ed the edge to compare to
			  //! @return true if both descriptors describe different edges
			virtual
			bool
			operator!=(const EdgeDescriptor& ed ) const;

			  //! Inequality comparison
			  //! @param ed the edge to compare to
			  //! @return true if both descriptors describe different edges
			virtual
			bool
			operator!=(const Graph_Interface::EdgeDescriptor& ed ) const;

			  //! Iterator support
			  //! @return the next EdgeDescriptor in the adjacency
			virtual
			EdgeDescriptor&
			operator++();

			  //! Create a heap copy of this object. NOTE: this has to be
			  //! removed by the calling function.
			  //! @return a new heap copy of this object
			virtual
			EdgeDescriptor *
			clone() const;

		};

		public:

			  //! Default Constructor
			ITS_Graph_Interface();

			/*!
			 * Construction of an ITS_Graph_Interface from explicit ITS node list. .
			 *
			 * @param fullGraph the original graph which this object extends
			 *
			 * @param ITSring vector contains ITS nodes
			 */
			ITS_Graph_Interface(const Graph_Interface & fullGraph
						, const std::vector<std::size_t>& ITSring );

			  //! Destruction of the ITS graph
			virtual ~ITS_Graph_Interface();

			  //! Access to the number of nodes of the original graph
			  //! @return the overall node number
			virtual
			size_t
			getNodeNumber(void) const;

			  //! Access to iteration begin for the edge in the adjacency list of
			  //! a specified node
			  //! @param i the index of the node of interest
			  //! @return the iterator to the first edge within the adjacency of i
			virtual
			OutEdge_iterator
			getOutEdgesBegin( const IndexType & i ) const;

			  //! Access to iteration end for the edge in the adjacency list of
			  //! a specified node
			  //! @param i the index of the node of interest
			  //! @return the iterator the end of the adjacency iteration of i
			virtual
			OutEdge_iterator
			getOutEdgesEnd( const IndexType & i ) const;

			  //! Access to the label of a specified node
			  //! @param i the index of the node of interest
			  //! @return a string representation of the node label
			virtual
			std::string
			getNodeLabel(const IndexType & i) const;

			  //! Access to the wrapper molecule graph
		   	  //! @return the internally wrapped molecule graph
			const sgm::Graph_Interface&
			getMoleculeGraph() const;

};

#endif /* ITS_GRAPH_INTERFACE_H_ */
