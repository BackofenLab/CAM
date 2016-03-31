#include "ITS_Graph_Interface.h"
#include <cassert>
#include <climits>

#include <ggl/chem/MoleculeUtil.hh>

	const std::string ITS_Graph_Interface::EdgeDescriptor::ITSedgeLabel = "ITS-edge";
	const std::string ITS_Graph_Interface::ITSnodeLabelAdd = "-ITS";

	ITS_Graph_Interface::
	ITS_Graph_Interface(const Graph_Interface& fullGraph_
				, const std::vector<size_t>& ITSring)
	 :	fullGraph(fullGraph_)
		, nodesITS(ITSring.begin(), ITSring.end())
		, neighInITS()
	{
		// store all edges within ITS
		for (size_t i=0; i < ITSring.size(); ++i) {

			assert(ITSring.at(i) < (int)fullGraph.getNodeNumber());

			size_t leftNeigh = i==0 ? *ITSring.rbegin() : ITSring.at(i-1);
			size_t rightNeigh = i==(ITSring.size()-1) ? *ITSring.begin() : ITSring.at(i+1);
			neighInITS[ITSring.at(i)].insert(leftNeigh);
			neighInITS[ITSring.at(i)].insert(rightNeigh);

			assert(neighInITS[ITSring.at(i)].size() == 2);
		}

		assert(neighInITS.size() <= fullGraph.getNodeNumber());
	}

	////////////////////////////////////////////////////////////////////////////

	ITS_Graph_Interface::~ITS_Graph_Interface()
	{
	}

	////////////////////////////////////////////////////////////////////////////

	ITS_Graph_Interface::
	OutEdge_iterator
	ITS_Graph_Interface::
	getOutEdgesBegin(const IndexType & i) const
	{
		return OutEdge_iterator(
				EdgeDescriptor(
						(int)i
						, fullGraph.getOutEdgesBegin(i)
						, fullGraph.getOutEdgesEnd(i)
						, *this
				) );
	}

	////////////////////////////////////////////////////////////////////////////

	ITS_Graph_Interface::
	OutEdge_iterator
	ITS_Graph_Interface::
	getOutEdgesEnd(const IndexType & i) const
	{

		// return end of according edge iteration
		return OutEdge_iterator(
				EdgeDescriptor(
						-1
						, fullGraph.getOutEdgesEnd(i)
						, fullGraph.getOutEdgesEnd(i)
						, *this
				) );
	}

	////////////////////////////////////////////////////////////////////////////

	size_t
	ITS_Graph_Interface::
	getNodeNumber(void) const
	{
		// return the number of ITS nodes in the underlying graph
		return fullGraph.getNodeNumber();
	}

	////////////////////////////////////////////////////////////////////////////

	std::string
	ITS_Graph_Interface::
	getNodeLabel(const IndexType & i) const
	{
		  // check if ITS node
		if (nodesITS.find(i) != nodesITS.end()) {
			  // return atom label only plus ITS node label postfix
			return ggl::chem::MoleculeUtil::getAtom(fullGraph.getNodeLabel(i)) + ITSnodeLabelAdd;
		}
		  // return normal label for non-ITS nodes or non-charged ITS nodes
		return fullGraph.getNodeLabel(i);
	}

	////////////////////////////////////////////////////////////////////////////

	const sgm::Graph_Interface&
	ITS_Graph_Interface::
	getMoleculeGraph() const
	{
		return fullGraph;
	}

	////////////////////////////////////////////////////////////////////////////
	/////////////////////////////  EdgeDescriptor  /////////////////////////////
	////////////////////////////////////////////////////////////////////////////

	ITS_Graph_Interface::
	EdgeDescriptor::
	EdgeDescriptor( const int fromIndex,
					const Graph_Interface::OutEdge_iterator& cur_edge_,
				   const Graph_Interface::OutEdge_iterator& list_end_,
				   const ITS_Graph_Interface & parent_)
	 : Graph_Interface::EdgeDescriptor()
		, curEdge(cur_edge_)
		, edgeEnd(list_end_)
		, parent(parent_)
		, returnVirtualBondNeigh(false)
	{
		if (fromIndex < 0) {
			  // no more edges
			from = UINT_MAX;
			to = UINT_MAX;
		} else {
			from = (IndexType)fromIndex;
			if(curEdge != edgeEnd) {
				// check whether a node belongs to an ITS
				from = curEdge->getFromIndex();
				bool sourceIsITSnode = parent.neighInITS.find(from) != parent.neighInITS.end();
				if (sourceIsITSnode) {

					const std::set< size_t > & fromNeighbors = parent.neighInITS.find(from)->second;

					// go to first edge not part of the ITS
					bool targetIsITSmember = fromNeighbors.find(curEdge->getToIndex()) != fromNeighbors.end();
					while (curEdge != edgeEnd && targetIsITSmember) {
						++curEdge;
						if (curEdge != edgeEnd) {
							targetIsITSmember = fromNeighbors.find(curEdge->getToIndex()) != fromNeighbors.end();
						}
					}
					// check if no neighbors available
					if (curEdge == edgeEnd) {
						// handle virtual bonds part of ITS
						if (setNextVirtualNeighbor()) {
							to   = *curVirtualBondNeigh;
						} else {
							  // no virtual bond available .. ??? should not happen!
							from = UINT_MAX;
							to = UINT_MAX;
						}
					} else {
						to = curEdge->getToIndex();
					}
				} else {
					// copy edge target information
					to = curEdge->getToIndex();
				}
			} else {

				// handle virtual bonds part of ITS
				if (setNextVirtualNeighbor()) {
					to = *curVirtualBondNeigh;
				} else {
					  // no more edges
					from = UINT_MAX;
					to = UINT_MAX;
				}
			}
		}
	}

	////////////////////////////////////////////////////////////////////////////

	ITS_Graph_Interface::
	EdgeDescriptor::
	~EdgeDescriptor()
	{
	}

	////////////////////////////////////////////////////////////////////////////

	const std::string&
	ITS_Graph_Interface::
	EdgeDescriptor::getEdgeLabel(void) const {

		// return ITSedgeLabel for dummy edges
		if (from != UINT_MAX && returnVirtualBondNeigh
				&& curVirtualBondNeigh != parent.neighInITS.find(from)->second.end())
		{
			return ITSedgeLabel;
		}

		return curEdge->getEdgeLabel();
	}

	////////////////////////////////////////////////////////////////////////////

	bool
	ITS_Graph_Interface::
	EdgeDescriptor::
	operator==(const EdgeDescriptor& ed) const
	{
		return this->from == ed.from && this->to == ed.to;
	}

	////////////////////////////////////////////////////////////////////////////

	bool
	ITS_Graph_Interface::
	EdgeDescriptor::
	operator!=(const EdgeDescriptor& ed) const
	{
		return this->from != ed.from || this->to != ed.to;
	}

	////////////////////////////////////////////////////////////////////////////

	bool
	ITS_Graph_Interface::
	EdgeDescriptor::
	operator!=(const Graph_Interface::EdgeDescriptor& ed) const
	{
		// check if type casting possible
		assert(dynamic_cast<const EdgeDescriptor*>(&ed) != NULL);
		// forward call via typecast
		return this->operator !=( static_cast<const EdgeDescriptor&>(ed) );
	}

	////////////////////////////////////////////////////////////////////////////

	bool
	ITS_Graph_Interface::
	EdgeDescriptor::
	setNextVirtualNeighbor()
	{
		  // security issue
		if (from == UINT_MAX) {
			return false;
		}
		  // check if part of ITS
		if (parent.neighInITS.find(from) == parent.neighInITS.end()) {
			return false;
		}

		// no dummy set so far
		if (returnVirtualBondNeigh == false) {
			// set first neighbor
			curVirtualBondNeigh = parent.neighInITS.find(from)->second.begin();
			returnVirtualBondNeigh = true;
		} else {
			// go to next dummy if available
			if (curVirtualBondNeigh != parent.neighInITS.find(from)->second.end()) {
				++curVirtualBondNeigh;
			}
		}

		return curVirtualBondNeigh != parent.neighInITS.find(from)->second.end();
	}


	////////////////////////////////////////////////////////////////////////////

	ITS_Graph_Interface::
	EdgeDescriptor&
	ITS_Graph_Interface::
	EdgeDescriptor::
	operator++()
	{
		 // end of list reached -> nothing left to be done
		if ( from == UINT_MAX) {
			return *this;
		}

		 // handle non-ITS source nodes
		if ( parent.neighInITS.find(from) == parent.neighInITS.end()) {
			if ( curEdge != edgeEnd ) {
				++curEdge;
			}
			if (curEdge == edgeEnd) {
				from = UINT_MAX;
				to = UINT_MAX;
			} else {
				to = curEdge->getToIndex();
			}
			return *this;
		}

		//////////  ITS SOURCE NODE !!!!!!!!!!!!!!

		// check if end of real list already reached
		if (curEdge == edgeEnd) {

			// handle virtual bonds part of ITS
			if (setNextVirtualNeighbor()) {
				to = *curVirtualBondNeigh;
			} else {
				from = UINT_MAX;
				to = UINT_MAX;
			}

			return *this;
		}

		 // short access to neighbors of current source node
		const std::set< size_t > & fromNeighbors = parent.neighInITS.find(from)->second;

		// proceed to next edge part of the graph
		do {
			++curEdge;
		} while (curEdge != edgeEnd
			&& fromNeighbors.find(curEdge->getToIndex()) != fromNeighbors.end());


		// check if end of list already reached
		if (curEdge == edgeEnd) {
			// handle virtual bonds part of ITS
			if (setNextVirtualNeighbor()) {
				to = *curVirtualBondNeigh;
			} else {
				from = UINT_MAX;
				to = UINT_MAX;
			}
		} else {
			// handle normal non-ITS target nodes
			to = curEdge->getToIndex();
		}

		return *this;
	}

	////////////////////////////////////////////////////////////////////////////

	ITS_Graph_Interface::EdgeDescriptor*
	ITS_Graph_Interface::
	EdgeDescriptor::
	clone() const
	{
		return new EdgeDescriptor(*this);
	}

	////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////// END //////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////
