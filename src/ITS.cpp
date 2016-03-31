
#include "ITS.h"

#include <cstdlib>
#include <stdexcept>
#include <boost/lexical_cast.hpp>
#include <ggl/chem/MoleculeUtil.hh>

using namespace ggl::chem;


//////////////////////////////////////////////////////////////////////////

ITS::ITS( const std::string & its )
 :	itsGraph()
	, itsChargeChange()
	, itsBondChange()
{
	////////////////////  INIT DATA FROM STRING  ////////////////////////

	  // prune leading and tailing whitespaces
	std::string tmp = its.substr(its.find_first_not_of(" \n\t"));
	tmp = tmp.substr(0, (its.find_last_not_of(" \n\t")+1));
	ggl::Graph::vertex_descriptor v, vPrev;
	ggl::Graph::edge_descriptor e;
	boost::property_map< ggl::Graph, ggl::PropNodeLabel >::type vLabel = boost::get( ggl::PropNodeLabel(), itsGraph);
	boost::property_map< ggl::Graph, ggl::PropEdgeLabel >::type eLabel = boost::get( ggl::PropEdgeLabel(), itsGraph);

	int lastBondChange = 0;
	size_t nodeIndex = 0;
	  // successive decomposition of the its string
	while ( ! tmp.empty() ) {
		   // get enclosing brackets
		if (tmp.at(0) != '[') { throw std::runtime_error("missing '['"); }
		size_t end = tmp.find_first_not_of("+-1234567890",1);
		if (end == std::string::npos || tmp.at(end) != ']') {
			throw std::runtime_error("missing ']'");
		}
		  // get charge change
		int chargeChange = atoi( tmp.substr(1, end-1).c_str() );
		  // store charge change
		itsChargeChange.push_back( chargeChange );

		  // get bond
		int bondChange = 0;
		if (end == tmp.size()-1) {throw std::runtime_error("last bond missing"); }
		switch (tmp.at(end+1)) {
			case '-' : bondChange = -1; break;
			case '+' : bondChange = +1; break;
			case '=' : bondChange = 0; break;
			default: throw std::runtime_error("missing bond change label");
		}
		  // store bond change
		itsBondChange.push_back( bondChange );

		  // create charge change node
		v = boost::add_vertex( itsGraph );
		vLabel[v] = boost::lexical_cast< std::string >(chargeChange);

		  // add bond to previous atom
		if (nodeIndex > 0) {
			  // add bond to previous node
			vPrev = boost::vertex( nodeIndex-1, itsGraph );
			e = boost::add_edge( v, vPrev, itsGraph ).first;
			eLabel[e] = lastBondChange == 0 ? "=" : ( lastBondChange < 0 ? "-" : "+" );
		}

		  // prepare next iteration
		nodeIndex++;
		lastBondChange = bondChange;
		tmp = tmp.substr(end+2);
	}
	  // add bond to previous atom
	if (nodeIndex > 0) {
		  // add bond to previous node
		vPrev = boost::vertex( 0, itsGraph );
		e = boost::add_edge( v, vPrev, itsGraph ).first;
		eLabel[e] = lastBondChange == 0 ? "=" : ( lastBondChange < 0 ? "-" : "+" );
	}


	////////////////  GET SYMMETY INFORMATION  ////////////////

	ggl::Graph_boost graph(itsGraph);
	sgm::Pattern pattern( graph );

	  // get order check constraints to break symmetries
	itsOrderConstraints = sgm::PA_OrderCheck::getGraphAutomorphism( pattern ).getCheckList();

	  // get list of symmetric matches
	  // --> calculate all symmetric matches of the graph on itself
	sgm::MR_StoringInsert mr(itsSymmMatch);
	sgm::GM_vf2 gm;
	gm.findMatches( pattern, pattern.getPatternGraph(), mr, UINT_MAX );

	  // remove identity isomorphism
	assert( !itsSymmMatch.empty() );
	sgm::Match identity(itsSymmMatch.begin()->size());
	for (size_t i=0 ; i<identity.size(); ++i) {
		identity[i] = i;
	}
	for (MatchSet::iterator it = itsSymmMatch.begin(); it!= itsSymmMatch.end(); ++it) {
		if (std::equal(identity.begin(), identity.end(), it->begin())) {
			itsSymmMatch.erase(it);
			break;
		}
	}

#ifdef MARTIN
	std::cout <<"\n ITS = \n" <<graph <<"\n";
	for (MatchSet::iterator curSym = itsSymmMatch.begin(); curSym!= itsSymmMatch.end(); ++curSym) {
		std::cout <<"s";
		for (sgm::Match::const_iterator it = curSym->begin(); it != curSym->end(); ++it) {
			std::cout <<"\t" <<*it;
		}
		std::cout <<std::endl;
	}
#endif

}

//////////////////////////////////////////////////////////////////////////

ITS::~ITS()
{
}

//////////////////////////////////////////////////////////////////////////

size_t
ITS::
getSize() const {
	return itsChargeChange.size();
}

//////////////////////////////////////////////////////////////////////////

const ITS::MatchSet &
ITS::
getSymmetricMatches() const {
	return itsSymmMatch;
}

//////////////////////////////////////////////////////////////////////////

const sgm::PA_OrderCheck::CheckList &
ITS::
getNeededOrderChecks() const {
	return itsOrderConstraints;
}

//////////////////////////////////////////////////////////////////////////

const std::vector< int > &
ITS::
getChargeChange() const {
	return itsChargeChange;
}

//////////////////////////////////////////////////////////////////////////

const std::vector< int > &
ITS::
getBondChange() const {
	return itsBondChange;
}


//////////////////////////////////////////////////////////////////////////

const ggl::Graph &
ITS::
getGraph() const {
	return itsGraph;
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////



