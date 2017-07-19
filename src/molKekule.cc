

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <exception>
#include <iterator>

#include <boost/algorithm/string.hpp>

#include "biu/OptionParser.hh"

#include <sgm/GM_vf2.hh>
#include <sgm/MR_Storing.hh>
#include <sgm/SubGraph.hh>
#include "sgm/RP_Hanser96.hh"

#include <ggl/Graph.hh>
#include <ggl/Graph_GML_writer.hh>
#include <ggl/chem/Molecule.hh>
#include <ggl/chem/SMILESwriter.hh>
#include <ggl/chem/SMILESparser.hh>
#include <ggl/chem/AromaticityPerception.hh>




	using namespace ggl;
	using namespace ggl::chem;


//////////////////////////////////////////////////////////////////////////
#ifndef ARGEXCEPTION_
#define ARGEXCEPTION_

	  /*! Exception class for exeptions thrown during argument and input parsing.
	   */
	class ArgException : public std::exception {
	public:
		  //! the error message
		std::string errorMsg;
		ArgException( std::string errorMsg_ ) : errorMsg(errorMsg_) {}
		virtual ~ArgException() throw (){}
		virtual const char* what() const throw() {
			return errorMsg.c_str();
		}
	};

#endif
//////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////


void
initAllowedArguments( biu::OptionMap & allowedArgs, std::string &infoText );

//////////////////////////////////////////////////////////////////////////

using namespace ggl::chem;

/*!
 *
 * Dummy aromaticity perception to access rings to predict
 *
 */
class AP_enumerate
	: public AromaticityPerception
{

protected:

	  //! defines an edge; NOTE: use as ordered pairs, i.e. first <= second
	typedef AromaticityPerception::Edge AromaticEdge;
	  //! defines an aromatic ring via the set of aromatic edges
	typedef AromaticityPerception::EdgeSet AromaticEdgeSet;

public:
	  //! ring describing class
	typedef AromaticityPerception::RingDescriptor RingDescriptor;

protected:
	  //! the aromatic ring container to fill
	using AromaticityPerception::aromaticEdges;


	  //! container that will hold all rings of the current molecule
	  //! -> this list is later pruned to the rings of interest
	using AromaticityPerception::allRings;

	const size_t maxRingSize;

public:

	AP_enumerate(const size_t maxRingSize)
		: maxRingSize(maxRingSize)
	{}

	AP_enumerate( const AP_enumerate & toCopy )
		: AromaticityPerception(toCopy)
		, maxRingSize(toCopy.maxRingSize)
	{}

	virtual ~AP_enumerate()
	{}

	const std::vector< RingDescriptor* > &
	enumerateRings( const Molecule & mol ) {
		  // clear temporary data
		clearData();
		  // identify all rings
		findAllRings( mol );

		  // prune rings that cannot be aromatic
		pruneNonSingleDoubleBondRings( allRings, mol );

		  // prune rings to remove fused representation
		pruneFusedRings( allRings );

		return allRings;
	}

	virtual
	void
	findAllRings( const Molecule & mol ) {

		  // store all rings up to the maximally predictable size
		sgm::RP_Hanser96 ringFinder;
		  // setup graph interface
		Molecule_Graph molGraph(mol);
		  // run ring perception
		ringFinder.findRings( molGraph, *this, maxRingSize );

	}

	void
	identifyAromaticEdges( const Molecule & mol )
	{}


	virtual
	AP_enumerate *
	clone() const
	{
		return new AP_enumerate(*this);
	}




	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
};

	typedef AP_enumerate::RingDescriptor RingDesc;
	typedef std::vector< AP_enumerate::RingDescriptor* > RingVec;
	typedef sgm::RingReporter::RingList RingList;

	typedef AromaticityPerception::Edge Edge;
	typedef AromaticityPerception::EdgeSet EdgeSet;

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

#include <gecode/int.hh>
#include <gecode/search.hh>

class BondEnumerate : public Gecode::Space {

protected:

	  //! array containing one variable for each bond to assign
	Gecode::IntVarArray bonds;

	std::vector< Edge > id2edge;

public:


	  /*!
	   * CSP setup based on ITS, educt, and product information.
	   * @param rings all rings of interest within the molecule
	   * @param mol the molecule the rings are from
	   * @param useBondSum whether or not to use the sum of ring bonds for each
	   *  atom to setup the constraints or to derive the information from the
	   *  valence and all connected non-ring bonds (the latter is needed if
	   *  aromatic edge assignments are given in mol)
	   */
	BondEnumerate(	const RingVec & rings
					, const Molecule & mol
					, const bool useBondSum
				)
		: bonds()
		, id2edge()
	{
		  // data structures
		typedef std::map< size_t, std::set< size_t > > NeighMap;
		NeighMap neighbors;
		std::map< Edge, size_t > edge2id;

		  // needed constraint information
		typedef std::map< size_t, int > IntMap;
		IntMap bondOrderSum;

		  // property maps for access
		boost::property_map<	Molecule , ggl::chem::PropNodeLabel >
			::const_type nodeLabel = boost::get( ggl::chem::PropNodeLabel(), mol );
		boost::property_map<	Molecule , ggl::chem::PropNodeIndex >
			::const_type nodeIndex = boost::get( ggl::chem::PropNodeIndex(), mol );
		boost::property_map<	Molecule , ggl::chem::PropEdgeLabel >
			::const_type edgeLabel = boost::get( ggl::chem::PropEdgeLabel(), mol );

		  // fill bond data structures
		bool aromaticEdgeFound = false;
		for (size_t i=0; i< rings.size(); ++i) {

			const EdgeSet & edges = rings.at(i)->edges;
			for (EdgeSet::const_iterator e = edges.begin(); e!= edges.end(); ++e) {
				if (edge2id.find(*e) != edge2id.end()) {
					continue;
				}
				const size_t newID = id2edge.size();
				id2edge.push_back(*e);
				edge2id[*e] = newID;
				neighbors[ e->first ].insert( e->second );
				neighbors[ e->second ].insert( e->first );

				  // update bond sum information
				if (useBondSum) {
					Molecule::edge_descriptor edge = boost::edge(
										  boost::vertex( e->first, mol )
										, boost::vertex( e->second, mol )
										, mol).first;
					const MoleculeUtil::BondLabelData * const bondData = MoleculeUtil::getBondData(edgeLabel[edge]);
					assert( bondData != NULL );

					  // check for aromatic bonds
					aromaticEdgeFound = bondData->valence != 0;

					  // add bond valence of all ring bonds
					if (bondOrderSum.find(e->first) == bondOrderSum.end()) {
						  // init
						bondOrderSum[e->first] = 0;
					}
					bondOrderSum[e->first] += (int) bondData->valence;
					if (bondOrderSum.find(e->second) == bondOrderSum.end()) {
						  // init
						bondOrderSum[e->second] = 0;
					}
					bondOrderSum[e->second] += (int) bondData->valence;
				}
			}

		}

		if ( useBondSum && aromaticEdgeFound) {
			throw std::runtime_error("bond sum to be use but aromatic ring bond found. please rerun without bond sum parameter.");
		}

		  // setup variables
		bonds = Gecode::IntVarArray(*this, (int)id2edge.size(), 1, 2);

		if (!useBondSum) {
			  // get available bond order sum for each ring atom
			for (NeighMap::const_iterator n = neighbors.begin(); n!= neighbors.end(); ++n) {
				Molecule::vertex_descriptor node = boost::vertex( n->first, mol );
				const MoleculeUtil::AtomLabelData* const atomData = MoleculeUtil::getAtomData( nodeLabel[node] );
				assert(atomData != NULL);

				  // init with overall available valence for bonds
				bondOrderSum[n->first] = (int)atomData->valence;

				  // take charge into account (take care if current atom is a proton!)
				const int curCharge = (int)MoleculeUtil::getCharge( nodeLabel[node] );
				bondOrderSum[n->first] += (MoleculeUtil::getAtom(nodeLabel[node])=="H"?-1:+1) * curCharge;

				  // check if an explicit maximal proton number was given
				if (atomData->valence != atomData->maxProtons) {
					  // reduce available valence
					  // TODO : might not be correct
					bondOrderSum[n->first] += (int)atomData->maxProtons - (int)atomData->valence;
				}

				  // subtract bond valence of all non-ring edges
				Molecule::out_edge_iterator eIt, eItEnd;
				for (boost::tie(eIt,eItEnd)= boost::out_edges(node, mol); eIt != eItEnd; ++eIt) {
					  // check if target node is unknown neighbor
					if (n->second.find(nodeIndex[boost::target(*eIt,mol)]) == n->second.end()) {
						const MoleculeUtil::BondLabelData * const bondData = MoleculeUtil::getBondData(edgeLabel[*eIt]);
						assert( bondData != NULL );
						  // substract bond valence of non-ring bond
						bondOrderSum[n->first] -= (int) bondData->valence;
						assert(bondData->isAromatic == 0);
					}
				}
				  // check if enough valence left for remaining ring bonds
				assert( bondOrderSum[n->first] >= (int)n->second.size() );
			}
		}

		  // setup constraints
		  // get available bond order sum for each ring atom
		for (NeighMap::const_iterator n = neighbors.begin(); n!= neighbors.end(); ++n) {
			  // bond set for this node to fill
			Gecode::IntVarArgs adjBonds;
			for (std::set<size_t>::const_iterator neigh = n->second.begin(); neigh != n->second.end(); ++neigh) {
				Edge edge = Edge(n->first, *neigh);
				if (edge.first > edge.second) {
					edge = Edge(*neigh, n->first);
				}
				assert( edge2id.find(edge) != edge2id.end() ); // has to be known edge
				  // add according variable to bond set
				adjBonds << bonds[(int) edge2id.find(edge)->second ];
			}
			assert(adjBonds.size() == (int)n->second.size());
			Gecode::linear(*this, adjBonds, Gecode::IRT_EQ, bondOrderSum[n->first], Gecode::ICL_DOM);
		}

	}


	virtual ~BondEnumerate()
	{}

	  /*!
	   * Copy constructor
	   * @param bool share states whether a shared copy is constructed
	   * @param toCopy is a reference to an object of the this class
	   */
	BondEnumerate( bool share, BondEnumerate & toCopy)
		: Space( share, toCopy )
		, bonds( toCopy.bonds)
		, id2edge( toCopy.id2edge )
	{
		bonds.update( *this, share, toCopy.bonds );
	}

	  /*!
	   * The function is responsible for space cloning. It calls the previously
	   * mentioned copy constructor.
	   * @param bool share states whether a shared copy is constructed
	   * @return a pointer to an updated copy of the underlying space
	   */
	virtual
	Space*
	copy( bool share )
	{
		  // call copy constructor
		return new BondEnumerate( share, *this );
	}

	Gecode::IntVarArray&
	getBonds()
	{
		return bonds;
	}

	size_t
	numOfAmbiguousEdges() const
	{
		  // get number of unassigned bond variables
		size_t unassigned = 0;
		for (int i = 0; i<bonds.size(); ++i) {
			if (!bonds[i].assigned()) {
				++unassigned;
			}
		}
		return unassigned;
	}

	size_t
	numOfAmbiguousRings( const RingVec & rings) const
	{
		  // get set of unassigned edges
		EdgeSet ambiEdges;
		for (int i = 0; i<bonds.size(); ++i) {
			if (!bonds[i].assigned()) {
				ambiEdges.insert(id2edge.at((size_t)i));
			}
		}
		  // find all rings overlapping the set of unassigned edges
		size_t ambiRings = 0;
		for (size_t i=0; i< rings.size(); ++i) {
			const EdgeSet & edges = rings.at(i)->edges;
			EdgeSet::const_iterator amiIt = ambiEdges.begin(), curIt = edges.begin();
			bool noOverlap = true;
			while( noOverlap && amiIt != ambiEdges.end() && curIt != edges.end()) {
				noOverlap = *amiIt != *curIt;
				if (noOverlap) {
					if ( *amiIt < *curIt ) {
						++amiIt;
					} else {
						++curIt;
					}
				}
			}
			  // check if overlap with ambiguous edges
			if (!noOverlap) {
				++ambiRings;
			}
		}
		  // return number of rings overlapping with unassigned edges
		return ambiRings;
	}


};



//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


int main( int argc, char** argv ) {

	//////////////////////////////////////////////////////////////
	// variables
	//////////////////////////////////////////////////////////////

	std::ostream* out = &std::cout;
	std::ofstream* outFile = NULL;

	//////////////////////////////////////////////////////////////
	// parameter parsing and checking
	//////////////////////////////////////////////////////////////

	biu::OptionMap allowedArgs;	//< map of allowed arguments
	std::string infoText;		//< info string of the program
	initAllowedArguments(allowedArgs,infoText);	// init

		// parse programm arguments
	biu::COptionParser opts = biu::COptionParser(	allowedArgs, argc,
													argv, infoText);
		// check arguments parseable and all mandatory arguments given
		// + help output if needed
	if (opts.argExist("help")) {
		opts.coutUsage();
		return 0;
	}
	if (!opts.noErrors()) {
		return -1;
	}

	int exitValue = 0;

	try {



		  // set SMILES input stream
		const std::string molSMILES = opts.getStrVal("smiles");
		if (molSMILES.empty()) {
			throw ArgException("no SMILES given");
		}
		  // set output stream
		if (opts.getStrVal("out").size() == 0) {
			throw ArgException("no output file given");
		} else if ( !boost::iequals(opts.getStrVal("out"),"STDOUT")) {
			outFile = new std::ofstream(	opts.getStrVal("out").c_str()
											, std::ofstream::out );
			if (!outFile->is_open()) {
				std::ostringstream oss;
				oss	<<"cannot open output file '" <<opts.getStrVal("out") <<"'";
				throw ArgException(oss.str());
			}
			out = outFile;
		}

		const int maxRingSize = opts.getIntVal("maxRingSize");
		if (maxRingSize < 3) {
			throw ArgException("'maxRingSize' has to be >= 3");
		}

		//////////////////////////////////////////////////////////////
		// parse molecules from input
		//////////////////////////////////////////////////////////////

		  // convert to molecule
		std::pair<Molecule,int> ret = ggl::chem::SMILESparser::parseSMILES( molSMILES );
		if (ret.second >= 0) {
			throw ArgException("cannot parse SMILES");
		}
		Molecule mol = ret.first;

		  // fill protons
		MoleculeUtil::fillProtons( mol );


		//////////////////////   RING ENUMERATION  /////////////////////////////////////

		AP_enumerate ringEnumerator(maxRingSize);
		const RingVec & allRings = ringEnumerator.enumerateRings( mol );

		boost::property_map<	Molecule , ggl::chem::PropEdgeLabel >
			::type edgeLabel = boost::get( ggl::chem::PropEdgeLabel(), mol );

		 // get number of aromatic rings
		size_t aromRings = 0;
		  // fill bond data structures
		for (size_t i=0; i< allRings.size(); ++i) {
			const EdgeSet & edges = allRings.at(i)->edges;
			bool allAromatic = true;
			for (EdgeSet::const_iterator e = edges.begin(); allAromatic && e!= edges.end(); ++e) {
				const MoleculeUtil::BondLabelData * const bondData
					= MoleculeUtil::getBondData( edgeLabel[
					      boost::edge(
					    		  boost::vertex( e->first, mol),
					    		  boost::vertex( e->second, mol),
					    		  mol ).first ] );
				allAromatic = bondData->isAromatic != 0;
			}
			  // check if all bond labels are aromatic
			if (allAromatic) {
				++aromRings;
			}
		}

		 // get number of overlapping rings
		size_t overlappingRings = 0;
		  // fill bond data structures
		for (size_t i=0; i< allRings.size(); ++i) {
			const EdgeSet & ri = allRings.at(i)->edges;
			for (size_t j=i+1; j< allRings.size(); ++j) {
				const EdgeSet & rj = allRings.at(j)->edges;
				EdgeSet::const_iterator iIt = ri.begin(), jIt = rj.begin();
				bool noOverlap = true;
				while( noOverlap && iIt != ri.end() && jIt != rj.end()) {
					noOverlap = *iIt != *jIt;
					if (noOverlap) {
						if ( *iIt < *jIt ) {
							++iIt;
						} else {
							++jIt;
						}
					}
				}
				if (!noOverlap) {
					++overlappingRings;
				}
			}
		}

		//////////////////////   MESOMER ENUMERATION  /////////////////////////////////////


		BondEnumerate bondEnumerate(allRings, mol, opts.argExist("useBondSum"));

		Gecode::StatusStatistics stats;
		bondEnumerate.status(stats);

		if (bondEnumerate.failed()) {
			std::cerr <<"CSP failed\n";
		}

		//////////////////////   SEARCH & OUTPUT  ////////////////////////////////

		std::cout <<" aromRings " <<aromRings;
		std::cout <<" ambiRings " <<bondEnumerate.numOfAmbiguousRings(allRings) <<" of " <<allRings.size();
		std::cout <<" ringOverlaps " <<overlappingRings;
		std::cout <<" ambiBonds " <<bondEnumerate.numOfAmbiguousEdges() <<" of " <<bondEnumerate.getBonds().size();

		  // branching on products
		Gecode::branch(bondEnumerate, bondEnumerate.getBonds(), Gecode::INT_VAR_SIZE_MIN(), Gecode::INT_VAL_MIN());

		size_t numOfMesomers = 0;
		Gecode::DFS<BondEnumerate> dfsITS( &bondEnumerate );
		while ( BondEnumerate* solution = dfsITS.next() ) {
			++numOfMesomers;
			delete solution;
		}

		std::cout <<" mesomers " <<numOfMesomers <<std::endl;



	} catch (std::exception& ex) {
		std::cerr <<"\n\n ERROR : " <<ex.what() <<"\n"<<std::endl;
		exitValue = -1;
	}


	//////////////////////////////////////////////////////////////
	// final stream handling
	//////////////////////////////////////////////////////////////

	out = &std::cout;
	if (outFile != NULL)	{ outFile->close(); delete outFile; }

	return exitValue;
}



//////////////////////////////////////////////////////////////////////////

/*!
 * Initialises allowed parameters and their default values and prepares an
 * information string for the tool.
 */
void
initAllowedArguments(biu::OptionMap & allowedArgs, std::string &infoText )
{
	infoText = "\n"
		"Enumerates all Kelkule structures for a given molecule.\n"
		;

	allowedArgs.push_back(biu::COption(
							"smiles", false, biu::COption::STRING,
							"The SMILES of the molecule to enumerate",
							"STDIN"));
	allowedArgs.push_back(biu::COption(
							"maxRingSize", true, biu::COption::INT,
							"The maximal ring size to be considered",
							"15"));
	allowedArgs.push_back(biu::COption(
							"useBondSunm", true, biu::COption::BOOL,
							"If given, the sum of ring bond valences is used for constraint setup. Otherwise this information is derived from the atom valence etc (required if aromatic edges are annotated).",
							"15"));
	allowedArgs.push_back(biu::COption(
							"out", true, biu::COption::STRING,
							"Output file name or 'STDOUT' when to write to standard output",
							"STDOUT"));
	allowedArgs.push_back(biu::COption(
							"help", true, biu::COption::BOOL,
							"Displays help on all parameters"));
}

//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
