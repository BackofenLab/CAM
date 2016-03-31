
#include <iostream>
#include <fstream>
#include <exception>

#include "biu/OptionParser.hh"
#include "biu/Timer.hh"

#include "ReactionSMILES_Extractor.h"
#include "ITS_Graph_Interface.h"
#include "ITS_CSP.h"
#include <gecode/search.hh>
#include <ggl/chem/SMILESwriter.hh>

#include <boost/lexical_cast.hpp>

////////////////////////////////////////////////////////////////////////////////////

class AtomMappingStatistics {

public:

	size_t solutionNumber;
	size_t validSolNumber;

	  //! Execution statistics
	Gecode::Search::Statistics stat;

	  // time consumption
	double timeCSP;
	double timeVF2;
	double timeAll;
	double time1st;

	  /*!
	   * Timers to calculate the time consumption for each
	   * procedure and then store timings in the above variables.
	   */
	Timer CSPtimer,
		  VF2timer,
		  allTimer,
		  firstTimer;

	AtomMappingStatistics()
	:
		solutionNumber(0)
		, validSolNumber(0)
		, timeCSP(0.0)
		, timeVF2(0.0)
		, timeAll(0.0)
		, time1st(0.0)
	{

	}
	//////////////////////////////////////////

	void
	update( const Gecode::Search::Statistics _stat )
	{
		stat = _stat;
	}
	//////////////////////////////////////////

	void
	incSolution( )
	{
		++solutionNumber;
	}
	//////////////////////////////////////////

	void
	incValidSol( )
	{
		++validSolNumber;
	}
	//////////////////////////////////////////

	void
	onTimerCSP()
	{
	    CSPtimer.start();
	}
	//////////////////////////////////////////

	void
	offTimerCSP()
	{
		timeCSP += CSPtimer.stop();
	}
	//////////////////////////////////////////

	void
	onTimerVF2()
	{
	    VF2timer.start();
	}
	//////////////////////////////////////////

	void
	offTimerVF2()
	{
		timeVF2 += VF2timer.stop();
	}
	//////////////////////////////////////////

	void
	onTimerAll()
	{
	    allTimer.start();
	}
	//////////////////////////////////////////

	void
	offTimerAll()
	{
		timeAll += allTimer.stop();
	}
	//////////////////////////////////////////

	void
	onTimer1st()
	{
	    firstTimer.start();
	}
	//////////////////////////////////////////

	void
	offTimer1st()
	{
		time1st += firstTimer.stop();
	}
	//////////////////////////////////////////

	void
	print() const {
		std::cout	<<"\n Statistics:\n"
			<<"\t number of ITS candidates = " << solutionNumber << "\n"
			<< "\t valid atom mappings = " << validSolNumber << "\n"
//			<< "\t propagations = " << stat.propagate << "\n"
//			<< "\t failures = " << stat.fail << "\n"
//			<< "\t nodes = " << stat.node << "\n"
//			<< "\t peak depth = " << stat.depth << "\n"
//			<< "\t peak memory = " << static_cast<int>((stat.memory+1023) / 1024) << " KB" << "\n"
			<< "\t timings: "<< "\n"
			<< "\t first solution = " << time1st<< " millisec" << "\n"
			<< "\t CSP = " << timeCSP << " millisec" << "\n"
			<< "\t VF2 = " << timeVF2 << " millisec" << "\n"
			<< "\t overall = " << timeAll << " millisec" << "\n";;
	}

};

////////////////////////////////////////////////////////////////////////////////////

/*!
 * Initialises allowed parameters and their default values and prepares an
 * information string for the tool.
 */
void
initAllowedArguments(biu::OptionMap & allowedArgs, std::string &infoText )
{
	infoText = "\n"
			"CAM - Constraint-based Atom-Atom Mapping - identifies atom maps that encode chemically feasible cyclic Imaginary Transition States (ITS).\n"
			"Beside using ITS derived from known reaction mechanisms, the user can provide new ITS encodings.\n"
		;

	allowedArgs.push_back(biu::COption(
							"reaction", false, biu::COption::STRING,
							 "reaction encoding in SMIRKS format, i.e. dot-separated educt SMILES followed by '>>' and the dot-separated product SMILES"));
	allowedArgs.push_back(biu::COption(
							"k", true, biu::COption::INT,
							"ITS ring size. If not specified all ring sizes are tested in increasing order until a match was found."));
	allowedArgs.push_back(biu::COption(
							"ITS", true, biu::COption::STRING,
							"The ITS layout string encoding to use. If not given, all ITS layouts are tested in increasing size."));;
	allowedArgs.push_back(biu::COption(
							"solPerITS", true, biu::COption::INT,
							"The number of atom mappings to enumerate for each ITS.",
							"99999"));
	allowedArgs.push_back(biu::COption(
							"allITS", true, biu::COption::BOOL,
							"Enumerate atom mappings for all known ITS"));
	allowedArgs.push_back(biu::COption(
							"noProtonFilling", true, biu::COption::BOOL,
							"Disable the adding of protons not part of the input reaction SMILES"));
	allowedArgs.push_back(biu::COption(
							"allowSubset", true, biu::COption::BOOL,
							"If given, the ITS mapping is allowed to cover only a subset of the given molecules. Otherwise, each input/output molecule has to be part of the ITS."));
	allowedArgs.push_back(biu::COption(
								"printGraphs", true, biu::COption::BOOL,
								"Display relabeled ITS graphs of educts and products"));
	allowedArgs.push_back(biu::COption(
							"v", true, biu::COption::BOOL,
							"Enable verbose output"));
}

////////////////////////////////////////////////////////////////////////////////////

 /*!
  * Checks whether or not a given ITS is compatible to the educts and products
  * participating in the ITS.
  * @param its the ITS of interest
  * @param educts the educt graph representation
  * @param products the product graph representation
  */
bool
areCompatible( const ITS &its
				, const ReactionGraph & educts
				, const ReactionGraph & products
) {

	  // if ITS contains charge change
	std::vector< int >::const_iterator it;
	for ( it = its.getChargeChange().begin(); it != its.getChargeChange().end(); ++it ) {
		if ( *it != 0 ) {
			if ( !educts.isHomovalent() )
				return true;

			if ( !products.isHomovalent() )
				return true;

			return false;
		}
	}
	return true;
}

////////////////////////////////////////////////////////////////////////////////////

const std::string
annotateITS( const sgm::Graph_Interface & graphITS
		   , const sgm::Match itsIndices
		   , const ITS & its )
{
	 ggl::chem::Molecule itsMolecule;
	 ggl::chem::MoleculeUtil::copy( graphITS, itsMolecule );
	 boost::property_map< ggl::chem::Molecule, ggl::chem::PropNodeLabel >::type
	 	 	 itsAtomLabel = boost::get( ggl::chem::PropNodeLabel(), itsMolecule );
	 std::string itsString = "";
	 for ( size_t i = 0; i < itsIndices.size(); i++ )
	 {
		 std::string label = ggl::chem::MoleculeUtil::getAtom(itsAtomLabel[boost::vertex(itsIndices.at(i), itsMolecule)]);
		 itsString += "[" + label;
		 if (its.getChargeChange().at(i) != 0) {
			 if (its.getChargeChange().at(i)>0) {
				 itsString += "+";
			 }
			 itsString +=  boost::lexical_cast< std::string >( its.getChargeChange().at(i) );
		 }
		 itsString += ":" + boost::lexical_cast< std::string >( itsIndices.at(i) + 1 ) + "]"
				 ;
		 switch( its.getBondChange().at(i) ) {
		 	 case  0 : itsString += "="; break;
		 	 case +1 : itsString += "+"; break;
		 	 case -1 : itsString += "-"; break;
		 	 default :
		 		 throw std::runtime_error("unexpected bond change");
		 }
	 }
	 return itsString;
}

////////////////////////////////////////////////////////////////////////////////////

class AtomMapping_Reporter : public sgm::Match_Reporter {

	  //! string encoding of the ITS ring only
	const std::string itsStr;

	  //! educts, products original graphs ( without relabeled edges )
	ggl::chem::Molecule_Graph eduGraph, proGraph;

	  //! symmetrical solutions set of products
	const ReactionGraph::SymmetricMatchSet proSymmSet;

	  //! list of already known matching to avoid duplicates
	std::set< std::string > & knownMappings;

public:

	AtomMapping_Reporter(const std::string _itsStr
			, const ggl::chem::Molecule_Graph _eduGraph
			, const ggl::chem::Molecule_Graph _proGraph
			, const ReactionGraph::SymmetricMatchSet _proSymmSet
			, std::set< std::string > & _knownMappings)
	:
		itsStr(_itsStr)
		, eduGraph(_eduGraph)
		, proGraph(_proGraph)
	    , proSymmSet(_proSymmSet)
		, knownMappings(_knownMappings)
    {

    }


	sgm::Match
	getSmallestSymmetry( const sgm::Match& curSol ) const
	{
		  // set of all unique symmetric solutions
		ReactionGraph::SymmetricMatchSet symSol;
		symSol.insert( curSol );

		  // generate all value symmetric solutions and store
		sgm::Match nextSym(curSol.size());
		for (ReactionGraph::SymmetricMatchSet::const_iterator
				curSym = proSymmSet.begin();
				curSym != proSymmSet.end(); ++curSym)
		{
			  // get symmetric match, ie. symmetric matching positions of current assignment
			for (size_t i=0; i<curSol.size(); ++i) {
				assert( curSol.at(i) < curSym->size() );
				nextSym[i] = curSym->at( curSol.at(i) );
			}
			  // insert to set (ensures uniqueness of symmetric solutions)
			symSol.insert(nextSym);
		}

		return *symSol.begin();
	}




	void
	reportHit ( const sgm::Pattern_Interface & eduPattern
			, const sgm::Graph_Interface & targetGraph
			, const sgm::Match & match )
	{
		  // get smallest symmetry of this match to avoid symmetric solutions
		sgm::Match nonSymmetricMatch = getSmallestSymmetry( match );

		const sgm::SubGraph mappedProGraph( proGraph,  nonSymmetricMatch );
		// const sgm::SubGraph mappedProGraph( proGraph,  match )  ;

		  // copy original educt and product graphs
		ggl::chem::Molecule educts, products;
		ggl::chem::MoleculeUtil::copy( eduGraph, educts );
		ggl::chem::MoleculeUtil::copy( mappedProGraph, products );
		assert(boost::num_vertices(educts) == boost::num_vertices(products));
		  // add IDs to atom labels
		  // retrieving label informations for every atom in the molecule
		boost::property_map< ggl::chem::Molecule, ggl::chem::PropNodeLabel >::type
						eLabel = boost::get( ggl::chem::PropNodeLabel(), educts),
						pLabel = boost::get( ggl::chem::PropNodeLabel(), products);
		const size_t numAtoms = boost::num_vertices(educts);
		for ( size_t i = 0; i < numAtoms; ++i ) {
		      // set educt class information
			std::string label = eLabel[ boost::vertex(i, educts) ];
		      // access to the class information for every atom label,
		      // no class must be present
			assert(ggl::chem::MoleculeUtil::getClass(label) == 0);
			  // adding (:number) to atoms labels
			eLabel[boost::vertex(i, educts)] = label + ":" + boost::lexical_cast< std::string >( i + 1 ) ;

		      // set product class information
			label = pLabel[ boost::vertex(i, products) ];
		      // access to the class information for every atom label,
			  // no class must be present
			assert(ggl::chem::MoleculeUtil::getClass(label) == 0);
			  // adding (:number) to atoms labels
			pLabel[boost::vertex(i, products)] = label + ":" + boost::lexical_cast< std::string >( i + 1 ) ;
		}

		  /*!
		   * checking the case whether the educt resp. products contains
		   * multiple molecules. Returning the component label accordingly
		   * i.e. in which molecule each component locates.
		   */
		sgm::Graph_Interface::CompLabel eduComplabel, proComplabel;
		size_t connectedEdu, connectedPro;
		connectedEdu = sgm::Graph_Interface::connectedComponents( eduGraph, eduComplabel );
		connectedPro = sgm::Graph_Interface::connectedComponents( mappedProGraph, proComplabel );
		std::string mappedEduStr = "", mappedProStr = "";

		std::multiset< std::string > sortedSMILES;
	      // educts: create SMILES for each component
		for ( size_t c = 0; c < connectedEdu; c++ ) {
			  // collect node IDs of this component
			sgm::SubGraph::NodeList compNodes;
			for( size_t i = 0; i < numAtoms; ++i ) {
				if ( eduComplabel.at(i) == c ) {
				   compNodes.push_back(i);
				}
			}
			  // create molecule
			ggl::chem::Molecule nextMol;
			ggl::chem::MoleculeUtil::copy(sgm::SubGraph(ggl::chem::Molecule_Graph(educts),compNodes), nextMol);
			  // append SMILES
			sortedSMILES.insert(ggl::chem::SMILESwriter::getSMILES(nextMol, true));
		}

		  // chain all sorted SMILES to one multi SMILES
		std::multiset< std::string >::const_iterator s = sortedSMILES.begin();
		assert(s != sortedSMILES.end());
		mappedEduStr = *s;
		++s;
		while(s != sortedSMILES.end()) {
			mappedEduStr += "." + *s;
			++s;
		}

		sortedSMILES.clear();
		  // products: create SMILES for each component
		for (size_t c = 0; c < connectedPro; c++) {
		     // collect node IDs of this component
		   sgm::SubGraph::NodeList compNodes;
		   for( size_t i  =0; i < numAtoms; ++i) {
			  if (proComplabel.at(i) == c) {
				 compNodes.push_back(i);
			  }
		   }
		     // create molecule
		   ggl::chem::Molecule nextMol;
		   ggl::chem::MoleculeUtil::copy(sgm::SubGraph(ggl::chem::Molecule_Graph(products),compNodes), nextMol);
			 // append SMILES
		   sortedSMILES.insert(ggl::chem::SMILESwriter::getSMILES(nextMol, true));
		}
		  // chain all sorted SMILES to one multi SMILES
		s = sortedSMILES.begin();
		assert(s != sortedSMILES.end());
		mappedProStr = *s;
		++s;
		while(s != sortedSMILES.end()) {
		  mappedProStr += "." + *s;
		  ++s;
		}


		  // final mapped SMILES string
		const std::string finalMapping = mappedEduStr + mappedProStr;

		  // check if matching is not already known
		if (knownMappings.find(finalMapping) == knownMappings.end()) {
			  // store
			knownMappings.insert(finalMapping);
			  // print
			std::cout << mappedEduStr + ">" + itsStr + ">" + mappedProStr << std::endl;
			std::cout << std::endl;
		}
	}
};

////////////////////////////////////////////////////////////////////////////////////

bool
findMappings( ITS_CSP & csp
			, AtomMappingStatistics&  statistics
			, const bool printGraphs
			, const size_t maxSolutionNumber
) {


	std::set< std::string > finalMappings;
	Gecode::DFS<ITS_CSP> dfsITS( &csp );
	double firstSol = true;
	statistics.onTimerAll();
	statistics.onTimerCSP();
	statistics.onTimer1st();
	const size_t initialSolNumber = statistics.validSolNumber;
	ITS_CSP* solution = NULL;
	while ( ((statistics.validSolNumber-initialSolNumber) < maxSolutionNumber)
			&& (solution = dfsITS.next()) )
	{
	    statistics.offTimerCSP();
	  	statistics.incSolution();
		// solution->print();
		sgm::Match eduNewITSIndices, proNewITSIndices;
		ggl::chem::Molecule eduMol, proMol;
		eduMol = solution->getEduGraph().getGraphToMatch( solution->getEduITSsolution(), eduNewITSIndices );
		proMol = solution->getProGraph().getGraphToMatch( solution->getProITSsolution(), proNewITSIndices );
		ggl::chem::Molecule_Graph eduMolGraph(eduMol), proMolGraph(proMol);

		ITS_Graph_Interface eduGraphITS( eduMolGraph, eduNewITSIndices );
		ITS_Graph_Interface proGraphITS( proMolGraph, proNewITSIndices );

		if ( printGraphs ) {
			std::cout << "Educt ITS:"   << eduGraphITS << std::endl;;
			std::cout << "Product ITS:" << proGraphITS << std::endl;
			std::cout << std::endl;
		}

		sgm::GM_vf2 matcher;

		  // get symmetries of current shrinked product graph to match and store them
		ReactionGraph::SymmetricMatchSet proMolGraphSymmetries;
		sgm::MR_StoringInsertT< ReactionGraph::SymmetricMatchSet > proMolGraphSymStore( proMolGraphSymmetries );
		matcher.findMatches( sgm::Pattern( proMolGraph ), proMolGraph, proMolGraphSymStore, UINT_MAX );

		std::string itsEncoding = annotateITS( eduGraphITS, eduNewITSIndices, solution->getITS() );
		AtomMapping_Reporter reporter( itsEncoding, eduMolGraph, proMolGraph,
				proMolGraphSymmetries, finalMappings );
		const size_t matchesSoFar = finalMappings.size();
		statistics.onTimerVF2(); // VF-2 Timer starts
		matcher.findMatches( sgm::Pattern( eduGraphITS ), proGraphITS, reporter, 1);
		if ( matchesSoFar < finalMappings.size() )
		{
			if ( firstSol ) {
				statistics.offTimer1st();
				firstSol = false;
			}
			statistics.incValidSol();
		}
		statistics.offTimerVF2(); // VF-2 Timer stops

	    delete solution;
	    statistics.onTimerCSP();
	}// while
	statistics.offTimerCSP();
	statistics.offTimerAll();
	statistics.update( dfsITS.statistics() );


	return initialSolNumber != statistics.validSolNumber;
}

////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]) {

	biu::OptionMap allowedArgs;	//< map of allowed arguments
	std::string infoText;		//< info string of the program
	initAllowedArguments(allowedArgs, infoText);	// init

	size_t minK = 3;
	size_t maxK = 8;
	size_t maxSolPerITS = 0;

	  // list of known ITS layouts in string encoding
	  // where the end of the list is indicated by an empty string
	const std::string knownITSLayouts[] = {
		// length 3
		"[+2]+[0]-[0]+",
		"[+1]+[0]-[-1]=",
		// length 4
		"[0]+[0]-[0]+[0]-",
		"[+1]+[0]-[0]+[+1]=",
		"[-1]-[0]+[0]-[+1]=",
//		"graph [ node [ id 2 label \"0\"] node [ id 1 label \"0\"] node [ id 22 label \"0\"] node [ id 23 label \"0\"] edge [ source 1 target 23 label \"+2\" ] edge [ source 22 target 23 label \"-2\" ] edge [ source 2 target 22 label \"+2\" ] edge [ source 1 target 2 label \"-2\" ] ]",
//		"graph [ node [ id 2 label \"0\"] node [ id 1 label \"0\"] node [ id 22 label \"0\"] node [ id 23 label \"0\"] edge [ source 1 target 23 label \"-2\" ] edge [ source 22 target 23 label \"+2\" ] edge [ source 2 target 22 label \"-2\" ] edge [ source 1 target 2 label \"+2\" ] ]",
		// length 5
		"[+1]+[0]-[0]+[0]-[-1]=",
		"[+2]+[0]-[0]+[0]-[0]+",
//		"graph [ node [ id 3 label \"0\"] node [ id 2 label \"0\"] node [ id 4 label \"0\"] node [ id 10 label \"0\"] node [ id 5 label \"0\"] edge [ source 4 target 5 label \"-\" ] edge [ source 3 target 10 label \"+\" ] edge [ source 3 target 4 label \"+\" ] edge [ source 2 target 3 label \"-2\" ] edge [ source 5 target 10 label \"-\" ] edge [ source 2 target 5 label \"+2\" ] ]",
		// length 6
		"[0]+[0]-[0]+[0]-[0]+[0]-",
		"[+1]+[0]-[0]+[0]-[0]+[+1]=",
		"[-1]-[0]+[0]-[0]+[0]-[-1]=",
//		"graph [ node [ id 17 label \"0\"] node [ id 1 label \"0\"] node [ id 18 label \"0\"] node [ id 9 label \"0\"] node [ id 14 label \"0\"] node [ id 19 label \"0\"] edge [ source 1 target 19 label \"-\" ] edge [ source 1 target 18 label \"-\" ] edge [ source 19 target 17 label \"+\" ] edge [ source 14 target 17 label \"-\" ] edge [ source 9 target 14 label \"-\" ] edge [ source 18 target 9 label \"+\" ] edge [ source 1 target 14 label \"+2\" ] ]",
//		"graph [ node [ id 17 label \"0\"] node [ id 1 label \"0\"] node [ id 18 label \"0\"] node [ id 9 label \"0\"] node [ id 14 label \"0\"] node [ id 19 label \"0\"] edge [ source 1 target 19 label \"+\" ] edge [ source 1 target 18 label \"+\" ] edge [ source 19 target 17 label \"-\" ] edge [ source 14 target 17 label \"+\" ] edge [ source 9 target 14 label \"+\" ] edge [ source 18 target 9 label \"-\" ] edge [ source 1 target 14 label \"-2\" ] ]",
		// length 7
		"[+1]+[0]-[0]+[0]-[0]+[0]-[-1]=",
		"[+2]+[0]-[0]+[0]-[0]+[0]-[0]+",
//		"graph [ node [ id 17 label \"0\"] node [ id 2 label \"0\"] node [ id 1 label \"0\"] node [ id 59 label \"0\"] node [ id 42 label \"0\"] node [ id 12 label \"0\"] node [ id 47 label \"0\"] edge [ source 1 target 12 label \"+2\" ] edge [ source 2 target 12 label \"-\" ] edge [ source 47 target 59 label \"-\" ] edge [ source 2 target 42 label \"+\" ] edge [ source 12 target 17 label \"-\" ] edge [ source 42 target 47 label \"+\" ] edge [ source 59 target 17 label \"+\" ] edge [ source 1 target 42 label \"-2\" ] ]",
		// length 8
		"[0]+[0]-[0]+[0]-[0]+[0]-[0]+[0]-",
		"[+1]+[0]-[0]+[0]-[0]+[0]-[0]+[+1]=",
		"[-1]-[0]+[0]-[0]+[0]-[0]+[0]-[-1]=",
//		"graph [ node [ id 2 label \"0\"] node [ id 1 label \"0\"] node [ id 30 label \"0\"] node [ id 27 label \"0\"] node [ id 31 label \"0\"] node [ id 39 label \"0\"] node [ id 38 label \"0\"] node [ id 34 label \"0\"] edge [ source 2 target 38 label \"-\" ] edge [ source 1 target 34 label \"-2\" ] edge [ source 39 target 27 label \"+\" ] edge [ source 1 target 30 label \"+\" ] edge [ source 2 target 34 label \"+\" ] edge [ source 31 target 39 label \"-\" ] edge [ source 38 target 1 label \"+\" ] edge [ source 27 target 30 label \"-\" ] edge [ source 31 target 34 label \"+\" ] ]",
//		"graph [ node [ id 2 label \"0\"] node [ id 1 label \"0\"] node [ id 30 label \"0\"] node [ id 27 label \"0\"] node [ id 31 label \"0\"] node [ id 39 label \"0\"] node [ id 38 label \"0\"] node [ id 34 label \"0\"] edge [ source 2 target 38 label \"+\" ] edge [ source 1 target 34 label \"+2\" ] edge [ source 39 target 27 label \"-\" ] edge [ source 1 target 30 label \"-\" ] edge [ source 2 target 34 label \"-\" ] edge [ source 31 target 39 label \"+\" ] edge [ source 38 target 1 label \"-\" ] edge [ source 27 target 30 label \"+\" ] edge [ source 31 target 34 label \"-\" ] ]",
		// end of list
		""
	};


		// parse programm arguments
	biu::COptionParser opts = biu::COptionParser(	allowedArgs, argc,
													argv, infoText);
	if (opts.noErrors()) {
		  // rules are officially not mandatory.. we have to request manually
		if (opts.getStrVal("reaction").size() == 0) {
			std::cerr <<"\n\n\tERROR : no reaction given\n\n";
			return -1;
		}
		  // check ITS size if given
		if (opts.argExist("k") && opts.getIntVal("k") < 3 ) {
			std::cerr <<"\n\n\tERROR : k has to be at least 3\n\n";
			return -1;
		}
		  // ensure either ITS size or layout given
		if (opts.argExist("k") && opts.argExist("ITS")) {
			std::cerr <<"\n\n\tERROR : either provide an ITS ring size or a layout; both is not possible\n\n";
			return -1;
		}
		  // ensure either ITS size or layout given
		if (opts.argExist("k") && opts.argExist("allITS")) {
			std::cerr <<"\n\n\tERROR : either provide an ITS ring size or request the enumeration of all ITS layouts; both is not possible\n\n";
			return -1;
		}
		if (opts.argExist("ITS")) {
			try {
				  // test ITS for parsability
				ITS testITS(opts.getStrVal("ITS"));
				  // test ITS size
				if (testITS.getSize() < 3) {
					std::cerr <<"\n\n\tERROR : given ITS layout has less than 3 nodes\n\n";
					return -1;
				}
			} catch (std::exception & ex) {
				std::cerr <<"\n\n\tERROR : cannot parse ITS layout:\n" <<ex.what() <<"\n";
				return -1;
			}
		}

		if (opts.getIntVal("solPerITS") < 1) {
			std::cerr <<"\n\n\tERROR : 'solPerITS' has to be >= 1\n";
			return -1;
		} else {
			maxSolPerITS = (size_t)opts.getIntVal("solPerITS");
		}

	} else {
		return -1;
	}

	const bool allowSubset = opts.argExist("allowSubset");

	  // get reaction SMILES
	std::string SMILES = opts.getStrVal("reaction");
	  // setup ITS size range
	if (opts.argExist("k")) {
		minK = (size_t)opts.getIntVal("k");
		maxK = (size_t)opts.getIntVal("k");
	}
	if (opts.argExist("allITS")) {
		maxK = 999999;
	}
	  // container for all atom mapping statistics of interest
	AtomMappingStatistics statistics;

	// Reaction SMILES object...
	ReactionSMILES_Extractor chemReaction( SMILES, !opts.argExist("noProtonFilling") );

	  // create final educt/product graphs for CSP setup
	ReactionGraph educts( ggl::chem::Molecule_Graph_V( chemReaction.getEduMolecules() ), !allowSubset, false );
	ReactionGraph products( ggl::chem::Molecule_Graph_V( chemReaction.getProMolecules() ), !allowSubset, false );
//	ReactionGraph educts( ggl::chem::Molecule_Graph_V( chemReaction.getEduMolecules() ), !allowSubset, opts.argExist("v") );
//	ReactionGraph products( ggl::chem::Molecule_Graph_V( chemReaction.getProMolecules() ), !allowSubset, opts.argExist("v") );

	if (opts.argExist("ITS")) {
		  // setup ITS layout
		ITS curITS(opts.getStrVal("ITS"));
		  // check if ITS is compatible to reaction
		if (areCompatible( curITS, educts, products )) {
			  // setup CSP
			ITS_CSP::SymmetryHandler symHandler;
			ITS_CSP curCSP( curITS, educts, products, symHandler);
			  // find all mappings for the given ITS layout
			findMappings( curCSP
						, statistics
						, opts.argExist("printGraphs")
						, maxSolPerITS );
		}
	} else {

		 // find first ITS compatible with minK
		size_t curITSid = 0;
		assert( ! knownITSLayouts[curITSid].empty() );
		ITS curITS(knownITSLayouts[curITSid]);
		  // skip all ITS that are too small
		while( curITS.getSize() < minK  && !knownITSLayouts[curITSid+1].empty() ) {
			++curITSid;
			curITS = ITS(knownITSLayouts[curITSid]);
		}
		  // find all mappings for ITS layouts with minK <= size < maxK
		while( !knownITSLayouts[curITSid].empty() && curITS.getSize() <= maxK ) {

			if (opts.argExist("v")) {
				std::cout <<"testing ITS '" <<knownITSLayouts[curITSid] <<"'"<<std::endl;
			}

			  // check if ITS is compatible to reaction
			if (areCompatible( curITS, educts, products )) {
				  // setup CSP
				ITS_CSP::SymmetryHandler symHandler;
				ITS_CSP curCSP( curITS, educts, products, symHandler);
				  // find all mappings for the given ITS layout
				if(findMappings( curCSP
							, statistics
							, opts.argExist("printGraphs")
							, maxSolPerITS )
					&& curITS.getSize() < maxK )
				{
					  // if solution was found set maximum k to current size
					maxK = curITS.getSize();
				}
			}
			  // next ITS
			++curITSid;
			if (!knownITSLayouts[curITSid].empty()) {
				  // generate next ITS from string
				curITS = ITS(knownITSLayouts[curITSid]);
			}
		}

	}

	  // print statistics if interested
	if (opts.argExist("v")) {
		statistics.print();
	}

    return 0;
}// end main

