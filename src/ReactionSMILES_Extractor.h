/**
 *  Created on: 2012
 *  	@author Feras Nahar <naharf@informatik.uni-freiburg.de>
 */

#ifndef REACTIONSMILES_EXTRACTOR_H_
#define REACTIONSMILES_EXTRACTOR_H_

#include <ggl/chem/SMILESparser.hh>
#include <ggl/chem/Molecule.hh>

/**
 *  This class parses Reaction-SMILES format and split it to educt and product strings.
 *  Afterwards they are converted them to molecule graphs.
 */

class ReactionSMILES_Extractor {

public:

	  //! Vector of pointer on Molecules
	typedef std::vector< const ggl::chem::Molecule* > MolVec;

	  //! SMILES-reactions sides
	struct ReactionInputOutput {
		public:
			  //! concatenated educts, the string part before >>
			std::string inputSMILES;
			  //! concatenated products, the string part after >>
			std::string outputSMILES;
	};


protected:


	  /*!
	   * Vector of educts, products molecules constructed out
	   * of vector<string>
	   */
	MolVec eductsVec, productsVec;

public:

	  //! Construction
	  //! @param SMILESinput SMILES string to parse
	  //! @param fillProtons whether or not proton filling and molecule sanity
	  //! checks are to be done
	ReactionSMILES_Extractor(const std::string& SMILESinput
							    , const bool fillProtons = true );

	  //! Default Destructor
	~ReactionSMILES_Extractor();


	  //! Access to the molecule representation of the educts
	  //! @return vector of educts molecules
	const MolVec & getEduMolecules() const;

	  //! Access to the molecule representation of the products
	  //! @return vector of products molecules
	const MolVec & getProMolecules() const;


private:

	  //! String tokenizer customized according to SMILES-format
	  //! @param input string
	  //! @param tokenization delemiter
	  //! @return vector of tokens of string
	std::vector<std::string> tokenize(std::string & str, const char delemiter);

	  /*! Parses SMILES using ggl::chem::SMILESparser and saves the returned
	   *  values (Molecule) in eductsVec, productsVec respectively
	   *
	   *  @param educts the list of educt molecule SMILES
	   *  @param products the list of product molecule SMILES
	   *  @param fillProtons whether or not proton filling is to be done
	   */
	void SMILESparsing(const std::vector<std::string> & educts
					, const std::vector<std::string> & products
					, const bool fillProtons = true);

	  //! Removes all educt and product graph instances
	void clearMoleculeData();

};

#endif /* REACTIONSMILES_EXTRACTOR_H_ */
