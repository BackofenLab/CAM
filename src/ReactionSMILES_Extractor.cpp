
#include "ReactionSMILES_Extractor.h"
#include <stdexcept>

	 std::vector<std::string>
	 ReactionSMILES_Extractor::
	 tokenize( std::string & str, const char delemiter )
	 {
		 std::vector<std::string> tokens;
		 int i = 0, pos = 0;
		   // loop until no delemiter exists in the string
		 while ((pos = str.find(delemiter)) != std::string::npos) {
			   // extract substring before the given delemiter saving this token in the vector
			 tokens.push_back(str.substr(0, pos));
			   /*!
			    * get rid of already extracted token by assigning the substring located after
			    * the delemiter. It is going to be tokenized in the next iteration step
			    */
			 str = str.substr(pos + 1, str.size());
			 i++;
		 }
		 tokens.push_back(str);
		 return tokens;
	 }


////////////////////////////////////////////////////////////////////////////////////////////

	 ReactionSMILES_Extractor::
	 ReactionSMILES_Extractor( const std::string & SMILESinput
			 	 	 	 	 	 , const bool fillProtons )
	 {

		   // SPLIT SMILES INTO EDUCT INPUT AN PRODUCT OUTPUT

		  // concatenated educts, the string part before >>
		std::string inputSMILES;
		  // concatenated products, the string part after >>
		std::string outputSMILES;

		 size_t arrowPos1 = SMILESinput.find_first_of('>');
		 size_t arrowPos2 = SMILESinput.find_last_of('>');

		   /*
		    * check for a multiple SMILES reaction i.e. when '>>' appears more than once
		    * in the string. If any --> not allowed
		    */
		 if (arrowPos2 == (arrowPos1 +1)) {
			 if (arrowPos1 != std::string::npos) {
				   // extract the educts
				 inputSMILES = SMILESinput.substr(0, arrowPos1);
				   // arrowPos + 2 ignores '>>' and extract the products
				 outputSMILES = SMILESinput.substr(arrowPos1 + 2, arrowPos1 + SMILESinput.size());
			 }
		 } else {
			 throw std::runtime_error("#err: multiple reaction check: multiple reaction SMILES is not allowed!");
		 }

		   // SPLIT INTO INDIVIDUAL MOLECULE SMILES

		 std::vector<std::string> educts = tokenize(inputSMILES, '.');
		 std::vector<std::string> products = tokenize(outputSMILES, '.');

		   // PARSE SMILES

		 SMILESparsing( educts, products, fillProtons );

	 }

////////////////////////////////////////////////////////////////////////////////////////////

	 ReactionSMILES_Extractor::
	 ~ReactionSMILES_Extractor( )
	 {
		   // remove molecule data if any present
		 clearMoleculeData();
	 }

////////////////////////////////////////////////////////////////////////////////////////////

	 void
	 ReactionSMILES_Extractor::
	 clearMoleculeData( )
	 {
		 for (size_t i=0; i < eductsVec.size(); ++i) {
			 delete eductsVec[i];
		 }
		 eductsVec.clear();

		 for (size_t i=0; i < productsVec.size(); ++i) {
			 delete productsVec[i];
		 }
		 productsVec.clear();
	 }

////////////////////////////////////////////////////////////////////////////////////////////

	void
	ReactionSMILES_Extractor::
	SMILESparsing(const std::vector<std::string> & educts
					, const std::vector<std::string> & products
					, const bool fillProtons)
	 {
		std::pair< ggl::chem::Molecule, int > parseResult;
		for (std::vector<std::string>::const_iterator it = educts.begin(); it != educts.end(); it++) {
			parseResult = ggl::chem::SMILESparser::parseSMILES(*it);
			  // check parsing result: -1 successful
			if (parseResult.second == -1) {
				if (fillProtons) {
					  // enable proton filling
					ggl::chem::MoleculeUtil::fillProtons(parseResult.first);
				}
					  // Consistency check whether the entered educts are valid
// TODO REENABLE CHECKING
//					if (ggl::chem::MoleculeUtil::isConsistent(parseResult.first) ==
//							ggl::chem::MoleculeUtil::C_Consistent) {
						eductsVec.push_back(new ggl::chem::Molecule(parseResult.first));
//					} else {
//						throw std::runtime_error("#err: consistency check: educts or products molecule(s) are not consistent!");
//					}
			} else {
				throw std::runtime_error("#err: parsing check: unsuccessful parsing of educt or products molecule(s)!");
			}
		}

		for (std::vector<std::string>::const_iterator it = products.begin(); it != products.end(); it++) {
			parseResult = ggl::chem::SMILESparser::parseSMILES(*it);
			  // check parsing result: -1 successful
			if (parseResult.second == -1) {
				if (fillProtons) {
					  // enable proton filling
					ggl::chem::MoleculeUtil::fillProtons(parseResult.first);
				}
					  // Consistency check whether the entered products are valid
// TODO REENABLE CHECKING
//					if (ggl::chem::MoleculeUtil::isConsistent(parseResult.first) ==
//										ggl::chem::MoleculeUtil::C_Consistent) {
						  // adding protons to a product molecule
						productsVec.push_back(new ggl::chem::Molecule(parseResult.first));
//					} else {
//						throw std::runtime_error("#err: consistency check: educts or products molecule(s) are not consistent!");
//					}
			} else {
				throw std::runtime_error("#err: parsing check: unsuccessful parsing of educts or product molecule(s)!");
			}
		}

	 }

////////////////////////////////////////////////////////////////////////////////////////////

	 const ReactionSMILES_Extractor::MolVec &
	 ReactionSMILES_Extractor::
	 getEduMolecules() const
	 {
		return eductsVec;
	 }

////////////////////////////////////////////////////////////////////////////////////////////

	 const ReactionSMILES_Extractor::MolVec &
	 ReactionSMILES_Extractor::
	 getProMolecules() const
	 {
	 	return productsVec;
	 }

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
