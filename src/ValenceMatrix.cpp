
#include "ValenceMatrix.h"

#include <ggl/chem/MoleculeUtil.hh>

ValenceMatrix::ValenceMatrix()
	: SuperClass(0)
{
}

ValenceMatrix::ValenceMatrix( const sgm::Graph_Interface & graph )
	: SuperClass( graph.getNodeNumber(), graph.getNodeNumber(), 0)
{

	  // set valence and charge information
	sgm::Graph_Interface::OutEdge_iterator e;
	for (size_t from = 0; from < graph.getNodeNumber(); from++) {
		  // set charge information
		this->at( from, from) = ggl::chem::MoleculeUtil::getCharge(graph.getNodeLabel(from));
		  // set bond valence information for each bond
		for ( e = graph.getOutEdgesBegin(from); e != graph.getOutEdgesEnd(from); ++e) {
			/*!
			 * for each edge, transform its bond order to a number in the sparse matrix
			 * check for the existence of aromatic bonds
			 */
			const ggl::chem::MoleculeUtil::BondLabelData * const bond = ggl::chem::MoleculeUtil::getBondData((*e).getEdgeLabel());
			assert(bond != NULL);
			if ( bond->isAromatic) {
				std::cout << "#err: aromaticity check: aromatic bond are not supported!"
						  << std::endl;
			} else {
				  // set bond valence
				this->at( from, e->getToIndex() ) = bond->valence;
			}
		}
	}

}

//////////////////////////////////////////////////////////////////////////////

ValenceMatrix::~ValenceMatrix() {
}

//////////////////////////////////////////////////////////////////////////////

void
ValenceMatrix::
printMatrix( ) const
{
	for (size_t r = 0; r < this->rows; r++) {
		for (size_t c = 0; c < this->cols; c++)
			std::cout << this->at(r, c) << "  ";
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
