#ifndef ITS_CSP_H_
#define ITS_CSP_H_

#include <gecode/int.hh>
#include "ITS.h"
#include "ReactionGraph.h"
#include "ValenceMatrix.h"
#include "EdgeDegree.h"
#include "ChargeChange.h"
#include "AlternatingCycle.h"
#include "MinEdgeValence.h"
#include "ConnectedCoverage.h"
#include "ITSAtomLabel.h"
#include "NeighborhoodCount.h"

/*!
 * Defines a CSP for a cyclic imaginary transition state encoding.
 */
class ITS_CSP : public Gecode::Space  {


protected:

	  //! the targeted ITS encoding
	const ITS & its;

	  //! the given educt graphs
	const ReactionGraph & educts;

	  //! the given product graphs
	const ReactionGraph & products;

	  //! array containing atom variables of educts to be mapped
	Gecode::IntVarArray eduITS;

	  //! array containing variables of products i.e. mapping result
	Gecode::IntVarArray proITS;

	  //! array containing variables that hold what atom label is assigned to each ITS position
	Gecode::IntVarArray itsAtomLabel;

public:


	  /*!
	   * CSP setup based on ITS, educt, and product information.
	   */
	ITS_CSP(	const ITS & its
				, const ReactionGraph & educts
				, const ReactionGraph & products
			);


	virtual ~ITS_CSP();

	  /*!
	   * Copy constructor: returns a copy of the current space and updates
	   * CSP result during search i.e. space must be capable to return a
	   * copy of itself.
	   * @param bool share states whether a shared copy is constructed
	   * @param its is a reference to an object of the this class
	   */
	ITS_CSP( bool share, ITS_CSP & itsCSP);

	  /*!
	   * The function is responsible for space cloning. It calls the previously
	   * mentioned copy constructor.
	   * @param bool share states whether a shared copy is constructed
	   * @return a pointer to an updated copy of the underlying space
	   */
	virtual
	Space*
	copy( bool share );

	  /*!
	   * Access to the current educt ITS variable values.
	   * @return the assignment of the educt ITS variables
	   */
	sgm::Match
	getEductITS() const;

	  /*!
	   * Access to the current product ITS variable values.
	   * @return the assignment of the product ITS variables
	   */
	sgm::Match
	getProductITS() const;



	   /*!
	    * print the solution of current CSP
	    */
	void
	print(void) const;

	   /*!
	    * Access to reaction graph of educts in the underlying CSP
	    * @return educts graph
	    */
	const ReactionGraph &
	getEduGraph() const ;

	   /*!
	    * Access to reaction graph of products in the underlying CSP
	    * @return products graph
	    */
	const ReactionGraph &
	getProGraph() const ;

	 /*!
	  * Access to the ITS instance represented by this CSP
	  * @return the ITS encoding
	  */
	const ITS&
	getITS() const;

protected:




	  /*!
	   * Constraint post function to ensure the degree condition for the ITS
	   * @param home space
	   * @param x0 represents an educt variable
	   * @param x1 represents a product variable
	   * @param maxDegreeDiff the maximally allowed difference in degree
	   */
	void
	edgeDegree( Gecode::Home home
			, Gecode::IntVar x0
			, Gecode::IntVar x1
			, const size_t maxDegreeDiff );

	  /*!
	   * Constraint post function to ensure a change in atom charge
	   * @param home space
	   * @param x0 represents an educt variable
	   * @param x1 represents a product variable
	   * @param chargeChangeVal describes the change in charged atoms between educts and products
	   */
	void
	chargeChange( Gecode::Home home
			, Gecode::IntVar x0
			, Gecode::IntVar x1
			, const int chargeChangeVal);

	  /*!
	   * Constraint post function to ensure if the found ITS constitutes a simple
	   * alternating cycle structure.
	   * @param home space
	   * @param x0 represents an educt variable (edge from)
	   * @param x1 represents an educt variable (edge to)
	   * @param x2 represents a product variable (edge from)
	   * @param x3 represents a product variable (edge to)
	   * @param valence difference between ITS atom pairs (+1, -1, 0)
	   */
	void
	alternateCycle( Gecode::Home home
			, Gecode::IntVar x0
			, Gecode::IntVar x1
			, Gecode::IntVar x2
			, Gecode::IntVar x3
			, const int valenceDiff);

	  /*!
	   * Constraint  function to ensure a minimal bond valence between the
	   * two variables.
	   * @param x0 the source of the bond
	   * @param x0 the target of the bond
	   * @param minEdgeValence the minimally allowed bond valence
	   * @param valence the valence information to be used
	   */
	void
	minEdgeValence( Gecode::Home home
			, Gecode::IntVar x0, Gecode::IntVar x1
			, const int minEdgeValence
			, const ValenceMatrix & valences );

      /*!
       * Constraint post function to ensure atom label preservation
       * @param home space
       * @param x0 represents an educt variable
       * @param x1 represents a product variable
       * @param x2 stores the information if the variables are constrained to
       * a single label
       */
	void
	ensureLabelITS(Gecode::Home home, Gecode::IntVar x0, Gecode::IntVar x1, Gecode::IntVar x2);

	  /*!
	   * Constraint to the atoms which have to appear in the ITS ring, according
	   * to number of occurrences of their neighborhood patterns.
	   * @param home space
	   * @param array of atoms represents educts and products respectively
	   * @param neighID vector of neighborhoods identifiers
	   * @param neighCount maps between those IDs and their numbers of occurrences
	   */
	void
	countNeighborhood( Gecode::Home home
			, const Gecode::IntVarArray & atomsArray
			, const std::vector < size_t > & neighID
			, const std::map< std::size_t, std::size_t > & neighCount);

	  /*!
	   * Ensure the coverage of all connected components of a reaction
	   * SMIELS in the final mappings.
	   * @param home space
	   * @param x represents either an educt or a product array
	   * @param complabel labels of connected components
	   * @param connNum the number of the corresponding connected components
	   */
	void
	coverConnectedComps( Gecode::Home home
			, const Gecode::IntVarArray & x
			, const sgm::Graph_Interface::CompLabel & compLabel
			, const size_t connNum );

};

#endif /* ITS_CSP_H_ */
