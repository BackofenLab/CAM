/*
 * Created on: 2012
 *  	@author Feras Nahar <naharf@informatik.uni-freiburg.de>
 */

#ifndef CONNECTEDCOVERAGE_H_
#define CONNECTEDCOVERAGE_H_

#include <gecode/int.hh>
#include <sgm/Graph_Interface.hh>

/*!
 * A constraint ensuring that all connected components in the reaction SMILES
 * are present in the final mapping.
 */

// Define a  NaryPropagator since the constraints could be more that three.
typedef Gecode::NaryPropagator< Gecode::Int::IntView,
								Gecode::Int::PC_INT_DOM > NPropagator;
class ConnectedCoverage : public NPropagator
{

protected:

	  /*!
	   * x array of variables of type IntView represents either educt ITS
	   * candidates or a product ITS candidates.
	   */
	using NPropagator::x;

	  /*!
	   * Labels of the connected components of the underlying educts
	   * or products graph
	   */
	const sgm::Graph_Interface::CompLabel & complabel;

	  //! Number of the connected components
	const size_t connectedNum;

public:

	  /*!
	   * Construction
	   * @param home space of the current propagator
	   * @param _x represents an array of either educt atoms or product once
	   * @param complabel labels of connected components
	   * @param connNum the number of the corresponding connected components
	   */
	ConnectedCoverage( Gecode::Space & home
			   , Gecode::ViewArray< Gecode::Int::IntView > & _x
			   , const sgm::Graph_Interface::CompLabel & _complabel
			   , const size_t connNum );
	  /*!
	   * Propagator post function of the edge degree constraint
	   * @param home space of the current propagator
	   * @param _x represents an array of either educt atoms or product once
	   * @param complabel labels of connected components
	   * @param connNum the number of the corresponding connected components
	   */
	static
	Gecode::ExecStatus
	post( Gecode::Space & home
			, Gecode::ViewArray< Gecode::Int::IntView > & _x
			, const sgm::Graph_Interface::CompLabel & _complabel
			, const size_t connNum );

      /*!
	   * Copy constructor: returns a copy of the current space and updates CSP result during search
	   * i.e. space must be capable to return a copy of itself
	   * @paramspace home home of the edge degree propagator
	   * @param share bool share states whether a shared copy is constructed
	   * @param connCover is a reference to an object of the this propagator
       */
	ConnectedCoverage( Gecode::Space & home
						, bool share
						, ConnectedCoverage & connCover );

	  /*!
	   * The function is responsible for propagator cloning. It calls the previously mentioned copy ctor
	   * @param home home space of the current propagator
	   * @param bool share states whether a shared copy is constructed
	   * @return a pointer to an updated copy of the underlying propagator
	   */
	virtual
	Propagator*
	copy( Gecode::Space & home, bool share );

	  /*!
	   * Cost function of the current propagator. The returned cost controls
	   * when the propagator will be scheduled.
	   * @param space home of the edge degree propagator
	   * @param modification event describe how the domain of x0, x1 has changed
	   * @return cost value of the current propagator
	   */
	virtual
	Gecode::PropCost
	cost ( const Gecode::Space &, const Gecode::ModEventDelta & ) const;

	   /*!
	    * Executes the current propagator and reports about the propagation
	    * @param space home of the edge degree propagator.
	    * @param modification event describe how the domain of x0, x1 has changed.
	    * @return execution status of the propagator; OK posting was successful or FAILED
	    */
	virtual
	Gecode::ExecStatus
	propagate( Gecode::Space & home, const Gecode::ModEventDelta & );

	  //! Destruction
	~ConnectedCoverage();

};


#endif /* CONNECTEDCOVERAGE_H_ */
