#ifndef ITS_H_
#define ITS_H_

#include <string>

#include <sgm/PA_OrderCheck.hh>

#include <ggl/Graph.hh>



/*!
 * Represents an imaginary transition state layout.
 */
class ITS {

public:

	  //! set of symmetric matches
	typedef std::set< sgm::Match > MatchSet;

protected:

	  //! graph representation of the ITS
	ggl::Graph itsGraph;

	  //! charge change for each atom
	std::vector< int > itsChargeChange;

	  //! bond valence change for each bond connecting atom i and i+1
	std::vector< int > itsBondChange;

	  //! list of order constraints to break all symmetries of the ITS
	sgm::PA_OrderCheck::CheckList itsOrderConstraints;

	  //! set of isomorphic matches of the ITS
	MatchSet itsSymmMatch;


public:

	  /*!
	   * Constructs an ITS from string encoding.
	   */
	ITS( const std::string & its );

	  //! destruction
	virtual ~ITS();

	  /*!
	   * Size of the ITS.
	   * @return the number of nodes of the ITS
	   */
	size_t
	getSize() const;

	  /*!
	   * Access to the ITS symmetries
	   * @return set of all symmetric matches of the ITS on itself
	   */
	const MatchSet &
	getSymmetricMatches() const ;

	  /*!
	   * Provides the list of required order checks to break all symmetries.
	   * @return list of needed order checks for symmetry breaking
	   */
	const sgm::PA_OrderCheck::CheckList &
	getNeededOrderChecks() const;

	  /*!
	   * Access to charge change of atoms in the underlying ITS graph
	   * @return corresponding charge change vector
	   */
	const std::vector< int > &
	getChargeChange() const ;

	  /*!
	   * Access to bond change of atoms in the underlying ITS graph
	   * @return corresponding bond change vector
	   */
	const std::vector< int > &
	getBondChange() const ;

	  /*!
	   * Access to the graph representation of the ITS.
	   * @return the graph representation of the ITS
	   */
	const ggl::Graph&
	getGraph() const;


protected:

	 /*!
	  * Constructs a graph representation of the given ITS string encoding
	  * @param its the ITS string encoding
	  * @return the graph encoding
	  */
	static
	ggl::Graph
	getGraph( const std::string & its );

};

#endif /* ITS_H_ */
