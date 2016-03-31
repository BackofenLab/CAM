#ifndef VALENCEMATRIX_H_
#define VALENCEMATRIX_H_


#include "MatrixSparse.hh"

#include <sgm/Graph_Interface.hh>


/*!
 *
 * Class holding the bond valence and atom charge information, e.g. from an
 * ReactionGraph instance retrieved by getGraphToMatch().
 *
 */
class ValenceMatrix : public MatrixSparseR< int >
{
public:

	typedef MatrixSparseR< int > SuperClass;

	ValenceMatrix();
	ValenceMatrix( const sgm::Graph_Interface & graph );
	virtual ~ValenceMatrix();
	void printMatrix( ) const ;
};

#endif /* VALENCEMATRIX_H_ */
