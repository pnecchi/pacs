//----------------------------------------------------------------------
// Description: Class Edge that stores two unsigned int
// Author:      Pierpaolo Necchi
// Email:       pierpaolo.necchi@gmail.com
// Date:        ven 22 apr 2016 23:14:42 CEST
//----------------------------------------------------------------------

#ifndef EDGE_H
#define EDGE_H

#include <utility>    /* pair */
#include <algorithm>  /* min, max */

/** \class Edge 
 *  \bried class that implement an undirected edge inheriting from 
 *  std::pair<unsigned int, unsigned int>.
 *  
 *  The class forces the first element to be <= than the second one. In this
 *  way, exploiting the standard comparison operators for pairs, 
 *  (u, v) = (v, u) 
 */

class Edge : public  std::pair<unsigned int, unsigned int>
{
public:
	typedef std::pair<unsigned int, unsigned int> Base;

	// Constructors
	Edge() : Base() {}
	Edge(const unsigned int u_, const unsigned int v_)
		: Base(std::min(u_, v_), std::max(u_, v_)) {}
	Edge(const Edge& other) : Base(other) {};
	
	// Assignement opeator
	Edge& operator=(const Edge&)=default;

	// Destructor
	virtual ~Edge ()=default;
};

#endif /* end of include guard: EDGE_H */
